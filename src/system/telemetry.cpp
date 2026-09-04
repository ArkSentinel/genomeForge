#include "system/telemetry.hpp"

#include <mach/mach_host.h>
#include <mach/mach_init.h>
#include <mach/processor_info.h>
#include <mach/vm_map.h>
#include <mach/vm_statistics.h>
#include <sys/sysctl.h>
#include <sys/types.h>
#include <unistd.h>
#include <cstring>
#include <algorithm>
#include <numeric>

namespace genome::system {

TelemetryCollector::TelemetryCollector() {
    init_platform_info();
    update(); // Primera muestra para inicializar ticks
}

TelemetryCollector::~TelemetryCollector() = default;

void TelemetryCollector::init_platform_info() {
    // 1. Modelo de Chip (ej. Apple M1/M2/M3/M4 Max/Pro)
    char buffer[256];
    size_t size = sizeof(buffer);
    if (sysctlbyname("machdep.cpu.brand_string", buffer, &size, nullptr, 0) == 0) {
        current_metrics_.chip_model = std::string(buffer);
    } else {
        current_metrics_.chip_model = "Apple Silicon (ARM64)";
    }

    // 2. Conteo de núcleos P (Performance) y E (Efficiency)
    int p_cores = 0;
    size = sizeof(p_cores);
    if (sysctlbyname("hw.perflevel0.logicalcpu", &p_cores, &size, nullptr, 0) == 0) {
        current_metrics_.p_core_count = p_cores;
    } else {
        current_metrics_.p_core_count = 8; // fallback razonable
    }

    int e_cores = 0;
    size = sizeof(e_cores);
    if (sysctlbyname("hw.perflevel1.logicalcpu", &e_cores, &size, nullptr, 0) == 0) {
        current_metrics_.e_core_count = e_cores;
    } else {
        current_metrics_.e_core_count = 4;
    }

    // 3. Conteo total de núcleos
    int total_cores = 0;
    size = sizeof(total_cores);
    if (sysctlbyname("hw.ncpu", &total_cores, &size, nullptr, 0) == 0) {
        num_cpu_cores_ = static_cast<uint32_t>(total_cores);
    } else {
        num_cpu_cores_ = current_metrics_.p_core_count + current_metrics_.e_core_count;
    }

    // 4. Memoria RAM total (Unified Memory)
    uint64_t memsize = 0;
    size = sizeof(memsize);
    if (sysctlbyname("hw.memsize", &memsize, &size, nullptr, 0) == 0) {
        current_metrics_.total_memory_bytes = memsize;
    }

    // 5. Inicializar lista de métricas por núcleo
    current_metrics_.core_metrics.resize(num_cpu_cores_);
    for (uint32_t i = 0; i < num_cpu_cores_; ++i) {
        current_metrics_.core_metrics[i].core_id = static_cast<int>(i);
        // En Apple Silicon, los P-Cores suelen estar en perflevel0
        current_metrics_.core_metrics[i].is_performance_core = 
            (static_cast<int>(i) < current_metrics_.p_core_count);
        current_metrics_.core_metrics[i].usage_percent = 0.0f;
    }
}

void TelemetryCollector::sample_cpu() {
    natural_t cpu_count = 0;
    processor_info_array_t cpu_info = nullptr;
    mach_msg_type_number_t num_cpu_info = 0;

    kern_return_t kr = host_processor_info(
        mach_host_self(),
        PROCESSOR_CPU_LOAD_INFO,
        &cpu_count,
        &cpu_info,
        &num_cpu_info
    );

    if (kr != KERN_SUCCESS || cpu_info == nullptr) {
        return;
    }

    auto* load_info = reinterpret_cast<processor_cpu_load_info_t>(cpu_info);
    uint32_t active_cores = std::min(static_cast<uint32_t>(cpu_count), num_cpu_cores_);

    if (!has_prev_ticks_) {
        prev_cpu_ticks_.resize(active_cores * CPU_STATE_MAX);
        for (uint32_t i = 0; i < active_cores; ++i) {
            for (int state = 0; state < CPU_STATE_MAX; ++state) {
                prev_cpu_ticks_[i * CPU_STATE_MAX + state] = load_info[i].cpu_ticks[state];
            }
        }
        has_prev_ticks_ = true;
        vm_deallocate(mach_task_self(), reinterpret_cast<vm_address_t>(cpu_info), num_cpu_info * sizeof(int));
        return;
    }

    float total_usage_sum = 0.0f;
    float p_core_sum = 0.0f;
    float e_core_sum = 0.0f;
    int p_sampled = 0;
    int e_sampled = 0;

    for (uint32_t i = 0; i < active_cores; ++i) {
        uint64_t delta_user = load_info[i].cpu_ticks[CPU_STATE_USER] - prev_cpu_ticks_[i * CPU_STATE_MAX + CPU_STATE_USER];
        uint64_t delta_sys  = load_info[i].cpu_ticks[CPU_STATE_SYSTEM] - prev_cpu_ticks_[i * CPU_STATE_MAX + CPU_STATE_SYSTEM];
        uint64_t delta_idle = load_info[i].cpu_ticks[CPU_STATE_IDLE] - prev_cpu_ticks_[i * CPU_STATE_MAX + CPU_STATE_IDLE];
        uint64_t delta_nice = load_info[i].cpu_ticks[CPU_STATE_NICE] - prev_cpu_ticks_[i * CPU_STATE_MAX + CPU_STATE_NICE];

        uint64_t total_delta = delta_user + delta_sys + delta_idle + delta_nice;
        float usage = 0.0f;
        if (total_delta > 0) {
            usage = static_cast<float>(delta_user + delta_sys + delta_nice) / static_cast<float>(total_delta) * 100.0f;
        }

        usage = std::clamp(usage, 0.0f, 100.0f);
        current_metrics_.core_metrics[i].usage_percent = usage;
        total_usage_sum += usage;

        if (current_metrics_.core_metrics[i].is_performance_core) {
            p_core_sum += usage;
            p_sampled++;
        } else {
            e_core_sum += usage;
            e_sampled++;
        }

        // Actualizar ticks previos
        for (int state = 0; state < CPU_STATE_MAX; ++state) {
            prev_cpu_ticks_[i * CPU_STATE_MAX + state] = load_info[i].cpu_ticks[state];
        }
    }

    if (active_cores > 0) {
        current_metrics_.total_cpu_usage = total_usage_sum / static_cast<float>(active_cores);
    }
    if (p_sampled > 0) {
        current_metrics_.p_core_avg_usage = p_core_sum / static_cast<float>(p_sampled);
    }
    if (e_sampled > 0) {
        current_metrics_.e_core_avg_usage = e_core_sum / static_cast<float>(e_sampled);
    }

    vm_deallocate(mach_task_self(), reinterpret_cast<vm_address_t>(cpu_info), num_cpu_info * sizeof(int));
}

void TelemetryCollector::sample_memory() {
    vm_size_t page_size = 0;
    mach_port_t mach_port = mach_host_self();
    host_page_size(mach_port, &page_size);

    vm_statistics64_data_t vm_stat;
    mach_msg_type_number_t host_size = sizeof(vm_statistics64_data_t) / sizeof(integer_t);

    kern_return_t kr = host_statistics64(
        mach_port,
        HOST_VM_INFO64,
        reinterpret_cast<host_info64_t>(&vm_stat),
        &host_size
    );

    if (kr == KERN_SUCCESS) {
        uint64_t free_mem = static_cast<uint64_t>(vm_stat.free_count) * page_size;
        uint64_t active_mem = static_cast<uint64_t>(vm_stat.active_count) * page_size;
        uint64_t inactive_mem = static_cast<uint64_t>(vm_stat.inactive_count) * page_size;
        uint64_t wired_mem = static_cast<uint64_t>(vm_stat.wire_count) * page_size;
        uint64_t compressed_mem = static_cast<uint64_t>(vm_stat.compressor_page_count) * page_size;

        uint64_t used_mem = active_mem + wired_mem + compressed_mem;
        current_metrics_.free_memory_bytes = free_mem + inactive_mem;
        current_metrics_.used_memory_bytes = used_mem;
        current_metrics_.wired_memory_bytes = wired_mem;

        if (current_metrics_.total_memory_bytes > 0) {
            current_metrics_.memory_used_percent = 
                (static_cast<float>(used_mem) / static_cast<float>(current_metrics_.total_memory_bytes)) * 100.0f;
        }
    }
}

void TelemetryCollector::update() {
    sample_cpu();
    sample_memory();
}

} // namespace genome::system
