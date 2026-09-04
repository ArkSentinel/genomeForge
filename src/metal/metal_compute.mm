#include "metal/metal_compute.hpp"

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include <chrono>
#include <iostream>

namespace genome::metal {

struct MetalComputeEngine::Impl {
    id<MTLDevice> device{nil};
    id<MTLCommandQueue> commandQueue{nil};
    id<MTLComputePipelineState> pipelineState{nil};
    id<MTLComputePipelineState> physicsPipelineState{nil};
};

MetalComputeEngine::MetalComputeEngine()
    : pimpl_(new Impl()) {
    initialize();
}

MetalComputeEngine::~MetalComputeEngine() {
    if (pimpl_) {
        delete pimpl_;
        pimpl_ = nullptr;
    }
}

bool MetalComputeEngine::initialize() {
    @autoreleasepool {
        pimpl_->device = MTLCreateSystemDefaultDevice();
        if (!pimpl_->device) {
            initialized_ = false;
            return false;
        }

        device_info_.device_name = [[pimpl_->device name] UTF8String];
        device_info_.has_unified_memory = [pimpl_->device hasUnifiedMemory];
        device_info_.max_buffer_size = [pimpl_->device maxBufferLength];

        pimpl_->commandQueue = [pimpl_->device newCommandQueue];
        if (!pimpl_->commandQueue) {
            initialized_ = false;
            return false;
        }

        // Metal Shading Language (MSL) kernels para MoE GEMM y Física Celular
        NSString* shaderSource = @R"(
            #include <metal_stdlib>
            using namespace metal;

            struct GPUCellParticle {
                float x;
                float y;
                float vx;
                float vy;
                float radius;
                float pad[3];
            };

            kernel void moe_gemm_kernel(
                device const float* X [[buffer(0)]],
                device const float* W [[buffer(1)]],
                device float* Out [[buffer(2)]],
                constant uint& N [[buffer(3)]],
                constant uint& K [[buffer(4)]],
                constant uint& M [[buffer(5)]],
                uint2 gid [[thread_position_in_grid]]
            ) {
                uint row = gid.y;
                uint col = gid.x;

                if (row >= N || col >= M) return;

                float sum = 0.0f;
                uint x_offset = row * K;
                uint w_offset = col * K;

                for (uint i = 0; i < K; ++i) {
                    sum += X[x_offset + i] * W[w_offset + i];
                }

                Out[row * M + col] = sum;
            }

            kernel void cell_physics_kernel(
                device GPUCellParticle* particles [[buffer(0)]],
                constant uint& count [[buffer(1)]],
                constant float& dt [[buffer(2)]],
                constant float& dish_radius [[buffer(3)]],
                uint id [[thread_position_in_grid]]
            ) {
                if (id >= count) return;

                GPUCellParticle self = particles[id];
                float fx = 0.0f;
                float fy = 0.0f;

                for (uint j = 0; j < count; ++j) {
                    if (id == j) continue;
                    float dx = self.x - particles[j].x;
                    float dy = self.y - particles[j].y;
                    float dist_sq = dx * dx + dy * dy;
                    float min_dist = self.radius + particles[j].radius;

                    if (dist_sq < min_dist * min_dist && dist_sq > 0.0001f) {
                        float dist = sqrt(dist_sq);
                        float overlap = min_dist - dist;
                        float force = overlap * 12.0f;
                        fx += (dx / dist) * force;
                        fy += (dy / dist) * force;
                    }
                }

                // Amortiguamiento de viscosidad en medio líquido
                fx -= self.vx * 2.0f;
                fy -= self.vy * 2.0f;

                self.vx += fx * dt;
                self.vy += fy * dt;
                self.x += self.vx * dt;
                self.y += self.vy * dt;

                // Límite circular de la placa de Petri
                float r_sq = self.x * self.x + self.y * self.y;
                float max_r = dish_radius - self.radius;
                if (r_sq > max_r * max_r) {
                    float r = sqrt(r_sq);
                    self.x = (self.x / r) * max_r;
                    self.y = (self.y / r) * max_r;
                    self.vx = -self.vx * 0.4f;
                    self.vy = -self.vy * 0.4f;
                }

                particles[id] = self;
            }
        )";

        NSError* error = nil;
        id<MTLLibrary> library = [pimpl_->device newLibraryWithSource:shaderSource options:nil error:&error];
        if (!library) {
            if (error) {
                std::cerr << "Metal shader compilation error: " << [[error localizedDescription] UTF8String] << std::endl;
            }
            initialized_ = false;
            return false;
        }

        id<MTLFunction> moeFunction = [library newFunctionWithName:@"moe_gemm_kernel"];
        if (!moeFunction) {
            initialized_ = false;
            return false;
        }

        pimpl_->pipelineState = [pimpl_->device newComputePipelineStateWithFunction:moeFunction error:&error];
        if (!pimpl_->pipelineState) {
            initialized_ = false;
            return false;
        }

        id<MTLFunction> physicsFunction = [library newFunctionWithName:@"cell_physics_kernel"];
        if (physicsFunction) {
            pimpl_->physicsPipelineState = [pimpl_->device newComputePipelineStateWithFunction:physicsFunction error:&error];
        }

        initialized_ = true;
        return true;
    }
}

bool MetalComputeEngine::run_batch_gemm(
    const float* X, int N, int K,
    const float* W, int M,
    float* Out,
    float* elapsed_us
) {
    if (!initialized_ || !pimpl_->device || !pimpl_->commandQueue || !pimpl_->pipelineState) {
        return false;
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    @autoreleasepool {
        size_t size_X = N * K * sizeof(float);
        size_t size_W = M * K * sizeof(float);
        size_t size_Out = N * M * sizeof(float);

        // Apple Silicon Unified Memory: MTLResourceStorageModeShared permite compartir
        // el búfer entre CPU y GPU sin copias de memoria intermedias (Zero-Copy)
        id<MTLBuffer> buf_X = [pimpl_->device newBufferWithBytes:X length:size_X options:MTLResourceStorageModeShared];
        id<MTLBuffer> buf_W = [pimpl_->device newBufferWithBytes:W length:size_W options:MTLResourceStorageModeShared];
        id<MTLBuffer> buf_Out = [pimpl_->device newBufferWithLength:size_Out options:MTLResourceStorageModeShared];

        uint32_t n_val = static_cast<uint32_t>(N);
        uint32_t k_val = static_cast<uint32_t>(K);
        uint32_t m_val = static_cast<uint32_t>(M);

        id<MTLCommandBuffer> commandBuffer = [pimpl_->commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];

        [encoder setComputePipelineState:pimpl_->pipelineState];
        [encoder setBuffer:buf_X offset:0 atIndex:0];
        [encoder setBuffer:buf_W offset:0 atIndex:1];
        [encoder setBuffer:buf_Out offset:0 atIndex:2];
        [encoder setBytes:&n_val length:sizeof(n_val) atIndex:3];
        [encoder setBytes:&k_val length:sizeof(k_val) atIndex:4];
        [encoder setBytes:&m_val length:sizeof(m_val) atIndex:5];

        MTLSize gridSize = MTLSizeMake(M, N, 1);
        NSUInteger threadGroupWidth = std::min((NSUInteger)M, pimpl_->pipelineState.maxTotalThreadsPerThreadgroup);
        NSUInteger threadGroupHeight = std::min((NSUInteger)N, pimpl_->pipelineState.maxTotalThreadsPerThreadgroup / threadGroupWidth);
        if (threadGroupHeight == 0) threadGroupHeight = 1;
        MTLSize threadGroupSize = MTLSizeMake(threadGroupWidth, threadGroupHeight, 1);

        [encoder dispatchThreads:gridSize threadsPerThreadgroup:threadGroupSize];
        [encoder endEncoding];

        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];

        // Leer resultado desde la memoria compartida (Zero-Copy UMA)
        float* result_ptr = static_cast<float*>([buf_Out contents]);
        std::memcpy(Out, result_ptr, size_Out);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    if (elapsed_us) {
        *elapsed_us = std::chrono::duration<float, std::micro>(end_time - start_time).count();
    }

    return true;
}

bool MetalComputeEngine::run_cell_physics(
    GPUCellParticle* particles, uint32_t count, float dt, float dish_radius, float* elapsed_us
) {
    if (!initialized_ || !pimpl_->device || !pimpl_->commandQueue || !pimpl_->physicsPipelineState || count == 0) {
        return false;
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    @autoreleasepool {
        size_t buffer_size = count * sizeof(GPUCellParticle);

        // Memoria unificada compartida Zero-Copy (UMA)
        id<MTLBuffer> particle_buffer = [pimpl_->device newBufferWithBytes:particles length:buffer_size options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> commandBuffer = [pimpl_->commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];

        [encoder setComputePipelineState:pimpl_->physicsPipelineState];
        [encoder setBuffer:particle_buffer offset:0 atIndex:0];
        [encoder setBytes:&count length:sizeof(count) atIndex:1];
        [encoder setBytes:&dt length:sizeof(dt) atIndex:2];
        [encoder setBytes:&dish_radius length:sizeof(dish_radius) atIndex:3];

        MTLSize gridSize = MTLSizeMake(count, 1, 1);
        NSUInteger threadGroupWidth = std::min((NSUInteger)count, pimpl_->physicsPipelineState.maxTotalThreadsPerThreadgroup);
        if (threadGroupWidth == 0) threadGroupWidth = 1;
        MTLSize threadGroupSize = MTLSizeMake(threadGroupWidth, 1, 1);

        [encoder dispatchThreads:gridSize threadsPerThreadgroup:threadGroupSize];
        [encoder endEncoding];

        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];

        GPUCellParticle* result_ptr = static_cast<GPUCellParticle*>([particle_buffer contents]);
        std::memcpy(particles, result_ptr, buffer_size);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    if (elapsed_us) {
        *elapsed_us = std::chrono::duration<float, std::micro>(end_time - start_time).count();
    }

    return true;
}

} // namespace genome::metal
