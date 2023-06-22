#include <cstdlib>
#include <cstdio>

#include <stdio.h>
#include "Interval.h"
#include "Functions.h"
#include "IntervalVector.h"
#include <iostream>
#include <vector>
#include "AutoDiff/IntervalAutoDiff.h"
#include "AutoDiff/AutoDiff.h"
#include "Matrix.h"


// Kernel function for vector sum
__global__ void vectorSumKernel(const float* input1, const float* input2, float* output, size_t size) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        output[idx] = input1[idx] + input2[idx];
    }
}

// Host function for vector sum using CUDA
void vectorSumCUDA(const std::vector<float>& vec1, const std::vector<float>& vec2, std::vector<float>& result) {
    INTERVAL EXP(1,1);
    size_t size = vec1.size();
    float* deviceInput1;
    float* deviceInput2;
    float* deviceOutput;

    // Allocate memory on the GPU
    cudaMalloc(&deviceInput1, size * sizeof(float));
    cudaMalloc(&deviceInput2, size * sizeof(float));
    cudaMalloc(&deviceOutput, size * sizeof(float));

    // Copy input vectors from host to device
    cudaMemcpy(deviceInput1, vec1.data(), size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceInput2, vec2.data(), size * sizeof(float), cudaMemcpyHostToDevice);

    // Launch the kernel
    int threadsPerBlock = 256;
    int numBlocks = (size + threadsPerBlock - 1) / threadsPerBlock;
    vectorSumKernel<<<numBlocks, threadsPerBlock>>>(deviceInput1, deviceInput2, deviceOutput, size);

    // Copy result from device to host
    cudaMemcpy(result.data(), deviceOutput, size * sizeof(float), cudaMemcpyDeviceToHost);

    // Free memory on the GPU
    cudaFree(deviceInput1);
    cudaFree(deviceInput2);
    cudaFree(deviceOutput);
}

int main() {
    std::vector<float> vec1 = {1.0f, 2.0f, 3.0f, 4.0f};
    std::vector<float> vec2 = {0.5f, 1.5f, 2.5f, 3.5f};
    std::vector<float> result(vec1.size());

    vectorSumCUDA(vec1, vec2, result);

    // Print the result
    std::cout << "Result: ";
    for (float value : result) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    return 0;
}