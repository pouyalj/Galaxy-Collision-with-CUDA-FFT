//Compile with: nvcc CUDAfft.c -I/home/phyd57/N_Body1/9.2/include -L/home/phyd57/N_Body1/9.2/lib64 -lcufft -o CUDAfftc.out

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cufft.h>

#include <stdlib.h>
#include <iostream>

#define N 256 // N is the sidelength of the image -> N^3 pixels in entire image
#define block_size_x 2 
#define block_size_y 2
#define block_size_z 2

__global__ void real2complex(cufftComplex *c, float *a, int n);
__global__ void complex2real_scaled(float *a, cufftComplex *c, float scale, int n);
__global__ void solve_poisson(cufftComplex *c, float *kx, float *ky, int n);


void FFT_poisson(void) //will have some real values latter
{
	float *kx, *ky, *kz, *den;
	kx = (float *)malloc(sizeof(float) * N);
	ky = (float *)malloc(sizeof(float) * N);
	ky = (float *)malloc(sizeof(float) * N);
	r = (float *)malloc(sizeof(float) * N * N);

	float *kx_d, *ky_d, *r_d;
	cufftComplex *r_complex_d;
	cudaMalloc((void **)&kx_d, sizeof(float) * N);
	cudaMalloc((void **)&ky_d, sizeof(float) * N);
	cudaMalloc((void **)&r_d, sizeof(float) * N * N);
	cudaMalloc((void **)&r_complex_d, sizeof(cufftComplex) * N * N);

	for (int y = 0; y < N; y++)
		for (int x = 0; x < N; x++)
			r[x + y * N] = rand() / (float)RAND_MAX;

	float* r_inital = (float *)malloc(sizeof(float) * N * N);
	for (int i = 0; i < N * N; i++)
		r_inital[i] = r[i];

	for (int i = 0; i < N; i++)
	{
		kx[i] = i - N / 2.0f; //centers kx values to be at center of image
		ky[i] = N / 2.0f - i; //centers ky values to be at center of image
	}

	cudaMemcpy(kx_d, kx, sizeof(float) * N, cudaMemcpyHostToDevice);
	cudaMemcpy(ky_d, ky, sizeof(float) * N, cudaMemcpyHostToDevice);
	cudaMemcpy(r_d, r, sizeof(float) * N * N, cudaMemcpyHostToDevice);

	cufftHandle plan;
	cufftPlan2d(&plan, N, N, CUFFT_C2C);

	/* Compute the execution configuration, block_size_x*block_size_y = number of threads */
	dim3 dimBlock(block_size_x, block_size_y);
	dim3 dimGrid(N / dimBlock.x, N / dimBlock.y);
	/* Handle N not multiple of block_size_x or block_size_y */
	if (N % block_size_x != 0) dimGrid.x += 1;
	if (N % block_size_y != 0) dimGrid.y += 1;

	real2complex << < dimGrid, dimBlock >> > (r_complex_d, r_d, N);

	cufftExecC2C(plan, r_complex_d, r_complex_d, CUFFT_FORWARD);
	solve_poisson << <dimGrid, dimBlock >> > (r_complex_d, kx_d, ky_d, N);
	cufftExecC2C(plan, r_complex_d, r_complex_d, CUFFT_INVERSE);

	float scale = 1.0f / (N * N);
	complex2real_scaled << <dimGrid, dimBlock >> > (r_d, r_complex_d, scale, N);

	cudaMemcpy(r, r_d, sizeof(float) * N * N, cudaMemcpyDeviceToHost);

	for (int i = 0; i < N * N; i++)
		std::cout << i << "\tr_initial: " << r_inital[i] << "\tr: " << r[i] << std::endl;
	system("pause");

	/* Destroy plan and clean up memory on device*/
	free(kx);
	free(ky);
	free(kz);
	free(r);
	free(r_inital);
	cufftDestroy(plan);
	cudaFree(r_complex_d);
	cudaFree(kx_d);
}

__global__ void real2complex(cufftComplex *c, float *a, int n)
{
    /* compute idx and idy, the location of the element in the original NxN array */
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;
	if (idx < n && idy < n)
	{
		int index = idx + idy * n;
		c[index].x = a[index];
		c[index].y = 0.0f;
	}
}

__global__ void complex2real_scaled(float *a, cufftComplex *c, float scale, int n)
{
	/* compute idx and idy, the location of the element in the original NxN array */
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;
	if (idx < n && idy < n)
	{
		int index = idx + idy * n;
		a[index] = scale * c[index].x;
	}
}


__global__ void solve_poisson(cufftComplex *c, float *kx, float *ky, int n)
{
	/* compute idx and idy, the location of the element in the original NxN array */
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;
	if (idx < n && idy < n)
	{
		int index = idx + idy * n;
		float scale = -(kx[idx] * kx[idx] + ky[idy] * ky[idy]) + 0.00001f;
		if (idx == n/2 && idy == n/2) scale = 1.0f;
		scale = 1.0f / scale;
		c[index].x *= scale;
		c[index].y *= scale;
	}
}
