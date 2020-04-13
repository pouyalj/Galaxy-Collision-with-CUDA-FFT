//Compile with: nvcc CUDAfft2.0.cu -I/home/phyd57/N_Body1/9.2/include -L/home/phyd57/N_Body1/9.2/lib64 -lcufft -o CUDAfftcu2.out

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cufft.h>

#include <stdlib.h>
#include<stdio.h>
#include <iostream>

//#include "dislin.h"

#define N 256 // N is the sidelength of the image -> N^3 pixels in entire image
#define block_size_x 2 
#define block_size_y 2
#define block_size_z 2

float den_array[N][N][N];
float grav_po[N][N][N];
float image[N][N];

__global__ void real2complex(cufftComplex *c, float *a, int n);
__global__ void complex2real_scaled(float *a, cufftComplex *c, float scale, int n);
__global__ void solve_poisson(cufftComplex *c, float *kx, float *ky, float *kz, int n);


void FFT_poisson(float den_array[N][N][N], float grav_po[N][N][N], float L)
{
	int x, y, z, i;

	float *kx, *ky, *kz, *den;
	kx = (float *)malloc(sizeof(float)*N);
	ky = (float *)malloc(sizeof(float)*N);
	kz = (float *)malloc(sizeof(float)*N);
	den = (float *)malloc(sizeof(float)*N*N*N);

	float *kx_d, *ky_d, *kz_d, *den_d;
	cufftComplex *den_complex_d;
	cudaMalloc((void **)&kx_d, sizeof(float) * N);
	cudaMalloc((void **)&ky_d, sizeof(float) * N);
	cudaMalloc((void **)&kz_d, sizeof(float) * N);
	cudaMalloc((void **)&den_d, sizeof(float) * N * N * N);
	cudaMalloc((void **)&den_complex_d, sizeof(cufftComplex) * N * N * N);

	printf("A bunch of Mallocs\n");

	#pragma omp for
	for (x = 0; x < N; x++)
		for (y = 0; y < N; y++)
			for (z = 0; z < N; z++)
				den[x + y*N + z*N*N] = den_array[x][y][z];

	printf("Den transfer\n");

	float* den_inital = (float *)malloc(sizeof(float) * N * N * N);
	for (i = 0; i < N * N; i++)
		den_inital[i] = den[i];

	for (i = 0; i < N; i++)
	{
		kx[i] = (i / N)* L * 2.0f * 3.14159f;
		ky[i] = (i / N)* L * 2.0f * 3.14159f;
		kz[i] = (i / N)* L * 2.0f * 3.14159f;
	}

	printf("ks made, Memcpy next\n");

	cudaMemcpy(kx_d, kx, sizeof(float)*N, cudaMemcpyHostToDevice);
	cudaMemcpy(ky_d, ky, sizeof(float)*N, cudaMemcpyHostToDevice);
	cudaMemcpy(kz_d, kz, sizeof(float)*N, cudaMemcpyHostToDevice);
	cudaMemcpy(den_d, den, sizeof(float)*N*N*N, cudaMemcpyHostToDevice);
	printf("Memcpy completed\n");

	cufftHandle plan;
	cufftPlan3d(&plan,N,N,N,CUFFT_C2C);
	printf("Plan created\n");

	/* Compute the execution configuration, block_size_x*block_size_y*block_size_z = number of threads */
	dim3 dimBlock(block_size_x, block_size_y, block_size_z);
	dim3 dimGrid(N/dimBlock.x, N/dimBlock.y, N/dimBlock.z);
	/* Handle N not multiple of block_size_x, block_size_y, or block_size_y */
	if (N % block_size_x != 0) dimGrid.x += 1;
	if (N % block_size_y != 0) dimGrid.y += 1;
	if (N % block_size_z != 0) dimGrid.z += 1;

	real2complex<<<dimGrid, dimBlock>>>(den_complex_d, den_d, N);

	cufftExecC2C(plan, den_complex_d, den_complex_d, CUFFT_FORWARD);

	solve_poisson<<<dimGrid, dimBlock>>>(den_complex_d, kx_d, ky_d, kz_d, N);

	cufftExecC2C(plan, den_complex_d, den_complex_d, CUFFT_INVERSE);

	float scale = 1.0f / (N*N*N);
	complex2real_scaled<<<dimGrid, dimBlock>>>(den_d, den_complex_d, scale, N);
	

	cudaMemcpy(den, den_d, sizeof(float)*N*N*N, cudaMemcpyDeviceToHost);

	#pragma omp for
	for (x = 0; x < N; x++)
		for (y = 0; y < N; y++)
			for (z = 0; z < N; z++)
				grav_po[x][y][z] = den[x + y*N + z*N*N];

	printf("Transfered back, next dystroy\n");

	/* Destroy plan and clean up memory on device*/
	cudaFree(kx);
	cudaFree(ky);
	cudaFree(kz);
	cudaFree(den);
	cudaFree(den_inital);
	cufftDestroy(plan);
	cudaFree(den_complex_d);
	cudaFree(den);
	cudaFree(kx_d);
	cudaFree(ky_d);
	cudaFree(kz_d);

	printf("Destruction complete\n");
}

/*
void make_image(float array[N][N][N])
{
	int x, y, z;
	float Max = -500.0, Min = 500.0
	
	#pragma omp for
	for (x = 0; x < N; x++)
		for (y = 0; y < N; y++)
			image[x][y] = 0.0;

	#pragma omp for
	for (x = 0; x < N; x++)
		for (y = 0; y < N; y++)
			for (z = 0; z < N; z++)
				image[x][y] += den_array[x][y][z];

	#pragma omp for
	for (x = 0; x < N; x++)
		for (y = 0; y < N; y++)
			if (image[x][y] > Max) Max = image[x][y];
			if (image[x][y] < Min) Min = image[x][y]

	metafl("cons");
	disini();
	pagera();
	hwfont();

	titlin("Potential map", 4);
	//titlin("anthing below", 2)

	name("X [kP]", "x");
	name("Y [kP]", "y");
	name("Potential in Z", "z");

	intax()	;
	autres(800,800);
	axspos(300,1850);
	ax3len(1400,1400,1400);
	
	call labdig(4, "YZ")
	call labdig(1, "x")
	graf3(-N/2, N/2, -N/2, N/10, -N/2, N/2, -N/2, N/10, Min, Max, Min, (Max - Min)/10.0);
	crvmat((float *)XY, 800, 800, 1, 1);

	height(50);
	title();
	disfin();
}
*/

int main()
{
	int i, j, k;

	#pragma omp parallel for
	for (i = 0; i < 256; i ++)
	{
		for (j = 0; j < 256; j++)
		{
			for (k = 0; k < 256; k++)
			{
				den_array[i][j][k] = 0.0;
				grav_po[i][j][k] = 0.0;
			}
		}
	}

	den_array[128][128][128] = 500.0;

	#pragma omp parallel for
	for (i = 143; i > 113; i --)
	{
		for (j = 113; j < 143; j++)
		{
			printf("%f,", den_array[j][i][128]);
		}
		
		printf("\n");

	}

	FFT_poisson(den_array, grav_po, 128);

	printf("z = 127:\n");
	#pragma omp parallel for
	for (i = 143; i > 113; i --)
	{
		for (j = 113; j < 143; j++)
		{
			printf("%f,", grav_po[j][i][127]);
		}
		
		printf("\n");

	}
	
	printf("\n\n");
	printf("z = 128:\n");
	#pragma omp parallel for
	for (i = 143; i > 113; i --)
	{
		for (j = 113; j < 143; j++)
		{
			printf("%f,", grav_po[j][i][128]);
		}
		
		printf("\n");

	}

	printf("\n\n");
	printf("z = 129:\n");
	#pragma omp parallel for
	for (j = 143; j > 113; j --)
	{
		for (i = 113; i < 143; i++)
		{
			printf("%.1f,", grav_po[i][j][129]);
		}
		
		printf("\n");

	}

	printf("\n%f\n", grav_po[0][128][128]);
	printf("%f\n", grav_po[255][128][128]);
	printf("%f\n", grav_po[128][0][128]);
	printf("%f\n", grav_po[128][255][128]);

	//make_image(grav_po);
}

__global__ void real2complex(cufftComplex *c, float *a, int n)
{
    /* compute idx, idy, and idz, the location of the element in the original NxNxN array */
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;
	int idz = blockIdx.z * blockDim.z + threadIdx.z;
	if (idx < n && idy < n && idz < n)
	{
		int index = idx + idy*n + idz*n*n;
		c[index].x = a[index];
		c[index].y = 0.0f;
	}
}

__global__ void complex2real_scaled(float *a, cufftComplex *c, float scale, int n)
{
	/* compute idx and idy, the location of the element in the original NxN array */
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;
	int idz = blockIdx.z * blockDim.z + threadIdx.z;
	if (idx < n && idy < n && idz < n)
	{
		int index = idx + idy*n + idz*n*n;
		a[index] = scale * c[index].x;
	}
}


__global__ void solve_poisson(cufftComplex *c, float *kx, float *ky, float *kz, int n)
{
	/* compute idx and idy, the location of the element in the original NxN array */
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;
	int idz = blockIdx.z * blockDim.z + threadIdx.z;
	if (idx < n && idy < n && idz < n)
	{
		int index = idx + idy*n + idz*n*n;
		float scale = -(kx[idx]*kx[idx] + ky[idy]*ky[idy] + kz[idz]*kz[idz]) + 0.00001f;
		if (idx == 0 && idy == 0 && idz == 0) scale = 1.0f;
		scale = 1.0f / scale;
		c[index].x *= scale;
		c[index].y *= scale;
	}
}
