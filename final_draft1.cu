//Compile with: nvcc CUDAfft2.0.cu -I/home/phyd57/N_Body1/9.2/include -L/home/phyd57/N_Body1/9.2/lib64 -lcufft -o CUDAfftcu2.out -I/usr/local/dislin -ldislin
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cufft.h>

#include <stdlib.h>
#include<stdio.h>
#include <iostream>

#include "dislin.h"

#define N 256 // N is the sidelength of the image -> N^3 pixels in entire image
#define M 100000000 //M is the number of particles.
#define block_size_x 2 
#define block_size_y 2
#define block_size_z 2

float den_array[N][N][N];
float grav_po[N][N][N];
float part_array[M][6];
float image[N/2][N/2];

__global__ void real2complex(cufftComplex *c, float *a, int n);
__global__ void complex2real_scaled(float *a, cufftComplex *c, float scale, int n);
__global__ void solve_poisson(cufftComplex *c, float *k_xyz, int n);


void FFT_poisson(float den_array[N][N][N], float grav_po[N][N][N])
{
	int x, y, z, i;

	float *k_xyz, *den;
	k_xyz = (float *)malloc(sizeof(float)*N);
	den = (float *)malloc(sizeof(float)*N*N*N);

	float *k_xyz_d, *den_d;
	cufftComplex *den_complex_d;
	cudaMalloc((void **)&k_xyz_d, sizeof(float) * N);
	cudaMalloc((void **)&den_d, sizeof(float) * N * N * N);
	cudaMalloc((void **)&den_complex_d, sizeof(cufftComplex) * N * N * N);

	#pragma omp for
	for (x = 0; x < N; x++)
		for (y = 0; y < N; y++)
			for (z = 0; z < N; z++)
				den[x + y*N + z*N*N] = den_array[x][y][z];

	float* den_inital = (float *)malloc(sizeof(float) * N * N * N);
	for (i = 0; i < N * N; i++)
		den_inital[i] = den[i];

	for (i = 0; i < N; i++)
	{
		if (i < N/2)
		{
			k_xyz[i] = i;
		}
		else
		{
			k_xyz[i] = i-N;
		}
	}

	cudaMemcpy(k_xyz_d, k_xyz, sizeof(float)*N, cudaMemcpyHostToDevice);
	cudaMemcpy(den_d, den, sizeof(float)*N*N*N, cudaMemcpyHostToDevice);

	cufftHandle plan;
	cufftPlan3d(&plan,N,N,N,CUFFT_C2C);

	/* Compute the execution configuration, block_size_x*block_size_y*block_size_z = number of threads */
	dim3 dimBlock(block_size_x, block_size_y, block_size_z);
	dim3 dimGrid(N/dimBlock.x, N/dimBlock.y, N/dimBlock.z);
	/* Handle N not multiple of block_size_x, block_size_y, or block_size_y */
	if (N % block_size_x != 0) dimGrid.x += 1;
	if (N % block_size_y != 0) dimGrid.y += 1;
	if (N % block_size_z != 0) dimGrid.z += 1;

	real2complex<<<dimGrid, dimBlock>>>(den_complex_d, den_d, N);

	cufftExecC2C(plan, den_complex_d, den_complex_d, CUFFT_FORWARD);

	solve_poisson<<<dimGrid, dimBlock>>>(den_complex_d, k_xyz_d, N);

	cufftExecC2C(plan, den_complex_d, den_complex_d, CUFFT_INVERSE);

	float scale = 1.0f / (N*N*N);
	complex2real_scaled<<<dimGrid, dimBlock>>>(den_d, den_complex_d, scale, N);
	

	cudaMemcpy(den, den_d, sizeof(float)*N*N*N, cudaMemcpyDeviceToHost);

	#pragma omp for
	for (x = 0; x < N; x++)
		for (y = 0; y < N; y++)
			for (z = 0; z < N; z++)
				grav_po[x][y][z] = den[x + y*N + z*N*N];

	/* Destroy plan and clean up memory on device*/
	cudaFree(k_xyz);
	cudaFree(den);
	cudaFree(den_inital);
	cufftDestroy(plan);
	cudaFree(den_complex_d);
	cudaFree(den);
	cudaFree(k_xyz_d);
}


void make_image(float array[N][N][N], const char *output_name)
{
	int x, y, z;
	float Max = -500.0, Min = 500.0;
	
	#pragma omp for
	for (x = 0; x < N/2; x++)
		for (y = 0; y < N/2; y++)
			image[x][y] = 0.0;

	#pragma omp for
	for (x = 0; x < N/2; x++)
		for (y = 0; y < N/2; y++)
			for (z = 0; z < N/2; z++)
				image[x][y] += array[x+N/4][y+N/4][z+N/4];

	#pragma omp for
	for (x = 0; x < N/2; x++)
	{
		for (y = 0; y < N/2; y++)
		{
			if (image[x][y] > Max)
			{
				Max = image[x][y];
			}
			if (image[x][y] < Min)
			{
				Min = image[x][y];
			}
		}
	}

	metafl("PNG");
	setfil(output_name);
	//metafl("CONS");
	disini();
	pagera();
	hwfont();

	titlin("Potential map", 4);
	//titlin("anthing below", 2)

	name("X [kP]", "x");
	name("Y [kP]", "y");
	name("Potential in Z", "z");

	intax()	;
	autres(N/2,N/2);
	axspos(300,1850);
	ax3len(1600,1600,1600);
	
	labdig(6, "Z");
	graf3(-N/4, N/4, -N/4, N/40, -N/4, N/4, -N/4, N/40, Min, Max, Min, (Max-Min)/10);
	crvmat((float *)image, N/2, N/2 , 1, 1);

	height(50);
	title();
	disfin();
}

void *CM_finder(int galaxy_ID, float xyz_array[3])
{
	/*
	Fills xyz_array with the z, y, and z values of the CM of
	a given galaxy, in that order.
	galaxy_ID is 1 for galaxy 1 and 2 for galaxy 2.
	*/
	int i, n;
	
	if (galaxy_ID == 1) n = 0;
	else n = M/2;
	
	#pragma omp for
	for (i = 0; i < M/2; i ++)
	{
		xyz_array[0] = part_array[i+n][0];
		xyz_array[1] = part_array[i+n][1];
		xyz_array[2] = part_array[i+n][2];
	}
	
	xyz_array[0] /= (float)M/2;
	xyz_array[1] /= (float)M/2;
	xyz_array[2] /= (float)M/2;
}

void initial_velocity(int galaxy_ID)
{
	float CM_array[3];
	CM_finder(galaxy_ID, xyz_array[3]);
	
	int i, n;
	float x, y, z, r, v;
	
	if (galaxy_ID == 1) n = 0;
	else n = M/2;
	
	#pragma omp for
	for (i = 0; i < M/2; i ++)
	{
		x = xyz_array[0] - part_array[i+n][0];
		y = xyz_array[1] - part_array[i+n][1];
		z = xyz_array[2] - part_array[i+n][2];
		r = x*x + y*y + z*z;
		r = pow(r, 0.5)
		
		v = //pow(G*m*M/r, 0.5); need the unit of time to know the value of G
		
		//from there I need the direction it moves from there.
	}
	
	pow(value, 0.5);
	
	//Also should add the 402000 km/h here
}

int main()
{
	int i, j, k;

	//initialize particle array without velocity.
	
	//make density array of G1.
	//call to find potential.
	//use potential to find initial velocity.
	//update particle array.
	
	//make density array of G2
	//call to find potential.
	//use potential to find initial velocity.
	//update particle array.
	
	///Repeat until finished.
	//Fill density array with both galaxies
	//Find potential
	//update particle with potential
	
	//end.
	
	return 0;

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


__global__ void solve_poisson(cufftComplex *c, float *k_xyz, int n)
{
	/* compute idx and idy, the location of the element in the original NxN array */
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;
	int idz = blockIdx.z * blockDim.z + threadIdx.z;
	if (idx < n && idy < n && idz < n)
	{
		int index = idx + idy*n + idz*n*n;
		float scale = -(k_xyz[idx]*k_xyz[idx] + k_xyz[idy]*k_xyz[idy] + k_xyz[idz]*k_xyz[idz]) + 0.00001f;
		if (idx == 0 && idy == 0 && idz == 0) scale = 1.0f;
		scale = 1.0f / scale;
		c[index].x *= scale;
		c[index].y *= scale;
	}
}