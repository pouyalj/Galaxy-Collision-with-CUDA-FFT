//Compile with: nvcc CUDAfft2.0.cu -I/home/phyd57/N_Body1/9.2/include -L/home/phyd57/N_Body1/9.2/lib64 -lcufft -o CUDAfftcu2.out -I/usr/local/dislin -ldislin
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cufft.h>

#include <stdlib.h>
#include<stdio.h>
#include <iostream>
#include<time.h>

#include "dislin.h"

#define N 256 // N is the sidelength of the image -> N^3 pixels in entire image
#define M 100000000 //M is the number of particles.
#define block_size_x 2 
#define block_size_y 2
#define block_size_z 2

float den_array[N][N][N];
float grav_po[N][N][N];
//float particleArray[M][7];
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
				den[x + y*N + z*N*N] = 4.0 * 3.14159 * 1.4006 * den_array[x][y][z];
				//Where 1.4006 is G in units kPc**3/solar_mass * 10kyears

	float* den_inital = (float *)malloc(sizeof(float) * N * N * N);

	#pragma omp for
	for (i = 0; i < N * N; i++)
		den_inital[i] = den[i];

	#pragma omp for
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


void make_image(float array[N][N][N], const char *output_name, const char *title)
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

	titlin(title, 4);
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
/*
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
}*/

void densArray(float **particleArray, float*** threedArray) {
	int i, j, k = 0;
    // dynamically allocate memory of size M*N*O
	// assign values to allocated memory
	/*
	for (i = 0; i < I; i++) {
		threedArray[i] = (float**)malloc(J * sizeof(float*));
        if (threedArray[i] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
        for (j = 0; j < J; j++) {
			threedArray[i][j] = (float*)malloc(K * sizeof(float));
            if (threedArray[i][j] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
	}*/

    // assign values to allocated memory
	#pragma omp for
	for (i=0; i < M; i++) {
        // printf("%d\n", i);
        // printf("%d\n", (int)floorf(particleArray[i][0]));
        // printf("%d %d %d \n", (int)floorf(particleArray[i][0]), (int)floorf(particleArray[i][1]), (int)floorf(particleArray[i][2]));
        threedArray[(int)floorf(particleArray[i][0])][(int)floorf(particleArray[i][1])][(int)floorf(particleArray[i][2])] =
        threedArray[(int)floorf(particleArray[i][0])][(int)floorf(particleArray[i][1])][(int)floorf(particleArray[i][2])] + 1;
	}
    // // print the 3D array
	// for (i = 0; i < I; i++)
	// {
	// 	for (j = 0; j < J; j++)
	// 	{
	// 		for (k = 0; k < K; k++)
	// 			printf("%f\n", threedArray[i][j][k]);
	//    	}
	// }
}

void center_diff(int xN, int yN, int zN, float*** grav_po, float **particleArray) {
    int i, j, k, l;
    float v_half, x, v;

    // float gx[I][J][K], gy[I][J][K], gz[I][J][K];
    // float (*g)[I][J][K];

    // for(i=1; i<xN; i++){
    //     for(j=1; j<yN-1; j++){
    //         for(k=1; k<zN-1; k++){
    //             gx[i][j][k] = (grav_po[i+1][j][k] - grav_po[i-1][j][k])/(2); // get g for each directions
    //             gy[i][j][k] = (grav_po[i][j+1][k] - grav_po[i][j-1][k])/(2);
    //             gz[i][j][k] = (grav_po[i][j][k+1] - grav_po[i][j][k-1])/(2);
    //         }
    //     }
    // }
    // printf("g force created\n");

	#pragma omp for
    for(i=0; i<M; i++){
        for(l=0; l<1; l++){
            v_half = particleArray[i][l+3] + 
            (grav_po[(int)round(particleArray[i][0])+1][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])]
                - grav_po[(int)round(particleArray[i][0])-1][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])])/(2);
            x = particleArray[i][l] + v_half;
            v = v_half + 
            (grav_po[(int)round(particleArray[i][0])+1][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])]
                - grav_po[(int)round(particleArray[i][0])-1][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])])/(2);
            particleArray[i][l+3] = v;
        }
        for(l=1; l<2; l++){
            v_half = particleArray[i][l+3] + 
            (grav_po[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])+1][(int)round(particleArray[i][2])] 
            - grav_po[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])-1][(int)round(particleArray[i][2])])/(2);
            x = particleArray[i][l] + v_half;
            v = v_half + 
            (grav_po[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])+1][(int)round(particleArray[i][2])] 
            - grav_po[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])-1][(int)round(particleArray[i][2])])/(2);
            particleArray[i][l+3] = v;
        }
        for(l=2; l<3; l++){
            v_half = particleArray[i][l+3] + 
            (grav_po[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])+1] 
            - grav_po[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])-1])/(2);
            x = particleArray[i][l] + v_half;
            v = v_half + 
            (grav_po[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])+1] 
            - grav_po[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])-1])/(2);
            particleArray[i][l+3] = v;
        }
            // move all particles
            // updater(particleArray[i][l+3], particleArray[i][l],
            // *g[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])]);
            // v_half = particleArray[i][l+3] +
            // *g[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])]/2+
            // 1190/(pow(particleArray[i][0],2)+pow(particleArray[i][1],2)+pow(particleArray[i][2],2));
            // x = particleArray[i][l] + v_half;
            // v = v_half +
            // *g[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])]/2;
            // particleArray[i][l] = x;
            // particleArray[i][l+3] = v;
    }

    // // update density array (TDB)
    // printf("density array updater initiated\n");
}

int main()
{
	//initialize particle array without velocity.
	int i, j, k, index, max_number, min_number, counter;
	float t, dt, X, Y, R, V;
	float *particleArray = (float *)malloc(M * sizeof(float *));
	
	t = 0.0;
	dt = 1.0;

	#pragma omp for
	for (i = 0; i < M; i++) {
		particleArray[i] = (float *)malloc(M * sizeof(float*));

		if (particleArray[i] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
	}

	// first galaxy population
	#pragma omp for
    for (i = 0; i < (int)(M*0.05/2); i++) {
        particleArray[i][0] = 2*1.41*cos((float)(rand()%629)/100) + 96.0;
        particleArray[i][1] = 2*1.41*sin((float)(rand()%629)/100) + 96.0;
        particleArray[i][2] = (float)(rand()%(50+1))/1000 + 128.0; // rand() % (max_number + 1 - minimum_number) + minimum_number
        X = particleArray[i][0] - 96;
        Y = particleArray[i][1] - 96;
        R = sqrt(pow(X,2) + pow(Y,2));
        V = sqrt(1190.0*R);
        particleArray[i][3] = Y/R*V + 0.04;
        particleArray[i][4] = X/R*V + 0.04;
        particleArray[i][5] = 0;
    }
	
	#pragma omp for
    for (index=1; index<11; index++){
        for (i = (int)(M*0.05/2+((index-1)*0.095*M/2)); i < (int)(M*0.05/2+((index)*0.095*M/2)); i++) {
            particleArray[i][0] = (2+index)*1.41*cos((float)(rand()%629)/100) + 96.0;
            particleArray[i][1] = (2+index)*1.41*sin((float)(rand()%629)/100) + 96.0;
            particleArray[i][2] = (float)(rand()%(50+1))/1000 + 128.0;
            X = particleArray[i][0] - 96;
            Y = particleArray[i][1] - 96;
            R = sqrt(pow(X,2) + pow(Y,2));
            V = sqrt(1190.0*R);
            particleArray[i][3] = Y/R*V + 0.04;
            particleArray[i][4] = X/R*V + 0.04;
			particleArray[i][5] = 0;
        }
    }
	
	#pragma omp for
	for (i = 0; i < (int)(M/2); i++) {
		for (j=6;j<7;j++){
			particleArray[i][j] = 0.0; // 0.0 is indicator for Milky Way
		}
	}

	// second galaxy population
	#pragma omp for
	for (i = (int)(M*0.05/2+((10)*0.095*M/2)); i < (int)(M*0.05/2+((10)*0.095*M/2))+(int)(M*0.05/2); i++) {
        particleArray[i][0] = 2*1.41*cos((float)(rand()%629)/100)  + 160.0;
        particleArray[i][1] = 2*1.41*sin((float)(rand()%629)/100)  + 160.0;
        particleArray[i][2] = (float)(rand()%(50+1))/1000 + 128.0;
        X = particleArray[i][0] - 96;
        Y = particleArray[i][1] - 96;
        R = sqrt(pow(X,2) + pow(Y,2));
        V = sqrt(1190.0*R);
        particleArray[i][3] = Y/R*V - 0.04;
        particleArray[i][4] = X/R*V - 0.04;
        particleArray[i][5] = 0;
    }

	#pragma omp for
	for (index=11; index<21; index++){
        for (i = (int)(M*0.05+((index-1)*0.095*M/2)); i < (int)(M*0.05+((index)*0.095*M/2)); i++) {
            particleArray[i][0] = (2+index-10)*1.41*cos((float)(rand()%629)/100)  + 160.0;
            particleArray[i][1] = (2+index-10)*1.41*sin((float)(rand()%629)/100)  + 160.0;
            particleArray[i][2] = (float)(rand()%(150+1))/1000 + 128.0;
            X = particleArray[i][0] - 96;
            Y = particleArray[i][1] - 96;
            R = sqrt(pow(X,2) + pow(Y,2));
            V = sqrt(1190*R);
            particleArray[i][3] = Y/R*V - 0.04;
            particleArray[i][4] = X/R*V - 0.04;
            particleArray[i][5] = 0;
        }
    }

	#pragma omp for
	for (i = (int)(M/2); i < M; i++) {
		particleArray[i][6] = 1.0; // 1.0 is indicator for Andromeda
	}
	
	//create initial velocity, for each array.
	
	///Repeat until finished.
	while (t < 500)
	{
		densArray(particleArray, den_array);
		FFT_poisson(den_array, grav_po);
		enter_diff(256, 256, 256, grav_po, particleArray);
		
		if (time == 0.0)
		{
			make_image(den_array, "Initial.png", "Initial density of the system");
		}
		
		if (time == 125.0)
		{
			make_image(den_array, "fourth.png", "Density of the system after 1,250,000 years");
		}
		
		if (time == 250.0)
		{
			make_image(den_array, "half.png", "Density of the system after 2,500,000 years");
		}
		
		if (time == 375.0)
		{
			make_image(den_array, "three_fourths.png", "Density of the system after 3,750,000 years");
		}
		
		time += dt;
	}
	
	//Fill density array with both galaxies
	//Find potential
	//update particle with potential
	
	//end.
	
	make_image(den_array, "final.png", "Density of the system after 5,000,000 years");
	
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