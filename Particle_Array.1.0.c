// Generate and populate 3d array
// --------------------------------
// Useful link: https://www.techiedelight.com/dynamically-allocate-memory-for-3d-array/
// ---------------------------------------------------------------------------------------
// Compile with gcc Particle_Array.1.0.c -lm -o Particle_Array.1.0.out -mcmodel=large
// --------------------------------------------------------------------------------------------------------
// Can be used given a design decision to be made by the group
// --------------------------------------------------------------



#include <stdio.h>
#include <math.h>
#include<time.h>
#include<stdlib.h>

// M x N x O matrix
#define M 100000000
#define N 8
// #define O 150


// Dynamically Allocate Memory for 3D Array
int main() {
	int i, j, k;
    // dynamically allocate memory of size M*N*O
    float **particleArray = (float **)malloc(M * sizeof(float *));
    

    // assign values to allocated memory
	for (i = 0; i < M; i++) {
		particleArray[i] = (float *)malloc(N * sizeof(float*));

		if (particleArray[i] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
	}

    // assign values to allocated memory
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			particleArray[i][j] = rand() % 100;
        }
    }


    // print the 2D array
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++) {
            printf("%f\n", particleArray[i][j]);

	   	}
		printf("\n");
	}

	// deallocate memory
	for (i = 0; i < M; i++) {
        free(particleArray[i]);
	}
    free(particleArray);
}