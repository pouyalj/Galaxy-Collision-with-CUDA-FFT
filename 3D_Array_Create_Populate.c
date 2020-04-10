// Generate and populate 3d array
// --------------------------------
// Useful link: https://www.techiedelight.com/dynamically-allocate-memory-for-3d-array/
// ---------------------------------------------------------------------------------------
// Compile with gcc Test.c -lm -o Tout.out -mcmodel=large
// ------------------------------------------------------------


#include <stdio.h>
#include <math.h>
#include<time.h>
#include<stdlib.h>

// M x N x O matrix
#define M 150
#define N 150
#define O 150


// Dynamically Allocate Memory for 3D Array
int main() {
    // dynamically allocate memory of size M*N*O
    float* threedArray = (float*) malloc(M * N * O* sizeof(float));

    // assign values to allocated memory
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < O; k++) {
				*(threedArray + i*N*O + j*O + k) = rand() % 100;
            }
        }
    }
    // print the 3D array
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < O; k++) {
				printf("%f\n", *(threedArray + i*N*O + j*O + k));
            }
        }
    }
}