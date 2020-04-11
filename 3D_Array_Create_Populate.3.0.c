// Generate and populate 3d array
// --------------------------------
// Useful link: https://www.techiedelight.com/dynamically-allocate-memory-for-3d-array/
// ---------------------------------------------------------------------------------------
// Compile with gcc 3D_Array_Create_Populate.3.0.c -lm -o 3D_Array_Create_Populate.3.0.out -mcmodel=large
// --------------------------------------------------------------------------------------------------------
// Can be used given a design decision to be made by the group
// --------------------------------------------------------------


#include <stdio.h>
#include <math.h>
#include<time.h>
#include<stdlib.h>

// M x N x O matrix
#define M 100
#define N 100
#define O 10


// Dynamically Allocate Memory for 3D Array
int main() {
	int i, j, k;
    // dynamically allocate memory of size M*N*O
    float threedArray[M][N][O];

    // assign values to allocated memory
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < O; k++) {
				threedArray[i][j][k] = rand() % 100;
            }
		}
	}
    // print the 3D array
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			for (k = 0; k < O; k++)
				printf("%f\n", threedArray[i][j][k]);
	   	}
	}
}