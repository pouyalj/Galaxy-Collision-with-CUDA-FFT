// Generate and populate 3d array
// --------------------------------
// Useful link: https://www.techiedelight.com/dynamically-allocate-memory-for-3d-array/
// ---------------------------------------------------------------------------------------
// Compile with gcc densA_particleA.2.0.c -lm -o densA_particleA.2.0.out -mcmodel=large
// --------------------------------------------------------------------------------------------------------
// Can be used given a design decision to be made by the group
// --------------------------------------------------------------


#include <stdio.h>
#include <math.h>
#include<time.h>
#include<stdlib.h>

// M x N x O matrix
// #define M 64
// // #define N 64
#define M 100000000
#define N 8

#define I 256
#define J 256
#define K 256


// Dynamically Allocate Memory for 3D Array
void densArray(float **particleArray, float*** threedArray) {
	int i, j, k = 0;
    // dynamically allocate memory of size M*N*O
	// assign values to allocated memory
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
	}
	
	}
    printf("density array intitiation complete\n");

    // assign values to allocated memory
	for (i=0; i < M; i++) {
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
    printf("Density Array completed\n");
}


// Dynamically Allocate Memory for 3D Array
void create_particle_Array(float **particleArray) {
	int i, j, k, index, max_number, min_number;
    // dynamically allocate memory of size M*N*O
    

    // assign values to allocated memory
	for (i = 0; i < M; i++) {
		particleArray[i] = (float *)malloc(N * sizeof(float*));

		if (particleArray[i] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
	}

    printf("Starting to populate the particle array\n");
    // assign values to allocated memory
    for (i = 0; i < (int)(M*0.05); i++) {
        for (j = 0; j < 2; j++) {
            particleArray[i][j] = (float)(rand()%1001)/1000; // rand() % (max_number + 1 - minimum_number) + minimum_number
        }
		for (j=2;j<3;j++){
			particleArray[i][j] = (float)(rand()%(50+1))/1000;
		}
    }
    for (index=1; index<13; index++){
        for (i = (int)(0.05*M+((index-1)*0.095)); i < (int)(0.05*M+((index)*0.095)); i++) {
            max_number = 1000*(index+1);
            min_number = 1000*index;
            for (j = 0; j < 2; j++) {
                particleArray[i][j] = ((float)(rand()%(max_number+1-min_number))+(float)min_number)/1000;
            }
			for (j=2;j<3;j++){
				particleArray[i][j] = (float)(rand()%(150+1))/1000;
			}
        }
    }


    // print the 2D array
	for (i = 0; i < 10; i++)
	{
		for (j = 0; j < 3; j++) {
            printf("%f\n", particleArray[i][j]);

	   	}
		printf("\n");
	}

	// // deallocate memory
	// for (i = 0; i < M; i++) {
    //     free(particleArray[i]);
	// }
    // free(particleArray);
}



void main(){
    int i;
    float **particleArray = (float **)malloc(M * sizeof(float *));
    float*** threedArray = (float***)malloc(I * sizeof(float**));
    create_particle_Array(particleArray);
    printf("Particle Array completed\n");
    densArray(particleArray, threedArray);
    // deallocate memory
	for (i = 0; i < M; i++) {
        free(particleArray[i]);
	}
    free(particleArray);
}