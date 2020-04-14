// Generate and populate 3d array
// --------------------------------
// Useful link: https://www.techiedelight.com/dynamically-allocate-memory-for-3d-array/
// ---------------------------------------------------------------------------------------
// Compile with gcc Numerical_Differentiator.2.0.c -lm -o Numerical_Differentiator.2.0.out -mcmodel=large
// ------------------------------------------------------------


#include <stdio.h>
#include <math.h>
#include<time.h>
#include<stdlib.h>
#include <assert.h>


#define M 100000000
#define N 8

#define I 256
#define J 256
#define K 256


// updater takes the initial values and moves everything one step forward in time.
void updater(float v0, float x0, float a){
    float v_half, x, v;
    v_half = v0 + a/2;
    x = x0 + v_half;
    v = v_half + a/2;
    x0 = x;
    v0 = v;
}

void center_diff(int xN, int yN, int zN, float grav_po[I][J][K],
 float den_array[I][J][K], float particleArray[M][N]) {
    int i, j, k, l;
    float v_half, x, v;

    float gx[I][J][K], gy[I][J][K], gz[I][J][K];
    float (*g)[I][J][K];
    for(i=1; i<xN; i++){
        for(j=1; j<yN-1; j++){
            for(k=1; k<zN-1; k++){
                gx[i][j][k] = (grav_po[i+1][j][k] - grav_po[i-1][j][k])/(2); // get g for each directions
                gy[i][j][k] = (grav_po[i][j+1][k] - grav_po[i][j-1][k])/(2);
                gz[i][j][k] = (grav_po[i][j][k+1] - grav_po[i][j][k-1])/(2);
            }
        }
    }
    printf("g force created\n");

    printf("updater function initiated\n");
    for(i=0; i<xN; i++){
        for(l=0; l<3; l++){
            if(l==0){
                g = &gx;
            }
            else if(l==1){
                g = &gy;
            }
            else{
                g = &gz;
            }
            // move all particles
            // updater(particleArray[i][l+3], particleArray[i][l],
            // *g[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])]);
            v_half = particleArray[i][l+3] +
            *g[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])]/2;
            x = particleArray[i][l] + v_half;
            v = v_half +
            *g[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])]/2;
            particleArray[i][l] = x;
            particleArray[i][l+3] = v;

        }
    }

    // update density array (TDB)
    printf("density array updater initiated\n");

    printf("all arrays updated\n");
}


// Dynamically Allocate Memory for 3D Array
void densArray(float **particleArray, float threedArray[I][J][K]) {
	int i, j, k = 0;
    // dynamically allocate memory of size M*N*O
    printf("density array intitiation complete\n");

    // assign values to allocated memory
	for (i; i < M; i++) {
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




void main() {
    float grav_po[I][J][K];
    float threedArray[I][J][K];
    float particleArray[M][8];
    center_diff(I, J, K, grav_po, threedArray, particleArray);
    densArray(particleArray, threedArray);

    for(int i=1; i<64; i++){
        for(int j=1; j<8; j++){
            printf("%f ", particleArray[i][j]);
        }
        printf("\n");
    }               
}