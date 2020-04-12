// Generate and populate 3d array
// --------------------------------
// Useful link: https://www.techiedelight.com/dynamically-allocate-memory-for-3d-array/
// ---------------------------------------------------------------------------------------
// Compile with gcc Numerical_Differentiator.1.0.c -lm -o Numerical_Differentiator.1.0.out -mcmodel=large
// ------------------------------------------------------------


#include <stdio.h>
#include <math.h>
#include<time.h>
#include<stdlib.h>
#include <assert.h>




// updater takes the initial values and moves everything one step forward in time.
void updater(float v0, float x0, float a){
    float v_half, x, v;
    v_half = v0 + a/2;
    x = x0 + v_half;
    v = v_half + a/2;
    x0 = x;
    v0 = v;
}

void center_diff(int xN, int yN, int zN, float grav_po[xN][yN][zN],
 float den_array[xN][yN][zN], float particle_array[xN][8]) {
    int i, j, k, l;
    float gx[64][64][8], gy[64][64][8], gz[64][64][8];
    float (*g)[64][64][8];
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
            updater(particle_array[i][l+3], particle_array[i][l],
            *g[(int)round(particle_array[i][0])][(int)round(particle_array[i][1])][(int)round(particle_array[i][2])]);
        }
    }

    // update density array (TDB)
    printf("density array updater initiated\n");

    printf("all arrays updated\n");
}




void main() {
    int xN, yN = 64;
    int zN = 8;
    float grav_po[64][64][8];
    float den_array[64][64][8];
    float particle_array[64][8];
    center_diff(64, 64, 8, grav_po, den_array, particle_array);
    for(int i=1; i<64; i++){
        for(int j=1; j<8; j++){
            printf("%f ", particle_array[i][j]);
        }
        printf("\n");
    }               
}