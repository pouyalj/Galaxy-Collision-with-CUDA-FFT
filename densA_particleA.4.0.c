// Generate and populate 3d array
// --------------------------------
// Useful link: https://www.techiedelight.com/dynamically-allocate-memory-for-3d-array/
// ---------------------------------------------------------------------------------------
// Compile with gcc densA_particleA.4.0.c -lm -o densA_particleA.4.0.out -mcmodel=large
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
#define N 7

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
    printf("Density Array completed\n");
}


// Dynamically Allocate Memory for 3D Array
void create_particle_Array(float **particleArray) {
	int i, j, k, index;
    float X, Y, R, V;
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

	// first galaxy population
    for (i = 0; i < (int)(M*0.05/2); i++) {
        particleArray[i][0] = 2*1.41*cos((float)(rand()%629)/100) + 96.0;
        particleArray[i][1] = 2*1.41*sin((float)(rand()%629)/100) + 96.0;
        particleArray[i][2] = (float)(rand()%(50+1))/1000 + 128.0; // rand() % (max_number + 1 - minimum_number) + minimum_number
        X = particleArray[i][0] - 96;
        Y = particleArray[i][1] - 96;
        R = sqrt(pow(X,2) + pow(Y,2));
        V = sqrt(1190*R);
        particleArray[i][3] = Y/R*V;
        particleArray[i][4] = X/R*V;
        particleArray[i][5] = 0;
    }
    for (index=1; index<11; index++){
        for (i = (int)(M*0.05/2+((index-1)*0.095*M/2)); i < (int)(M*0.05/2+((index)*0.095*M/2)); i++) {
            particleArray[i][0] = (2+index)*1.41*cos((float)(rand()%629)/100) + 96.0;
            particleArray[i][1] = (2+index)*1.41*sin((float)(rand()%629)/100) + 96.0;
            particleArray[i][2] = (float)(rand()%(50+1))/1000 + 128.0;
            X = particleArray[i][0] - 96;
            Y = particleArray[i][1] - 96;
            R = sqrt(pow(X,2) + pow(Y,2));
            V = sqrt(1190*R);
            particleArray[i][3] = Y/R*V;
            particleArray[i][4] = X/R*V;
        }
    }
	for (i = 0; i < (int)(M/2); i++) {
		for (j=6;j<7;j++){
			particleArray[i][j] = 0.0; // 0.0 is indicator for Milky Way
		}
	}

	// second galaxy population
	for (i = (int)(M*0.05/2+((10)*0.095*M/2)); i < (int)(M*0.05/2+((10)*0.095*M/2))+(int)(M*0.05/2); i++) {
        particleArray[i][0] = 2*1.41*cos((float)(rand()%629)/100)  + 160.0;
        particleArray[i][1] = 2*1.41*sin((float)(rand()%629)/100)  + 160.0;
        particleArray[i][2] = (float)(rand()%(50+1))/1000 + 128.0;
        X = particleArray[i][0] - 96;
        Y = particleArray[i][1] - 96;
        R = sqrt(pow(X,2) + pow(Y,2));
        V = sqrt(1190*R);
        particleArray[i][3] = Y/R*V;
        particleArray[i][4] = X/R*V;
        particleArray[i][5] = 0;
    }

	for (index=11; index<21; index++){
        for (i = (int)(M*0.05+((index-1)*0.095*M/2)); i < (int)(M*0.05+((index)*0.095*M/2)); i++) {
            particleArray[i][0] = (2+index-10)*1.41*cos((float)(rand()%629)/100)  + 160.0;
            particleArray[i][1] = (2+index-10)*1.41*sin((float)(rand()%629)/100)  + 160.0;
            particleArray[i][2] = (float)(rand()%(150+1))/1000 + 128.0;
            X = particleArray[i][0] - 96;
            Y = particleArray[i][1] - 96;
            R = sqrt(pow(X,2) + pow(Y,2));
            V = sqrt(1190*R);
            particleArray[i][3] = Y/R*V;
            particleArray[i][4] = X/R*V;
            particleArray[i][5] = 0;
        }
    }

	for (i = (int)(M/2); i < M; i++) {
		for (j=6;j<7;j++){
			particleArray[i][j] = 1.0; // 1.0 is indicator for Andromeda
		}
	}
    // printf("%f\n", particleArray[M-1][0]);
    // print the 2D array
	// for (i = 0; i < M; i+=100000)
	// {
	// 	for (j = 0; j < 3; j++) {
    //         printf("%f\n", particleArray[i][j]);

	//    	}
	// 	printf("\n");
	// }

	// // deallocate memory
	// for (i = 0; i < M; i++) {
    //     free(particleArray[i]);
	// }
    // free(particleArray);
}



void cm(float **particleArray){
    float cm[3];
    int i;
    for(i=0;i<M;i++){
        cm[0] = particleArray[i][0];
        cm[1] = particleArray[i][1];
        cm[2] = particleArray[i][2];
    }
    cm[0] = cm[0]/M;
    cm[1] = cm[1]/M;
    cm[2] = cm[2]/M;
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

    printf("updater function initiated\n");
    for(i=0; i<M; i++){
        for(l=0; l<1; l++){
            v_half = particleArray[i][l+3] + 595 +
            (grav_po[(int)round(particleArray[i][0])+1][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])]
                - grav_po[(int)round(particleArray[i][0])-1][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])])/(2);
            x = particleArray[i][l] + v_half;
            v = v_half + 595 +
            (grav_po[(int)round(particleArray[i][0])+1][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])]
                - grav_po[(int)round(particleArray[i][0])-1][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])])/(2);
            particleArray[i][l+3] = v;
        }
        for(l=1; l<2; l++){
            v_half = particleArray[i][l+3] + 595 +
            (grav_po[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])+1][(int)round(particleArray[i][2])] 
            - grav_po[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])-1][(int)round(particleArray[i][2])])/(2);
            x = particleArray[i][l] + v_half;
            v = v_half + 595 +
            (grav_po[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])+1][(int)round(particleArray[i][2])] 
            - grav_po[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])-1][(int)round(particleArray[i][2])])/(2);
            particleArray[i][l+3] = v;
        }
        for(l=2; l<3; l++){
            v_half = particleArray[i][l+3] + 595 +
            (grav_po[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])+1] 
            - grav_po[(int)round(particleArray[i][0])][(int)round(particleArray[i][1])][(int)round(particleArray[i][2])-1])/(2);
            x = particleArray[i][l] + v_half;
            v = v_half + 595 +
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

    printf("all arrays updated\n");
}





// void densArray(float **particleArray, float threedArray[I][J][K]) {
// 	int i, j, k = 0;
//     // dynamically allocate memory of size M*N*O
//     printf("density array intitiation complete\n");

//     // assign values to allocated memory
// 	for (i; i < M; i++) {
//         threedArray[(int)floorf(particleArray[i][0])][(int)floorf(particleArray[i][1])][(int)floorf(particleArray[i][2])] =
//         threedArray[(int)floorf(particleArray[i][0])][(int)floorf(particleArray[i][1])][(int)floorf(particleArray[i][2])] + 1;
// 	}
//     // // print the 3D array
// 	// for (i = 0; i < I; i++)
// 	// {
// 	// 	for (j = 0; j < J; j++)
// 	// 	{
// 	// 		for (k = 0; k < K; k++)
// 	// 			printf("%f\n", threedArray[i][j][k]);
// 	//    	}
// 	// }
//     printf("Density Array completed\n");
// }





void main(){
    // int i;
	float*** grav_po = (float***)malloc(I * sizeof(float**));
    float **particleArray = (float **)malloc(M * sizeof(float *));
    float*** threedArray = (float***)malloc(I * sizeof(float**));
    int i, j, k = 0;
    // dynamically allocate memory of size M*N*O
	// assign values to allocated memory
	for (i = 0; i < I; i++) {
		grav_po[i] = (float**)malloc(J * sizeof(float*));
        if (grav_po[i] == NULL) {
			fprintf(stderr, "Out of memory");
			exit(0);
		}
        for (j = 0; j < J; j++) {
			grav_po[i][j] = (float*)malloc(K * sizeof(float));
            if (grav_po[i][j] == NULL) {
				fprintf(stderr, "Out of memory");
				exit(0);
			}
	}
	
	}
    create_particle_Array(particleArray);
    printf("Particle Array completed\n");
    densArray(particleArray, threedArray);
	center_diff(I, J, K, grav_po, particleArray);
    // densArray(particleArray, threedArray);
    // deallocate memory
	for (i = 0; i < M; i++) {
        free(particleArray[i]);
	}
    free(particleArray);
}