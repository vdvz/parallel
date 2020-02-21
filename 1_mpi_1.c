#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <string.h> 
#define N 5

int main(int argc, char *argv[]){
    double tau = 0.01;
    const double epsilon = 0.00001;
    //const int N = 5;
    double A [N*N];
    double x_n [N];
    double b [N];
    double x_n1[N];
    

    for(int i =0;i<N;i++){
        for(int j =0; j<N;j++){
            *(A + i*N + j) = 1.0;
            if(i==j){
                *( A + i*N + j) = 2.0;
            }
        }
    }

    for(int i =0;i<N;i++){
        b[i] = N+1;
    }
 
    for(int i =0;i<N;i++){
        x_n[i] = 0;
    }

    /*for(int i =0;i<N;i++){
        for(int j =0; j<N;j++){
            printf("%f ",*( A + i*N + j));
        }
        printf("\n");
    }
    
    for(int i =0;i<N;i++){
        printf("%f ", b[i]);
    }
    printf("\n");

    
    for(int i =0;i<N;i++){
        printf("%f ", x_n[i]);
    }
    printf("\n");
    */

    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double norm_Axn_b = 0;
    double norm_b = 0; 
    if(rank==0){
        for (int i = 0; i < floor(N/2); i++){    
           norm_b += b[i]*b[i];
        }
        MPI_Send(&norm_b, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
        MPI_Recv(&norm_b, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("From process: %d get %f", rank, &norm_b);
    }
    if(rank==1){
        double sub_norm;
        for (int i = ceil(N/2); i<N; i++){    
           norm_b += b[i]*b[i];
        }
        MPI_Recv(&sub_norm, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        norm_b += sub_norm;    
        printf("From process: %d send %f", rank, &norm_b);
        MPI_Send(&norm_b, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    }

    norm_b = sqrt(norm_b);
   
    char converged = 0;
    int itter = 0;
    while(!converged){
        if(rank==0){
            for (int i = 0; i<N; i++){
                double sub = 0;
                for(int j = 0; j<floor(N/2); j++){
                    sub += *( A + i*N + j) * x_n[j];
                }
                x_n1[i] = x_n[i] - tau*(sub - b[i]);
            }
            MPI_Send(x_n1, floor(N/2), MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
            MPI_Recv(x_n[(int)ceil(N/2)], ceil(N/2), MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if(rank==1){
            for (int i = 0; i<N; i++){
                double sub = 0;
                for(int j = ceil(N/2); j<N; j++){
                    sub += *( A + i*N + j) * x_n[j];
                }
                x_n1[i] = x_n[i] - tau*(sub - b[i]);
            }
            MPI_Recv(x_n, floor(N/2), MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(x_n1[(int)ceil(N/2)], ceil(N/2), MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        }


    /*
        for (int i = 0; i<N; i++){
            double sub = 0;
            for(int j = 0; j<N; j++){
                sub += *( A + i*N + j) * x_n[j];
            }
            x_n1[i] = x_n[i] - tau*(sub - b[i]);
        }

        memcpy(x_n, x_n1, N*sizeof(double));
    */
    
        double sub_norm_Axn_b = 0;  
        for (int i = 0; i<N; i++){
            double sub = 0;
            for(int j = 0; j<N; j++){
                sub += (*( A + i*N + j)) * x_n[j];
            }
            sub_norm_Axn_b += (sub - b[i])*(sub - b[i]);    
        }
        norm_Axn_b = sqrt(sub_norm_Axn_b);
        
        if(norm_Axn_b/norm_b<epsilon){
            converged = 1;
        }

        itter++;
    }

    printf("Itters: %i\n", itter);

    for(int i =0;i<N;i++){
        printf("%f ", x_n[i]);
    }

    MPI_Finalize();
/*
    free(x_n);
    free(x_n1);
    free(b);
    free(A);
*/  
    return 0;
}