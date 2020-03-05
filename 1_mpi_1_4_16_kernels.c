#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <string.h> 


#define kernels 15
#define N 1000


int main(int argc, char *argv[]){
    double tau = 0.001;
    const double epsilon = 0.000001;
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
    MPI_Status status;
    
    for (int i = 0; i < N; i++){    
           norm_b += b[i]*b[i];
    }

    norm_b = sqrt(norm_b);
    
    int converged = 0;
    int itter = 0;

    int first_element_in_cur_offset_include = rank * (int)(N/kernels);
    int last_elem_in_cur_offset_not_include = (rank+1) * (int)(N/kernels);    

    while(!converged){

        for (int i = first_element_in_cur_offset_include; i<last_elem_in_cur_offset_not_include; i++){
            double sub = 0;
            for(int j = 0; j<N; j++){
                sub += *( A + i*N + j) * x_n[j];
            }
            x_n1[i] = x_n[i] - tau*(sub - b[i]);
            
        }
        
        for(int thread = 0; thread<kernels; thread++){
            if(thread != rank){
                MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                            last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                    MPI_DOUBLE, thread, 1, MPI_COMM_WORLD);
            }
        }


        for(int thread = 0; thread<kernels; thread++){
            if(thread != rank){
                MPI_Recv(&x_n1[thread*(int)N/kernels],
                        (thread+1)*(int)N/kernels - thread*(int)N/kernels, MPI_DOUBLE, thread, 1, MPI_COMM_WORLD, &status);
            }
        }


        memcpy(x_n, x_n1, N*sizeof(double));


    /*  printf("rank: %i\n", rank);
        for(int i =0;i<N;i++){
            printf("%f ", x_n[i]);
        }
        printf("\n");
    */

        double sub_norm_Axn_b = 0;  
        for (int i = first_element_in_cur_offset_include; i<last_elem_in_cur_offset_not_include; i++){
            double sub = 0;
            for(int j = 0; j<N; j++){
                sub += (*( A + i*N + j)) * x_n[j];
            }
            sub_norm_Axn_b += (sub - b[i])*(sub - b[i]);    
        }
        
//      if(rank==0){
         
        double norm_from_another_stream;

        for(int thread = 0; thread<kernels; thread++){
            if(thread != rank){
                MPI_Send(&sub_norm_Axn_b, 1, MPI_DOUBLE, thread, 1, MPI_COMM_WORLD);
            }
        }

        for(int thread = 0; thread<kernels; thread++){
            if(thread != rank){
                MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, thread, 1, MPI_COMM_WORLD, &status);
                sub_norm_Axn_b += norm_from_another_stream;
            }
        }


        norm_Axn_b = sqrt(sub_norm_Axn_b);    

        if(norm_Axn_b/norm_b<epsilon){
            converged = 1;
        }    

    /*
        for(int thread = 0; thread<kernels; thread++){
            if(thread != rank){
                MPI_Send(&converged, 1, MPI_INT, thread, 1, MPI_COMM_WORLD);
            }
        }
    */

//        }

        /*if(rank!=0){

            MPI_Send(&sub_norm_Axn_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        
            MPI_Recv(&converged, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

            if(converged == 1) break;

        }*/
    
        itter++;
        if(itter == 100) break;

    }

    if(rank==0){
        printf("Itters: %i\n", itter);
        for(int i =0;i<N;i++){
            printf("%f ", x_n[i]);
        }
    }

    MPI_Finalize();
 
    return 0;
}