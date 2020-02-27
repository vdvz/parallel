#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <string.h> 


#define kernels 2
#define N 10


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

         if(rank==0){
            
            for (int i = first_element_in_cur_offset_include; i<last_elem_in_cur_offset_not_include; i++){
                double sub = 0;
                for(int j = 0; j<N; j++){
                    sub += *( A + i*N + j) * x_n[j];
                }
                x_n1[i] = x_n[i] - tau*(sub - b[i]);
                
            }

           MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);

            
            #if  kernels>2
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 2, 1, MPI_COMM_WORLD);
         

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 3, 1, MPI_COMM_WORLD);

            

            #endif
            
            #if kernels>4

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 4, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 5, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 6, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 7, 1, MPI_COMM_WORLD);

            #endif

            #if kernels > 8

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 8, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 9, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 10, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 11, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 12, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 13, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 14, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 15, 1, MPI_COMM_WORLD);

            #endif
            
            MPI_Recv(&x_n1[1*(int)N/kernels],
                            2*(int)N/kernels - 1*(int)N/kernels, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
            

            #if kernels > 2

            MPI_Recv(&x_n1[2*(int)N/kernels],
                            3*(int)N/kernels - 2*(int)N/kernels, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[3*(int)N/kernels],
                            4*(int)N/kernels - 3*(int)N/kernels, MPI_DOUBLE, 3, 1, MPI_COMM_WORLD, &status);
            
            #endif

            
            #if kernels > 4

            MPI_Recv(&x_n1[4*(int)N/kernels],
                            5*(int)N/kernels - 4*(int)N/kernels, MPI_DOUBLE, 4, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[5*(int)N/kernels],
                            6*(int)N/kernels - 5*(int)N/kernels, MPI_DOUBLE, 5, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[6*(int)N/kernels],
                            7*(int)N/kernels - 6*(int)N/kernels, MPI_DOUBLE, 6, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[7*(int)N/kernels],
                            8*(int)N/kernels - 7*(int)N/kernels, MPI_DOUBLE, 7, 1, MPI_COMM_WORLD, &status);

            #endif

            #if kernels > 8

            MPI_Recv(&x_n1[8*(int)N/kernels],
                            9*(int)N/kernels - 8*(int)N/kernels, MPI_DOUBLE, 8, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[9*(int)N/kernels],
                             10*(int)N/kernels - 9*(int)N/kernels, MPI_DOUBLE, 9, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[10*(int)N/kernels],
                            11*(int)N/kernels - 10*(int)N/kernels, MPI_DOUBLE, 10, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[11*(int)N/kernels],
                            12*(int)N/kernels - 11*(int)N/kernels, MPI_DOUBLE, 11, 1, MPI_COMM_WORLD, &status);
            
            MPI_Recv(&x_n1[12*(int)N/kernels],
                            13*(int)N/kernels - 12*(int)N/kernels, MPI_DOUBLE, 12, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[13*(int)N/kernels],
                            14*(int)N/kernels - 13*(int)N/kernels, MPI_DOUBLE, 13, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[14*(int)N/kernels],
                            15*(int)N/kernels - 14*(int)N/kernels, MPI_DOUBLE, 14, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[15*(int)N/kernels],
                            16*(int)N/kernels - 15*(int)N/kernels, MPI_DOUBLE, 15, 1, MPI_COMM_WORLD, &status);

            #endif

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

            double norm_from_another_stream;
    
            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;   
            
            #if kernels > 2
            
            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;
            
            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 3, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;
            #endif

            #if kernels > 4
                        
            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 4, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;
            
            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 5, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;
            
            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 6, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;
            
            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 7, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;

            #endif
            
            #if kernels > 8

            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 8, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;
            
            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 9, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;
            
            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 10, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;
            
            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 11, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;
            
            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 12, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;
            
            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 13, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;
            
            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 14, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;
            
            MPI_Recv(&norm_from_another_stream, 1, MPI_DOUBLE, 15, 1, MPI_COMM_WORLD, &status);
            sub_norm_Axn_b += norm_from_another_stream;
            
            #endif
            
            norm_Axn_b = sqrt(sub_norm_Axn_b);    
            
            
            if(norm_Axn_b/norm_b<epsilon){
                converged = 1;
            }    
            
            MPI_Send(&converged, 1, MPI_INT, 1, 1, MPI_COMM_WORLD);
            
            
            #if kernels > 2
            
            MPI_Send(&converged, 1, MPI_INT, 2, 1, MPI_COMM_WORLD);
            
            MPI_Send(&converged, 1, MPI_INT, 3, 1, MPI_COMM_WORLD);
            
            #endif

            #if kernels > 4 
            
            MPI_Send(&converged, 1, MPI_INT, 4, 1, MPI_COMM_WORLD);
            
            MPI_Send(&converged, 1, MPI_INT, 5, 1, MPI_COMM_WORLD);
            
            MPI_Send(&converged, 1, MPI_INT, 6, 1, MPI_COMM_WORLD);
            
            MPI_Send(&converged, 1, MPI_INT, 7, 1, MPI_COMM_WORLD);
        
            #endif

            #if kernels > 8

            MPI_Send(&converged, 1, MPI_INT, 8, 1, MPI_COMM_WORLD);

            MPI_Send(&converged, 1, MPI_INT, 9, 1, MPI_COMM_WORLD);
            
            MPI_Send(&converged, 1, MPI_INT, 10, 1, MPI_COMM_WORLD);
            
            MPI_Send(&converged, 1, MPI_INT, 11, 1, MPI_COMM_WORLD);
            
            MPI_Send(&converged, 1, MPI_INT, 12, 1, MPI_COMM_WORLD);
        
            MPI_Send(&converged, 1, MPI_INT, 13, 1, MPI_COMM_WORLD);
            
            MPI_Send(&converged, 1, MPI_INT, 14, 1, MPI_COMM_WORLD);
            
            MPI_Send(&converged, 1, MPI_INT, 15, 1, MPI_COMM_WORLD);
        
            #endif

        }

        if(rank==1){

            for (int i = first_element_in_cur_offset_include; i<last_elem_in_cur_offset_not_include; i++){
                double sub = 0;
                for(int j = 0; j<N; j++){
                    sub += *( A + i*N + j) * x_n[j];
                }
                x_n1[i] = x_n[i] - tau*(sub - b[i]);
            }
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

            #if  kernels>2
                
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 2, 1, MPI_COMM_WORLD);

            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 3, 1, MPI_COMM_WORLD);

            #endif
            
            #if kernels>4

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 4, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 5, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 6, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 7, 1, MPI_COMM_WORLD);

            #endif

            #if kernels > 8

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 8, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 9, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 10, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 11, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 12, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 13, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 14, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 15, 1, MPI_COMM_WORLD);

            #endif

            MPI_Recv(&x_n1[0*(int)N/kernels],
                            1*(int)N/kernels - 0*(int)N/kernels, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

            
            #if kernels > 2

            MPI_Recv(&x_n1[2*(int)N/kernels],
                            3*(int)N/kernels - 2*(int)N/kernels, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[3*(int)N/kernels],
                            4*(int)N/kernels - 3*(int)N/kernels, MPI_DOUBLE, 3, 1, MPI_COMM_WORLD, &status);
            
            #endif

            #if kernels > 4

            MPI_Recv(&x_n1[4*(int)N/kernels],
                            5*(int)N/kernels - 4*(int)N/kernels, MPI_DOUBLE, 4, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[5*(int)N/kernels],
                            6*(int)N/kernels - 5*(int)N/kernels, MPI_DOUBLE, 5, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[6*(int)N/kernels],
                            7*(int)N/kernels - 6*(int)N/kernels, MPI_DOUBLE, 6, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[7*(int)N/kernels],
                            8*(int)N/kernels - 7*(int)N/kernels, MPI_DOUBLE, 7, 1, MPI_COMM_WORLD, &status);

            #endif

            #if kernels > 8

            MPI_Recv(&x_n1[8*(int)N/kernels],
                            9*(int)N/kernels - 8*(int)N/kernels, MPI_DOUBLE, 8, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[9*(int)N/kernels],
                             10*(int)N/kernels - 9*(int)N/kernels, MPI_DOUBLE, 9, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[10*(int)N/kernels],
                            11*(int)N/kernels - 10*(int)N/kernels, MPI_DOUBLE, 10, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[11*(int)N/kernels],
                            12*(int)N/kernels - 11*(int)N/kernels, MPI_DOUBLE, 11, 1, MPI_COMM_WORLD, &status);
            
            MPI_Recv(&x_n1[12*(int)N/kernels],
                            13*(int)N/kernels - 12*(int)N/kernels, MPI_DOUBLE, 12, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[13*(int)N/kernels],
                            14*(int)N/kernels - 13*(int)N/kernels, MPI_DOUBLE, 13, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[14*(int)N/kernels],
                            15*(int)N/kernels - 14*(int)N/kernels, MPI_DOUBLE, 14, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[15*(int)N/kernels],
                            16*(int)N/kernels - 15*(int)N/kernels, MPI_DOUBLE, 15, 1, MPI_COMM_WORLD, &status);

            #endif

            memcpy(x_n, x_n1, N*sizeof(double));

       /*   printf("rank: %i\n", rank);
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

            MPI_Send(&sub_norm_Axn_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            
            MPI_Recv(&converged, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);


        }

    #if kernels > 2
         if(rank==2){

            for (int i = first_element_in_cur_offset_include; i<last_elem_in_cur_offset_not_include; i++){
                double sub = 0;
                for(int j = 0; j<N; j++){
                    sub += *( A + i*N + j) * x_n[j];
                }
                x_n1[i] = x_n[i] - tau*(sub - b[i]);
            }
                        
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

            #if  kernels>2
                
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);

            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 3, 1, MPI_COMM_WORLD);

            #endif
            
            #if kernels>4

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 4, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 5, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 6, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 7, 1, MPI_COMM_WORLD);

            #endif

            #if kernels > 8

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 8, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 9, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 10, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 11, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 12, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 13, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 14, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 15, 1, MPI_COMM_WORLD);

            #endif
            
            
            MPI_Recv(&x_n1[0*(int)N/kernels],
                            1*(int)N/kernels - 0*(int)N/kernels, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);


            #if kernels > 2

            MPI_Recv(&x_n1[1*(int)N/kernels],
                            2*(int)N/kernels - 1*(int)N/kernels, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[3*(int)N/kernels],
                            4*(int)N/kernels - 3*(int)N/kernels, MPI_DOUBLE, 3, 1, MPI_COMM_WORLD, &status);
            
            #endif

            #if kernels > 4

            MPI_Recv(&x_n1[4*(int)N/kernels],
                            5*(int)N/kernels - 4*(int)N/kernels, MPI_DOUBLE, 4, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[5*(int)N/kernels],
                            6*(int)N/kernels - 5*(int)N/kernels, MPI_DOUBLE, 5, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[6*(int)N/kernels],
                            7*(int)N/kernels - 6*(int)N/kernels, MPI_DOUBLE, 6, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[7*(int)N/kernels],
                            8*(int)N/kernels - 7*(int)N/kernels, MPI_DOUBLE, 7, 1, MPI_COMM_WORLD, &status);

            #endif

            #if kernels > 8

            MPI_Recv(&x_n1[8*(int)N/kernels],
                            9*(int)N/kernels - 8*(int)N/kernels, MPI_DOUBLE, 8, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[9*(int)N/kernels],
                             10*(int)N/kernels - 9*(int)N/kernels, MPI_DOUBLE, 9, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[10*(int)N/kernels],
                            11*(int)N/kernels - 10*(int)N/kernels, MPI_DOUBLE, 10, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[11*(int)N/kernels],
                            12*(int)N/kernels - 11*(int)N/kernels, MPI_DOUBLE, 11, 1, MPI_COMM_WORLD, &status);
            
            MPI_Recv(&x_n1[12*(int)N/kernels],
                            13*(int)N/kernels - 12*(int)N/kernels, MPI_DOUBLE, 12, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[13*(int)N/kernels],
                            14*(int)N/kernels - 13*(int)N/kernels, MPI_DOUBLE, 13, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[14*(int)N/kernels],
                            15*(int)N/kernels - 14*(int)N/kernels, MPI_DOUBLE, 14, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[15*(int)N/kernels],
                            16*(int)N/kernels - 15*(int)N/kernels, MPI_DOUBLE, 15, 1, MPI_COMM_WORLD, &status);

            #endif

            memcpy(x_n, x_n1, N*sizeof(double));

       /*   printf("rank: %i\n", rank);
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

            MPI_Send(&sub_norm_Axn_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            MPI_Recv(&converged, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);


        }
    
        if(rank==3){

            for (int i = first_element_in_cur_offset_include; i<last_elem_in_cur_offset_not_include; i++){
                double sub = 0;
                for(int j = 0; j<N; j++){
                    sub += *( A + i*N + j) * x_n[j];
                }
                x_n1[i] = x_n[i] - tau*(sub - b[i]);
            }
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

            #if  kernels>2
                
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);

            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 2, 1, MPI_COMM_WORLD);

            #endif
            
            #if kernels>4

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 4, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 5, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 6, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 7, 1, MPI_COMM_WORLD);

            #endif

            #if kernels > 8

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 8, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 9, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 10, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 11, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 12, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 13, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 14, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 15, 1, MPI_COMM_WORLD);

            #endif
            
            MPI_Recv(&x_n1[0*(int)N/kernels],
                            1*(int)N/kernels - 0*(int)N/kernels, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

            #if kernels > 2

            MPI_Recv(&x_n1[1*(int)N/kernels],
                            2*(int)N/kernels - 1*(int)N/kernels, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[2*(int)N/kernels],
                            3*(int)N/kernels - 2*(int)N/kernels, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, &status);
            
            #endif

            #if kernels > 4

            MPI_Recv(&x_n1[4*(int)N/kernels],
                            5*(int)N/kernels - 4*(int)N/kernels, MPI_DOUBLE, 4, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[5*(int)N/kernels],
                            6*(int)N/kernels - 5*(int)N/kernels, MPI_DOUBLE, 5, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[6*(int)N/kernels],
                            7*(int)N/kernels - 6*(int)N/kernels, MPI_DOUBLE, 6, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[7*(int)N/kernels],
                            8*(int)N/kernels - 7*(int)N/kernels, MPI_DOUBLE, 7, 1, MPI_COMM_WORLD, &status);

            #endif

            #if kernels > 8

            MPI_Recv(&x_n1[8*(int)N/kernels],
                            9*(int)N/kernels - 8*(int)N/kernels, MPI_DOUBLE, 8, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[9*(int)N/kernels],
                             10*(int)N/kernels - 9*(int)N/kernels, MPI_DOUBLE, 9, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[10*(int)N/kernels],
                            11*(int)N/kernels - 10*(int)N/kernels, MPI_DOUBLE, 10, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[11*(int)N/kernels],
                            12*(int)N/kernels - 11*(int)N/kernels, MPI_DOUBLE, 11, 1, MPI_COMM_WORLD, &status);
            
            MPI_Recv(&x_n1[12*(int)N/kernels],
                            13*(int)N/kernels - 12*(int)N/kernels, MPI_DOUBLE, 12, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[13*(int)N/kernels],
                            14*(int)N/kernels - 13*(int)N/kernels, MPI_DOUBLE, 13, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[14*(int)N/kernels],
                            15*(int)N/kernels - 14*(int)N/kernels, MPI_DOUBLE, 14, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[15*(int)N/kernels],
                            16*(int)N/kernels - 15*(int)N/kernels, MPI_DOUBLE, 15, 1, MPI_COMM_WORLD, &status);

            #endif

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

            MPI_Send(&sub_norm_Axn_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            
            MPI_Recv(&converged, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

        
        }
    #endif
    #if kernels>4
        if(rank==4){

            for (int i = first_element_in_cur_offset_include; i<last_elem_in_cur_offset_not_include; i++){
                double sub = 0;
                for(int j = 0; j<N; j++){
                    sub += *( A + i*N + j) * x_n[j];
                }
                x_n1[i] = x_n[i] - tau*(sub - b[i]);
            }
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

            #if  kernels>2
                
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);

            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 2, 1, MPI_COMM_WORLD);

            #endif
            
            #if kernels>4

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 3, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 5, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 6, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 7, 1, MPI_COMM_WORLD);

            #endif

            #if kernels > 8

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 8, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 9, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 10, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 11, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 12, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 13, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 14, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 15, 1, MPI_COMM_WORLD);

            #endif
            
            MPI_Recv(&x_n1[0*(int)N/kernels],
                            1*(int)N/kernels - 0*(int)N/kernels, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

            #if kernels > 2

            MPI_Recv(&x_n1[1*(int)N/kernels],
                            2*(int)N/kernels - 1*(int)N/kernels, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[2*(int)N/kernels],
                            3*(int)N/kernels - 2*(int)N/kernels, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, &status);
            
            #endif

            #if kernels > 4

            MPI_Recv(&x_n1[3*(int)N/kernels],
                            4*(int)N/kernels - 3*(int)N/kernels, MPI_DOUBLE, 3, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[5*(int)N/kernels],
                            6*(int)N/kernels - 5*(int)N/kernels, MPI_DOUBLE, 5, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[6*(int)N/kernels],
                            7*(int)N/kernels - 6*(int)N/kernels, MPI_DOUBLE, 6, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[7*(int)N/kernels],
                            8*(int)N/kernels - 7*(int)N/kernels, MPI_DOUBLE, 7, 1, MPI_COMM_WORLD, &status);

            #endif

            #if kernels > 8

            MPI_Recv(&x_n1[8*(int)N/kernels],
                            9*(int)N/kernels - 8*(int)N/kernels, MPI_DOUBLE, 8, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[9*(int)N/kernels],
                             10*(int)N/kernels - 9*(int)N/kernels, MPI_DOUBLE, 9, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[10*(int)N/kernels],
                            11*(int)N/kernels - 10*(int)N/kernels, MPI_DOUBLE, 10, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[11*(int)N/kernels],
                            12*(int)N/kernels - 11*(int)N/kernels, MPI_DOUBLE, 11, 1, MPI_COMM_WORLD, &status);
            
            MPI_Recv(&x_n1[12*(int)N/kernels],
                            13*(int)N/kernels - 12*(int)N/kernels, MPI_DOUBLE, 12, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[13*(int)N/kernels],
                            14*(int)N/kernels - 13*(int)N/kernels, MPI_DOUBLE, 13, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[14*(int)N/kernels],
                            15*(int)N/kernels - 14*(int)N/kernels, MPI_DOUBLE, 14, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[15*(int)N/kernels],
                            16*(int)N/kernels - 15*(int)N/kernels, MPI_DOUBLE, 15, 1, MPI_COMM_WORLD, &status);

            #endif

            memcpy(x_n, x_n1, N*sizeof(double));

       /*   printf("rank: %i\n", rank);
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

            MPI_Send(&sub_norm_Axn_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            
            MPI_Recv(&converged, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

        }
        if(rank==5){

            for (int i = first_element_in_cur_offset_include; i<last_elem_in_cur_offset_not_include; i++){
                double sub = 0;
                for(int j = 0; j<N; j++){
                    sub += *( A + i*N + j) * x_n[j];
                }
                x_n1[i] = x_n[i] - tau*(sub - b[i]);
            }
            

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

            #if  kernels>2
                
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);

            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 2, 1, MPI_COMM_WORLD);

            #endif
            
            #if kernels>4

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 3, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 4, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 6, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 7, 1, MPI_COMM_WORLD);

            #endif

            #if kernels > 8

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 8, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 9, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 10, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 11, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 12, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 13, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 14, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 15, 1, MPI_COMM_WORLD);

            #endif
            MPI_Recv(&x_n1[0*(int)N/kernels],
                            1*(int)N/kernels - 0*(int)N/kernels, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

            #if kernels > 2

            MPI_Recv(&x_n1[1*(int)N/kernels],
                            2*(int)N/kernels - 1*(int)N/kernels, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[2*(int)N/kernels],
                            3*(int)N/kernels - 2*(int)N/kernels, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, &status);
            
            #endif

            #if kernels > 4

            MPI_Recv(&x_n1[3*(int)N/kernels],
                            4*(int)N/kernels - 3*(int)N/kernels, MPI_DOUBLE, 3, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[4*(int)N/kernels],
                            5*(int)N/kernels - 4*(int)N/kernels, MPI_DOUBLE, 4, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[6*(int)N/kernels],
                            7*(int)N/kernels - 6*(int)N/kernels, MPI_DOUBLE, 6, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[7*(int)N/kernels],
                            8*(int)N/kernels - 7*(int)N/kernels, MPI_DOUBLE, 7, 1, MPI_COMM_WORLD, &status);

            #endif

            #if kernels > 8

            MPI_Recv(&x_n1[8*(int)N/kernels],
                            9*(int)N/kernels - 8*(int)N/kernels, MPI_DOUBLE, 8, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[9*(int)N/kernels],
                             10*(int)N/kernels - 9*(int)N/kernels, MPI_DOUBLE, 9, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[10*(int)N/kernels],
                            11*(int)N/kernels - 10*(int)N/kernels, MPI_DOUBLE, 10, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[11*(int)N/kernels],
                            12*(int)N/kernels - 11*(int)N/kernels, MPI_DOUBLE, 11, 1, MPI_COMM_WORLD, &status);
            
            MPI_Recv(&x_n1[12*(int)N/kernels],
                            13*(int)N/kernels - 12*(int)N/kernels, MPI_DOUBLE, 12, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[13*(int)N/kernels],
                            14*(int)N/kernels - 13*(int)N/kernels, MPI_DOUBLE, 13, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[14*(int)N/kernels],
                            15*(int)N/kernels - 14*(int)N/kernels, MPI_DOUBLE, 14, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[15*(int)N/kernels],
                            16*(int)N/kernels - 15*(int)N/kernels, MPI_DOUBLE, 15, 1, MPI_COMM_WORLD, &status);

            #endif
            

            memcpy(x_n, x_n1, N*sizeof(double));

       /*   printf("rank: %i\n", rank);
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

            MPI_Send(&sub_norm_Axn_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            
            MPI_Recv(&converged, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

        }
        if(rank==6){

            for (int i = first_element_in_cur_offset_include; i<last_elem_in_cur_offset_not_include; i++){
                double sub = 0;
                for(int j = 0; j<N; j++){
                    sub += *( A + i*N + j) * x_n[j];
                }
                x_n1[i] = x_n[i] - tau*(sub - b[i]);
            }
            
       
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

            #if  kernels>2
                
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);

            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 2, 1, MPI_COMM_WORLD);

            #endif
            
            #if kernels>4

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 3, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 4, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 5, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 7, 1, MPI_COMM_WORLD);

            #endif

            #if kernels > 8

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 8, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 9, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 10, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 11, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 12, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 13, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 14, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 15, 1, MPI_COMM_WORLD);

            #endif

           MPI_Recv(&x_n1[0*(int)N/kernels],
                            1*(int)N/kernels - 0*(int)N/kernels, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

            #if kernels > 2

            MPI_Recv(&x_n1[1*(int)N/kernels],
                            2*(int)N/kernels - 1*(int)N/kernels, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[2*(int)N/kernels],
                            3*(int)N/kernels - 2*(int)N/kernels, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, &status);
            
            #endif

            #if kernels > 4

            MPI_Recv(&x_n1[3*(int)N/kernels],
                            4*(int)N/kernels - 3*(int)N/kernels, MPI_DOUBLE, 3, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[4*(int)N/kernels],
                            5*(int)N/kernels - 4*(int)N/kernels, MPI_DOUBLE, 4, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[5*(int)N/kernels],
                            6*(int)N/kernels - 5*(int)N/kernels, MPI_DOUBLE, 5, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[7*(int)N/kernels],
                            8*(int)N/kernels - 7*(int)N/kernels, MPI_DOUBLE, 7, 1, MPI_COMM_WORLD, &status);

            #endif

            #if kernels > 8

            MPI_Recv(&x_n1[8*(int)N/kernels],
                            9*(int)N/kernels - 8*(int)N/kernels, MPI_DOUBLE, 8, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[9*(int)N/kernels],
                             10*(int)N/kernels - 9*(int)N/kernels, MPI_DOUBLE, 9, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[10*(int)N/kernels],
                            11*(int)N/kernels - 10*(int)N/kernels, MPI_DOUBLE, 10, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[11*(int)N/kernels],
                            12*(int)N/kernels - 11*(int)N/kernels, MPI_DOUBLE, 11, 1, MPI_COMM_WORLD, &status);
            
            MPI_Recv(&x_n1[12*(int)N/kernels],
                            13*(int)N/kernels - 12*(int)N/kernels, MPI_DOUBLE, 12, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[13*(int)N/kernels],
                            14*(int)N/kernels - 13*(int)N/kernels, MPI_DOUBLE, 13, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[14*(int)N/kernels],
                            15*(int)N/kernels - 14*(int)N/kernels, MPI_DOUBLE, 14, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[15*(int)N/kernels],
                            16*(int)N/kernels - 15*(int)N/kernels, MPI_DOUBLE, 15, 1, MPI_COMM_WORLD, &status);

            #endif
         

            memcpy(x_n, x_n1, N*sizeof(double));

       /*   printf("rank: %i\n", rank);
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

            MPI_Send(&sub_norm_Axn_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            
            MPI_Recv(&converged, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

        }

        if(rank==7){

            for (int i = first_element_in_cur_offset_include; i<last_elem_in_cur_offset_not_include; i++){
                double sub = 0;
                for(int j = 0; j<N; j++){
                    sub += *( A + i*N + j) * x_n[j];
                }
                x_n1[i] = x_n[i] - tau*(sub - b[i]);
            }
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

            #if  kernels>2
                
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);

            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 2, 1, MPI_COMM_WORLD);

            #endif
            
            #if kernels>4

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 3, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 4, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 5, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 6, 1, MPI_COMM_WORLD);

            #endif

            #if kernels > 8

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 8, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 9, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 10, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 11, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 12, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 13, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 14, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 15, 1, MPI_COMM_WORLD);

            #endif

            MPI_Recv(&x_n1[0*(int)N/kernels],
                            1*(int)N/kernels - 0*(int)N/kernels, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

            #if kernels > 2

            MPI_Recv(&x_n1[1*(int)N/kernels],
                            2*(int)N/kernels - 1*(int)N/kernels, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[2*(int)N/kernels],
                            3*(int)N/kernels - 2*(int)N/kernels, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, &status);
            
            #endif

            #if kernels > 4

            MPI_Recv(&x_n1[3*(int)N/kernels],
                            4*(int)N/kernels - 3*(int)N/kernels, MPI_DOUBLE, 3, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[4*(int)N/kernels],
                            5*(int)N/kernels - 4*(int)N/kernels, MPI_DOUBLE, 4, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[5*(int)N/kernels],
                            6*(int)N/kernels - 5*(int)N/kernels, MPI_DOUBLE, 5, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[6*(int)N/kernels],
                            7*(int)N/kernels - 6*(int)N/kernels, MPI_DOUBLE, 6, 1, MPI_COMM_WORLD, &status);

            #endif

            #if kernels > 8

            MPI_Recv(&x_n1[8*(int)N/kernels],
                            9*(int)N/kernels - 8*(int)N/kernels, MPI_DOUBLE, 8, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[9*(int)N/kernels],
                             10*(int)N/kernels - 9*(int)N/kernels, MPI_DOUBLE, 9, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[10*(int)N/kernels],
                            11*(int)N/kernels - 10*(int)N/kernels, MPI_DOUBLE, 10, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[11*(int)N/kernels],
                            12*(int)N/kernels - 11*(int)N/kernels, MPI_DOUBLE, 11, 1, MPI_COMM_WORLD, &status);
            
            MPI_Recv(&x_n1[12*(int)N/kernels],
                            13*(int)N/kernels - 12*(int)N/kernels, MPI_DOUBLE, 12, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[13*(int)N/kernels],
                            14*(int)N/kernels - 13*(int)N/kernels, MPI_DOUBLE, 13, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[14*(int)N/kernels],
                            15*(int)N/kernels - 14*(int)N/kernels, MPI_DOUBLE, 14, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[15*(int)N/kernels],
                            16*(int)N/kernels - 15*(int)N/kernels, MPI_DOUBLE, 15, 1, MPI_COMM_WORLD, &status);

            #endif
    

            memcpy(x_n, x_n1, N*sizeof(double));

       /*   printf("rank: %i\n", rank);
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

            MPI_Send(&sub_norm_Axn_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            
            MPI_Recv(&converged, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

        }
        
    #endif

    #if kernels>8

            if(rank==8){

            for (int i = first_element_in_cur_offset_include; i<last_elem_in_cur_offset_not_include; i++){
                double sub = 0;
                for(int j = 0; j<N; j++){
                    sub += *( A + i*N + j) * x_n[j];
                }
                x_n1[i] = x_n[i] - tau*(sub - b[i]);
            }
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

            #if  kernels>2
                
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);

            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 2, 1, MPI_COMM_WORLD);

            #endif
            
            #if kernels>4

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 3, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 4, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 5, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 6, 1, MPI_COMM_WORLD);

            #endif

            #if kernels > 8

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 7, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 9, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 10, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 11, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 12, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 13, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 14, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 15, 1, MPI_COMM_WORLD);

            #endif

            MPI_Recv(&x_n1[0*(int)N/kernels],
                            1*(int)N/kernels - 0*(int)N/kernels, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

            #if kernels > 2

            MPI_Recv(&x_n1[1*(int)N/kernels],
                            2*(int)N/kernels - 1*(int)N/kernels, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[2*(int)N/kernels],
                            3*(int)N/kernels - 2*(int)N/kernels, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, &status);
            
            #endif

            #if kernels > 4

            MPI_Recv(&x_n1[3*(int)N/kernels],
                            4*(int)N/kernels - 3*(int)N/kernels, MPI_DOUBLE, 3, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[4*(int)N/kernels],
                            5*(int)N/kernels - 4*(int)N/kernels, MPI_DOUBLE, 4, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[5*(int)N/kernels],
                            6*(int)N/kernels - 5*(int)N/kernels, MPI_DOUBLE, 5, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[6*(int)N/kernels],
                            7*(int)N/kernels - 6*(int)N/kernels, MPI_DOUBLE, 6, 1, MPI_COMM_WORLD, &status);

            #endif

            #if kernels > 8

            MPI_Recv(&x_n1[7*(int)N/kernels],
                            8*(int)N/kernels - 7*(int)N/kernels, MPI_DOUBLE, 7, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[9*(int)N/kernels],
                             10*(int)N/kernels - 9*(int)N/kernels, MPI_DOUBLE, 9, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[10*(int)N/kernels],
                            11*(int)N/kernels - 10*(int)N/kernels, MPI_DOUBLE, 10, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[11*(int)N/kernels],
                            12*(int)N/kernels - 11*(int)N/kernels, MPI_DOUBLE, 11, 1, MPI_COMM_WORLD, &status);
            
            MPI_Recv(&x_n1[12*(int)N/kernels],
                            13*(int)N/kernels - 12*(int)N/kernels, MPI_DOUBLE, 12, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[13*(int)N/kernels],
                            14*(int)N/kernels - 13*(int)N/kernels, MPI_DOUBLE, 13, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[14*(int)N/kernels],
                            15*(int)N/kernels - 14*(int)N/kernels, MPI_DOUBLE, 14, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[15*(int)N/kernels],
                            16*(int)N/kernels - 15*(int)N/kernels, MPI_DOUBLE, 15, 1, MPI_COMM_WORLD, &status);

            #endif
    

            memcpy(x_n, x_n1, N*sizeof(double));

       /*   printf("rank: %i\n", rank);
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

            MPI_Send(&sub_norm_Axn_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            
            MPI_Recv(&converged, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

        }

                if(rank==9){

            for (int i = first_element_in_cur_offset_include; i<last_elem_in_cur_offset_not_include; i++){
                double sub = 0;
                for(int j = 0; j<N; j++){
                    sub += *( A + i*N + j) * x_n[j];
                }
                x_n1[i] = x_n[i] - tau*(sub - b[i]);
            }
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

            #if  kernels>2
                
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);

            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 2, 1, MPI_COMM_WORLD);

            #endif
            
            #if kernels>4

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 3, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 4, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 5, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 6, 1, MPI_COMM_WORLD);

            #endif

            #if kernels > 8

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 7, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 8, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 10, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 11, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 12, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 13, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 14, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 15, 1, MPI_COMM_WORLD);

            #endif

            MPI_Recv(&x_n1[0*(int)N/kernels],
                            1*(int)N/kernels - 0*(int)N/kernels, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

            #if kernels > 2

            MPI_Recv(&x_n1[1*(int)N/kernels],
                            2*(int)N/kernels - 1*(int)N/kernels, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[2*(int)N/kernels],
                            3*(int)N/kernels - 2*(int)N/kernels, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, &status);
            
            #endif

            #if kernels > 4

            MPI_Recv(&x_n1[3*(int)N/kernels],
                            4*(int)N/kernels - 3*(int)N/kernels, MPI_DOUBLE, 3, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[4*(int)N/kernels],
                            5*(int)N/kernels - 4*(int)N/kernels, MPI_DOUBLE, 4, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[5*(int)N/kernels],
                            6*(int)N/kernels - 5*(int)N/kernels, MPI_DOUBLE, 5, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[6*(int)N/kernels],
                            7*(int)N/kernels - 6*(int)N/kernels, MPI_DOUBLE, 6, 1, MPI_COMM_WORLD, &status);

            #endif

            #if kernels > 8

            MPI_Recv(&x_n1[7*(int)N/kernels],
                            8*(int)N/kernels - 7*(int)N/kernels, MPI_DOUBLE, 7, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[8*(int)N/kernels],
                             9*(int)N/kernels - 8*(int)N/kernels, MPI_DOUBLE, 8, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[10*(int)N/kernels],
                            11*(int)N/kernels - 10*(int)N/kernels, MPI_DOUBLE, 10, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[11*(int)N/kernels],
                            12*(int)N/kernels - 11*(int)N/kernels, MPI_DOUBLE, 11, 1, MPI_COMM_WORLD, &status);
            
            MPI_Recv(&x_n1[12*(int)N/kernels],
                            13*(int)N/kernels - 12*(int)N/kernels, MPI_DOUBLE, 12, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[13*(int)N/kernels],
                            14*(int)N/kernels - 13*(int)N/kernels, MPI_DOUBLE, 13, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[14*(int)N/kernels],
                            15*(int)N/kernels - 14*(int)N/kernels, MPI_DOUBLE, 14, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[15*(int)N/kernels],
                            16*(int)N/kernels - 15*(int)N/kernels, MPI_DOUBLE, 15, 1, MPI_COMM_WORLD, &status);

            #endif
    

            memcpy(x_n, x_n1, N*sizeof(double));

       /*   printf("rank: %i\n", rank);
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

            MPI_Send(&sub_norm_Axn_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            
            MPI_Recv(&converged, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

        }

                if(rank==10){

            for (int i = first_element_in_cur_offset_include; i<last_elem_in_cur_offset_not_include; i++){
                double sub = 0;
                for(int j = 0; j<N; j++){
                    sub += *( A + i*N + j) * x_n[j];
                }
                x_n1[i] = x_n[i] - tau*(sub - b[i]);
            }
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

            #if  kernels>2
                
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);

            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 2, 1, MPI_COMM_WORLD);

            #endif
            
            #if kernels>4

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 3, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 4, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 5, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 6, 1, MPI_COMM_WORLD);

            #endif

            #if kernels > 8

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 7, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 8, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 9, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 11, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 12, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 13, 1, MPI_COMM_WORLD);

            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 14, 1, MPI_COMM_WORLD);
            
            MPI_Send(&x_n1[first_element_in_cur_offset_include], 
                        last_elem_in_cur_offset_not_include - first_element_in_cur_offset_include,
                                MPI_DOUBLE, 15, 1, MPI_COMM_WORLD);

            #endif

            MPI_Recv(&x_n1[0*(int)N/kernels],
                            1*(int)N/kernels - 0*(int)N/kernels, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

            #if kernels > 2

            MPI_Recv(&x_n1[1*(int)N/kernels],
                            2*(int)N/kernels - 1*(int)N/kernels, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[2*(int)N/kernels],
                            3*(int)N/kernels - 2*(int)N/kernels, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, &status);
            
            #endif

            #if kernels > 4

            MPI_Recv(&x_n1[3*(int)N/kernels],
                            4*(int)N/kernels - 3*(int)N/kernels, MPI_DOUBLE, 3, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[4*(int)N/kernels],
                            5*(int)N/kernels - 4*(int)N/kernels, MPI_DOUBLE, 4, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[5*(int)N/kernels],
                            6*(int)N/kernels - 5*(int)N/kernels, MPI_DOUBLE, 5, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[6*(int)N/kernels],
                            7*(int)N/kernels - 6*(int)N/kernels, MPI_DOUBLE, 6, 1, MPI_COMM_WORLD, &status);

            #endif

            #if kernels > 8

            MPI_Recv(&x_n1[7*(int)N/kernels],
                            8*(int)N/kernels - 7*(int)N/kernels, MPI_DOUBLE, 7, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[8*(int)N/kernels],
                             9*(int)N/kernels - 8*(int)N/kernels, MPI_DOUBLE, 8, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[9*(int)N/kernels],
                            10*(int)N/kernels - 9*(int)N/kernels, MPI_DOUBLE, 9, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[11*(int)N/kernels],
                            12*(int)N/kernels - 11*(int)N/kernels, MPI_DOUBLE, 11, 1, MPI_COMM_WORLD, &status);
            
            MPI_Recv(&x_n1[12*(int)N/kernels],
                            13*(int)N/kernels - 12*(int)N/kernels, MPI_DOUBLE, 12, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[13*(int)N/kernels],
                            14*(int)N/kernels - 13*(int)N/kernels, MPI_DOUBLE, 13, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[14*(int)N/kernels],
                            15*(int)N/kernels - 14*(int)N/kernels, MPI_DOUBLE, 14, 1, MPI_COMM_WORLD, &status);

            MPI_Recv(&x_n1[15*(int)N/kernels],
                            16*(int)N/kernels - 15*(int)N/kernels, MPI_DOUBLE, 15, 1, MPI_COMM_WORLD, &status);

            #endif
    

            memcpy(x_n, x_n1, N*sizeof(double));

       /*   printf("rank: %i\n", rank);
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

            MPI_Send(&sub_norm_Axn_b, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            
            MPI_Recv(&converged, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

        }




    #endif

        itter++;
        if(itter == 10000) break;

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