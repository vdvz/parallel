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