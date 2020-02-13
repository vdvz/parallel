#include <stdio.h>
#include <math.h>
#include <string.h> 

int main(){
    double tau = 0.01;
    const double epsilon = 0.00001;
    int N = 5;
    double * A = (double *)malloc(sizeof(double) * N * N);
    double * x_n = (double *)malloc(sizeof(double) * N);
    double * b = (double *)malloc(sizeof(double) * N);
    
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

    double norm_Axn_b = 0;
    double norm_b = 0; 
    for (int i = 0; i<N; i++){
       norm_b += b[i]*b[i];
    }
    norm_b = sqrt(norm_b);
   
    char converged = 0;
    int itter = 0;
    double * x_n1 =(double *)malloc(sizeof(double) * N);
    while(!converged){
      
        for (int i = 0; i<N; i++){
            double sub = 0;
            for(int j = 0; j<N; j++){
                sub += *( A + i*N + j) * x_n[j];
            }
            x_n1[i] = x_n[i] - tau*(sub - b[i]);
        }

        memcpy(x_n, x_n1, N*sizeof(double));

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

    free(x_n);
    free(x_n1);
    free(b);
    free(A);

    return 0;
}