#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void getMatrixFromFile(float* matrixA, int N, FILE* fileIn) {
    float* buffer = (float*)malloc(N * N * sizeof(float));
    fread(buffer, sizeof(float), N * N, fileIn);
    int n = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrixA[i * N + j] = buffer[n++];
        }
    }
    free(buffer);
}

void initVectorX0(float* vectorX0, int N) {
    for (int i = 0; i < N; i++) {
        vectorX0[i] = 1;
    }
}

void calcVectorR0(float* vectorR0, float* vectorB, float* vectorX0, float* matrixA, int N, float* vectorTmp) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            vectorTmp[i] += matrixA[i * N + j] * vectorX0[j];
        }
    }
    for (int i = 0; i < N; i++) {
        vectorR0[i] = vectorB[i] - vectorTmp[i];
    }
}
void copyVectorR0ToZ0(float* vectorR0, float* vectorZ0, int N) {
    int i;
    #pragma omp parallel for
    for (i = 0; i < N; i++) {
        vectorZ0[i] = vectorR0[i];
    }
}
void calcVectorXNPlus1(float* vectorX, float alphaN1, float* vectorZ, int N) {
    int i;
    #pragma omp for
    for (i = 0; i < N; i++) {
        vectorX[i] += (vectorZ[i] * alphaN1);
    }
}
void calcVectorRNPlus1(float* vectorR, float alphaN1, float* matrixA, float* vectorZ, int N, float* vectorTmp) {
    int i, j;
    #pragma omp for
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            vectorTmp[i] += matrixA[i * N + j] * vectorZ[j];
        }
    }
    #pragma omp for
    for (i = 0; i < N; i++) {
        vectorR[i] -= (alphaN1 * vectorTmp[i]);
    }
}
void calcVectorZNPlus1(float* vectorZ, float* vectorR, float bettaN1, int N) {
    int i;
    #pragma omp for
    for (i = 0; i < N; i++) {
        vectorZ[i] = bettaN1 * vectorZ[i] + vectorR[i];
    }
}
double check(float* response, float* vectorX, float* vectorTmp, int N) {
    double tmpDev = 0, tmpDiv = 0;
    int i;
    for (i = 0; i < N; i++) {
        vectorTmp[i] = fabs(response[i]) - fabs(vectorX[i]);
    }
    for (i = 0; i < N; i++) {
        tmpDev += vectorTmp[i] * vectorTmp[i];
        tmpDiv += response[i] * response[i];
    }
    double eps = sqrt(tmpDev) / sqrt(tmpDiv);
    return eps;
}
void printVectorXToFile(float* vectorX, FILE* fileOut, int N) {
    fwrite(vectorX, sizeof(float), N, fileOut);
}
void initVectorTmp(float* vectorTmp, int N) {
    int i;
    #pragma omp for
    for (i = 0; i < N; i++) {
        vectorTmp[i] = 0;
    }
}
float calcForFractionForCheckCondition(float* vectorB, int N) {
    double tmp = 0;
    int i;
    for (i = 0; i < N; i++) {
        tmp += (vectorB[i] * vectorB[i]);
    }
    return sqrt(tmp);
}
int main(int agrc, char* argv[]) {
    FILE* fileInForMatrix = fopen(argv[1], "rb");
    fseek(fileInForMatrix, 0, SEEK_END);
    int64_t fileSizeForMatrix = 0;
    fileSizeForMatrix = ftell(fileInForMatrix); //узнаю размер файла
    fseek(fileInForMatrix, 0, SEEK_SET);
    int N = sqrt(fileSizeForMatrix / sizeof(float)); //Размерность матрицы
    float *matrixA = (float*)malloc(N * N * sizeof(float));

    getMatrixFromFile(matrixA, N, fileInForMatrix);

    fclose(fileInForMatrix);

    FILE* fileInForVector = fopen(argv[2], "rb");
    float* vectorB = (float*)malloc(N * sizeof(float)); //вектор правых частей
    fread(vectorB, sizeof(float), N, fileInForVector);
    fclose(fileInForVector);

    float* vectorX = (float*)malloc(N * sizeof(float)); //произвольное начальное приближение решения
    initVectorX0(vectorX, N);

    float* vectorR = (float*)malloc(N * sizeof(float));
    float* vectorZ = (float*)malloc(N * sizeof(float));

    float* vectorTmp = (float*)malloc(N * sizeof(float));
    initVectorTmp(vectorTmp, N);
    calcVectorR0(vectorR, vectorB, vectorX, matrixA, N, vectorTmp);
    initVectorTmp(vectorTmp, N);
    copyVectorR0ToZ0(vectorR, vectorZ, N);

    float epsilon = 0.000000001;
    float dividerMain = calcForFractionForCheckCondition(vectorB, N);
    int counter = 0;
    struct timespec start, end;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    double devidend = 0; //делимое
    double divider = 0;
    double tmp = 0;
    float alphaN1, bettaN1;
    int i, j;

    #pragma omp parallel
    {
    
        //calcForFractionForCheckCondition
        #pragma omp for reduction(+:tmp)
        for (i = 0; i < N; i++) {
            tmp += (vectorR[i] * vectorR[i]);
        }
        //
        while (sqrt(tmp) /
        dividerMain > epsilon) {
            //calcAlpha
            #pragma omp for reduction(+:devidend)
            for (i = 0; i < N; i++) {
            devidend += vectorR[i] * vectorR[i];
            }
            #pragma omp for
            for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
            vectorTmp[i] += matrixA[i * N + j] * vectorZ[j];
            }
            }
            #pragma omp for reduction(+:divider)
            for (i = 0; i < N; i++) {
                divider += vectorTmp[i] * vectorZ[i];
            }
            #pragma omp single
                {
                    alphaN1 = devidend / divider;
                    devidend = 0; divider = 0; tmp = 0;
                }
            initVectorTmp(vectorTmp, N);
            //
            calcVectorXNPlus1(vectorX, alphaN1, vectorZ, N);
            //calcForBN1
            #pragma omp for reduction(+:divider)
            for (i = 0; i < N; i++) {
                divider += vectorR[i] * vectorR[i];
            }
            //
            calcVectorRNPlus1(vectorR, alphaN1, matrixA, vectorZ, N, vectorTmp);
            initVectorTmp(vectorTmp, N);
            //calc bettaN1
            #pragma omp for reduction(+:devidend)
            for (i = 0; i < N; i++) {
                devidend += vectorR[i] * vectorR[i];
            }
            #pragma omp single
                {
                    bettaN1 = devidend / divider;
                }
            //
            calcVectorZNPlus1(vectorZ, vectorR, bettaN1, N);
            #pragma omp for reduction(+:tmp)
            for (i = 0; i < N; i++) {
                tmp += (vectorR[i] * vectorR[i]);
            }
            #pragma omp single
            {
            devidend = 0; divider = 0;
            counter++;
            }
        }
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    double time = end.tv_sec - start.tv_sec + 0.000000001 * (end.tv_nsec - start.tv_nsec);

    printf("count: %d, time: %lf\n", counter, time);

    FILE* fileInForX = fopen(argv[4], "rb");
    float* response = (float*)malloc(N * sizeof(float));
    fread(response, sizeof(float), N, fileInForX);
    fclose(fileInForX);

    initVectorTmp(vectorTmp, N);
    double flag = check(response, vectorX, vectorTmp, N);
    printf("epsilon = %e\n\n", flag);

    FILE* fileOut = fopen(argv[3], "wb");
    printVectorXToFile(vectorX, fileOut, N);
    fclose(fileOut);

    free(matrixA);
    free(vectorX);
    free(vectorR);
    free(vectorZ);
    free(vectorTmp);
    free(response);
    return 0;
}