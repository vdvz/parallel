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
float calcAlphaN1(float* vectorR, float* vectorZ, float* matrixA, int N, float* vectorTmp) {
    double devidend = 0; //делимое
    int i, j;
    #pragma omp parallel for reduction(+:devidend)
    for (i = 0; i < N; i++) {
        devidend += vectorR[i] * vectorR[i];
    }
    #pragma omp parallel for
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            vectorTmp[i] += matrixA[i * N + j] * vectorZ[j];
        }
    }
    double divider = 0;
    
    #pragma omp parallel for reduction(+:divider)
    for (i = 0; i < N; i++) {
        divider += vectorTmp[i] * vectorZ[i];
    }
    return devidend / divider;
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
    #pragma omp parallel for
    for (i = 0; i < N; i++) {
        vectorX[i] += (vectorZ[i] * alphaN1);
    }
}
void calcVectorRNPlus1(float* vectorR, float alphaN1, float* matrixA, float* vectorZ, int N, float* vectorTmp) {
    int i, j;
    #pragma omp parallel for
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            vectorTmp[i] += matrixA[i * N + j] * vectorZ[j];
        }
    }
    #pragma omp parallel for
    for (i = 0; i < N; i++) {
        vectorR[i] -= (alphaN1 * vectorTmp[i]);
    }
}
float calcForBN1(float* vectorR, int N) {
    double res = 0;
    int i;
    #pragma omp parallel for reduction(+:res)
    for (i = 0; i < N; i++) {
        res += vectorR[i] * vectorR[i];
    }
    return res;
}
void calcVectorZNPlus1(float* vectorZ, float* vectorR, float bettaN1, int N) {
    int i;
    #pragma omp parallel for
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
float calcForFractionForCheckCondition(float* vector, int N) {
    double res = 0;
    int i;
    #pragma omp parallel for reduction(+:res)
    for (i = 0; i < N; i++) {
        res += (vector[i] * vector[i]);
    }
    return sqrt(res);
}
void printVectorXToFile(float* vectorX, FILE* fileOut, int N) {
    fwrite(vectorX, sizeof(float), N, fileOut);
}
void initVectorTmp(float* vectorTmp, int N) {
    int i;
    #pragma omp parallel for
    for (i = 0; i < N; i++) {
        vectorTmp[i] = 0;
    }
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
    float divider = calcForFractionForCheckCondition(vectorB, N);
    int counter = 0;
    float tmp = 0;
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    while ((tmp = calcForFractionForCheckCondition(vectorR, N) / divider) > epsilon) {
        float alphaN1 = calcAlphaN1(vectorR, vectorZ, matrixA, N, vectorTmp);
        initVectorTmp(vectorTmp, N);
        calcVectorXNPlus1(vectorX, alphaN1, vectorZ, N);
        float dividerForCalcBN1 = calcForBN1(vectorR, N);
        calcVectorRNPlus1(vectorR, alphaN1, matrixA, vectorZ, N, vectorTmp);
        initVectorTmp(vectorTmp, N);
        float bettaN1 = calcForBN1(vectorR, N) / dividerForCalcBN1;
        calcVectorZNPlus1(vectorZ, vectorR, bettaN1, N);
        counter++;
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