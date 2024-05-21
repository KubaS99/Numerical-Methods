#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <float.h>

using namespace std;

void cleanMatrix(float** A, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] = 0;
        }
    }
}
void initializeMatrix(float** M, int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        M[i] = new float[m];
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            M[i][j] = 0;
        }
    }
}
int findMax(float** A, int n, int row)
{
   float max = -FLT_MAX;
   int col = 0;
   for (int i = 0; i < n; i++)
   {
       if (fabs(A[row][i]) > max)
       {
           col = i;
           max = fabs(A[row][i]);
       }
   }
   return col;
}
void swapCols(float** A, int n, int fcol, int scol)
{
    float* tmp = new float[n];
    for (int i = 0; i < n; i++)
    {
        tmp[i] = A[i][fcol];
        A[i][fcol] = A[i][scol];
    }
    for (int i = 0; i < n; i++)
    {
        A[i][scol] = tmp[i];
    }
}
void identityMatrix(float** A, int n)
{
    for (int i = 0; i < n; i++)
    {
        A[i][i] = 1;
    }
}

float** multiplyMatrix(float** A, float** B,int an, int bn)
{
    float** Result = new float* [an];
    initializeMatrix(Result, an, bn);
    for (int i = 0; i < an; i++)
    {
        for (int j = 0; j < bn; j++)
        {
            for (int k = 0; k < an; k++)
            {
                Result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return Result;
}

float** calculateLUWithMainElement(float** M, float** L, float** U, int n)
{
    float** Permutation = new float* [n];
    initializeMatrix(Permutation, n, n);
    identityMatrix(Permutation, n);

    identityMatrix(L, n);

    float tmp;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            tmp = 0;
            for (int k = 0; k < i; k++)
            {
                tmp += L[i][k] * U[k][j];
            }
            U[i][j] = M[i][j] - tmp;
        }
        int col = findMax(U, n, i);
        swapCols(M, n, i, col);
        swapCols(U, n, i, col);
        swapCols(Permutation, n, i, col);
        for (int j = i + 1; j < n; j++)
        {
            tmp = 0;
            for (int k = 0; k < i; k++)
            {
                tmp += L[j][k] * U[k][i];
            }
            L[j][i] = (M[j][i] - tmp) / U[i][i];
        }
    }
    return Permutation;
}

void calculateLU(float** M, float** L, float** U, int n)
{
    identityMatrix(L, n);
    float tmp;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            tmp = 0;
            for (int k = 0; k < i; k++)
            {
                tmp += L[i][k] * U[k][j];
            }
            U[i][j] = M[i][j] - tmp;
        }
        for (int j = i + 1; j < n; j++)
        {
            tmp = 0;
            for (int k = 0; k < i; k++)
            {
                tmp += L[j][k] * U[k][i];
            }
            L[j][i] = (M[j][i] - tmp) / U[i][i];
        }
    }
}

void printMatrix(float** M, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%e ", M[i][j]);
        }
        printf("\n");
    }
}


float** readMatrix(int n, FILE* f)
{
    float** M = new float* [n];
    initializeMatrix(M, n,n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fscanf(f, "%f", &M[i][j]);
        }
    }
    fclose(f);
    return M;
}


int main(int argc, char* argv[])
{
    int n;
    if (argc == 1)
    {
        printf("You need to pass a data file as program argument\n");
        return 1;
    }
    FILE* f;
    if ((f = fopen(argv[1], "r")) == NULL)
    {
        printf("Couldn't open %s\n", argv[1]);
        return 1;
    }
    fscanf(f, "%d", &n);

    float** M;
    float **L = new float* [n];
    float **U = new float* [n];
    initializeMatrix(L, n,n);
    initializeMatrix(U, n,n);

    M = readMatrix(n, f);
    printf("------Metoda Doolittle'a------\n");
    calculateLU(M, L, U, n);
    printf("A:\n");
    printMatrix(M, n);
    printf("\n\n");
    printf("L:\n");
    printMatrix(L, n);
    printf("\n\n");
    printf("U:\n");
    printMatrix(U, n);
    printf("\n\n");
    float** Res;
    Res = multiplyMatrix(L, U, n, n);
    printf("A = L * U\n");
    printMatrix(Res, n);
    printf("\n\n");

    printf("------Czesciowy wybor el. glownego - metoda Doolittle'a------\n");
    cleanMatrix(L, n);
    cleanMatrix(U, n);
    printf("A:\n");
    printMatrix(M, n);
    printf("\n\n");
    float** Permutation = calculateLUWithMainElement(M, L, U, n);
    printf("L:\n");
    printMatrix(L, n);
    printf("\n\n");
    printf("U:\n");
    printMatrix(U, n);
    printf("\n\n");
    printf("Permutation:\n");
    printMatrix(Permutation, n);
    printf("\n\n");
    Res = multiplyMatrix(L, U, n, n);
    printf("A = L * U:\n");
    printMatrix(Res, n);
    printf("\n\n");
    Res = multiplyMatrix(Res, Permutation, n, n);
    printf("A = L * U * P\n");
    printMatrix(Res, n);
}


