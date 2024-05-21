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
float* multiplyByVector(float** A, float* V, int n)
{
    float* res = new float[n];
    for (int i = 0; i < n; i++)
    {
        res[i] = 0;
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res[i] += A[i][j] * V[j];
        }
    }
    return res;
}
float** multiplyMatrix(float** A, float** B, int an, int bn)
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
float** copyMatrix(float** A, int n)
{
    float** copy = new float* [n];
    initializeMatrix(copy, n,n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            copy[i][j] = A[i][j];
        }
    }
    return copy;
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
void printVector(float* V, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("%e ", V[i]);
    }
    printf("\n");
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

float* readVector(int n, FILE* f)
{
    float* V = new float[n];
    for (int i = 0; i < n; i++)
    {
        fscanf(f, "%f", &V[i]);
    }
    return V;
}

float** readMatrix(int n, FILE* f)
{
    float** M = new float* [n];
    initializeMatrix(M, n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fscanf(f, "%f", &M[i][j]);
        }
    }
    return M;
}
float** transpose(float** A, int n)
{
    float** res = new float* [n];
    initializeMatrix(res, n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res[j][i] = A[i][j];
        }
    }
    return res;
}
float* solveWithLU(float** A, float* b, int n)
{
    float** ACopy = copyMatrix(A, n);
    float* y = new float[n];
    float** L = new float* [n];
    float** U = new float* [n];
    initializeMatrix(L, n, n);
    initializeMatrix(U, n, n);
    float** P = calculateLUWithMainElement(A, L, U, n);
    y[0] = b[0] / L[0][0];
    float tmp;
    for (int i = 1; i < n; i++)
    {
        tmp = 0;
        for (int j = 0; j < i; j++)
        {
            tmp += L[i][j] * y[j];
        }
        y[i] = (b[i] - tmp)/L[i][i];
    }
    float* x = new float[n];
    for (int i = n - 1; i >= 0; i--)
    {
        tmp = 0;
        for (int j = i + 1; j < n; j++)
        {
            tmp += U[i][j] * x[j];
        }
        x[i] = (y[i] - tmp) / U[i][i];
    }
    float** Pt = transpose(P, n);
    float* res = multiplyByVector(Pt, x,n);
    float* test = multiplyByVector(ACopy, res, n);
    printf("A:\n");
    printMatrix(A, n);
    printf("b:\n");
    printVector(b, n);
    printf("L:\n");
    printMatrix(L, n);
    printf("U:\n");
    printMatrix(U, n);
    printf("P:\n");
    printMatrix(P, n);
    printf("y:\n");
    printVector(y, n);
    printf("x':\n");
    printVector(x, n);
    printf("x:\n");
    printVector(res, n);
    printf("A*x:\n");
    printVector(test, n);
    return res;
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
    float** A;
    A = readMatrix(n, f);
    float* b;
    b = readVector(n, f);
    solveWithLU(A,b, n);

}


