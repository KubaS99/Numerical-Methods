#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <float.h>


using namespace std;
char prim = ' ';

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
float rowSum(float** A, int n, int row)
{
    float sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += fabs(A[row][i]);
    }
    return sum;
}
float colSum(float** A, int n, int col)
{
    float sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += fabs(A[i][col]);
    }
    return sum;
}
float findRowMax(float** A, int n, int row)
{
    float max = -FLT_MAX;
    for (int i = 0; i < n; i++)
    {
        if (fabs(A[row][i]) > max)
        {
            max = fabs(A[row][i]);
        }
    }
    return max;
}
float findMax(float** A, int n, int row)
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
    initializeMatrix(copy, n, n);
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
void insertColumn(float** A, float* v, int index, int n)
{
    for (int i = 0; i < n; i++)
    {
        A[i][index] = v[i];
    }
}
float* getColumn(float** A, int col, int n)
{
    float* res = new float[n];
    for (int i = 0; i < n; i++)
    {
        res[i] = A[i][col];
    }
    return res;
}
float matrixNormInf(float** A, int n)
{
    float max = -FLT_MAX;
    for (int i = 0; i < n; i++)
    {
        float tmp = colSum(A, n, i);
        if (tmp > max)
        {
            max = tmp;
        }
    }
    return max;
}
float matrixNormZero(float** A, int n)
{
    float max = -FLT_MAX;
    for (int i = 0; i < n; i++)
    {
        float tmp = rowSum(A, n, i);
        if (tmp > max)
        {
            max = tmp;
        }
    }
    return max;
}
float vectorNormInf(float* v, int n)
{
    float max = -FLT_MAX;
    for (int i = 0; i < n; i++)
    {
        if (fabs(v[i]) > max)
        {
            max = fabs(v[i]);
        }
    }
    return max;
}
float vectorNormZero(float* v, int n)
{
    float sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += fabs(v[i]);
    }
    return sum;
}
void InverseMatrix(float** L, float** U,float** P,float** A,float* x, float* b, int n)
{
    float** Linv = new float* [n];
    float** Uinv = new float* [n];
    float** I = new float* [n];
    initializeMatrix(Linv, n, n);
    initializeMatrix(Uinv, n, n);
    initializeMatrix(I, n, n);
    identityMatrix(I, n);
    float tmp;
    float* e;
    float* y = new float[n];
    float* z = new float[n];
    for (int i = 0; i < n; i++)
    {
        z[i] = 0;
    }
    for (int k = 0; k < n; k++)
    {
        e = getColumn(I, k, n);
        for (int i = 0; i < n; i++)
        {
            tmp = 0;
            for (int j = 0; j < i; j++)
            {
                tmp += L[i][j] * y[j];
            }
            y[i] = (e[i] - tmp) / L[i][i];
        }
        for (int i = n - 1; i >= 0; i--)
        {
            tmp = 0;
            for (int j = i + 1; j < n; j++)
            {
                tmp += U[i][j] * z[j];
            }
            z[i] = (e[i] - tmp) / U[i][i];
        }
        insertColumn(Linv, y, k, n);
        insertColumn(Uinv, z, k, n);
    }

    printf("L^(-1): \n");
    printMatrix(Linv, n);
    printf("U^(-1): \n");
    printMatrix(Uinv, n);
    float** Pinv = transpose(P, n);
    printf("P^(-1): \n");
    printMatrix(Pinv, n);
    float** Ainv = multiplyMatrix(Pinv, Uinv, n, n);
    Ainv = multiplyMatrix(Ainv, Linv, n, n);
    printf("A^(-1) = P^(-1)*U^(-1)*L^(-1): \n");
    printMatrix(Ainv, n);
    float** test = multiplyMatrix(A, Ainv, n, n);
    printf("A*A^(-1): \n");
    printMatrix(test, n);
    printf("cond(A): \n");
    float Anorm = matrixNormZero(A, n);
    float AinvNorm = matrixNormZero(Ainv, n);
    printf("%e\n", Anorm * AinvNorm);
    printf("~cond(A): \n");
    printf("%e\n", Anorm * (vectorNormZero(x, n) / vectorNormZero(b, n)));
}
float* solveWithLU(float** A, float* b, int n)
{
    float** ACopy = copyMatrix(A, n);
    float* y = new float[n];
    float** L = new float* [n];
    float** U = new float* [n];
    float* bcopy = new float[n];
    for (int i = 0; i < n; i++)
    {
        bcopy[i] = b[i];
    }
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
        y[i] = (b[i] - tmp) / L[i][i];
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
    float* res = multiplyByVector(Pt, x, n);
    float* test = multiplyByVector(ACopy, res, n);
    printf("A%c:\n",prim);
    printMatrix(ACopy, n);
    printf("b%c:\n",prim);
    printVector(b, n);
    printf("L:\n");
    printMatrix(L, n);
    printf("U:\n");
    printMatrix(U, n);
    printf("P:\n");
    printMatrix(P, n);
    printf("x:\n");
    printVector(res, n);
    printf("A*x:\n");
    printVector(test, n);
    InverseMatrix(L, U, P, ACopy,x,bcopy, n);
    return res;
}
void ScaleLinEq(float** A,float** Aprim, float* b,float* bprim, int n)
{
    float max = 0;
    for (int i = 0; i < n; i++)
    {
        max = findRowMax(A, n, i);
        for (int j = 0; j < n; j++)
        {
            Aprim[i][j] = A[i][j]/max;
        }
        bprim[i] = b[i]/max;
    }
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
    float** Acopy = copyMatrix(A, n);
    float* bcopy = new float[n];
    for (int i = 0; i < n; i++)
    {
        bcopy[i] = b[i];
    }
    solveWithLU(A, b, n);
    float** Aprim = new float* [n];
    initializeMatrix(Aprim, n, n);
    float* bprim = new float[n];
    printf("====================================================================================================\n");
    prim = '\'';
    ScaleLinEq(Acopy, Aprim, bcopy, bprim, n);
    solveWithLU(Aprim, bprim, n);

}


