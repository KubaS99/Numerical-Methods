#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <cmath>

double** AlphaInv;



double Function(int x)
{
    return 4 + 3 * x + 2 * (x * x) + (x * x * x);
}

double derivative(int x1, int x2)
{
    return std::abs((Function(x2) - Function(x1)) / (x2 - x1));
}
void initializeMatrix(double** M, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        M[i] = new double[n];
    }
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            M[i][j] = 0;
        }
    }
}

double* multiplyByVector(double** A, double* v, int m, int n)
{
    double* res = new double[m];
    for (int i = 0; i < m; i++)
    {
        res[i] = 0;
    }
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res[i] += A[i][j] * v[j];
        }
    }
    return res;
}

double** multiplyMatrix(double** A, double** B, int m, int n)
{
    double** res = new double*[n];
    initializeMatrix(res, n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < m; k++)
            {
                res[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return res;
}

double** transpose(double** A, int m, int n)
{
    double** res = new double* [n];
    initializeMatrix(res, n, m);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res[j][i] = A[i][j];
        }
    }
    return res;
}
void printVector(double* V, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("%f ", V[i]);
    }
    printf("\n\n");
}

void printMatrix(double** M, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%f ", M[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

double** copyMatrix(double** A, int m, int n)
{
    double** copy = new double* [m];
    initializeMatrix(copy, m, n);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            copy[i][j] = A[i][j];
        }
    }
    return copy;
}

double Base(int index, int x)
{
    return pow(x,index);
}

double ChiSquare(double* y, double* a, double* noise,int n, int m)
{
    double sum = 0;
    double tmp = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            tmp += a[j] * Base(j, i + 1);
        }
        sum += pow((y[i] - tmp) / noise[i],2);
        tmp = 0;
    }
    return sum;
}

void testFill(double** A, int m, int n, double* noise)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] = Base(j, i);
        }
    }
}
void fillMatrixA(double** A, int m, int n,double* noise)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] = Base(j, i + 1) / fabs(noise[i]);
        }
    }
}
int findMax(double** A, int n, int row)
{
    double max = -DBL_MAX;
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
void swapCols(double** A, int n, int fcol, int scol)
{
    double* tmp = new double[n];
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
void identityMatrix(double** A, int n)
{
    for (int i = 0; i < n; i++)
    {
        A[i][i] = 1;
    }
}

double rowSum(double** A, int n, int row)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += fabs(A[row][i]);
    }
    return sum;
}
void insertColumn(double** A, double* v, int index, int n)
{
    for (int i = 0; i < n; i++)
    {
        A[i][index] = v[i];
    }
}

double* getColumn(double** A, int col, int n)
{
    double* res = new double[n];
    for (int i = 0; i < n; i++)
    {
        res[i] = A[i][col];
    }
    return res;
}

double vectorNormZero(double* v, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += fabs(v[i]);
    }
    return sum;
}

double matrixNormZero(double** A, int n)
{
    double max = -DBL_MAX;
    for (int i = 0; i < n; i++)
    {
        double tmp = rowSum(A, n, i);
        if (tmp > max)
        {
            max = tmp;
        }
    }
    return max;
}

void InverseMatrix(double** L, double** U, double** P, double** A, double* x, double* b, int n)
{
    double** Linv = new double* [n];
    double** Uinv = new double* [n];
    double** I = new double * [n];
    initializeMatrix(Linv, n, n);
    initializeMatrix(Uinv, n, n);
    initializeMatrix(I, n, n);
    identityMatrix(I, n);
    double tmp;
    double* e;
    double* y = new double[n];
    double* z = new double[n];
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


    double** Pinv = transpose(P, n, n);

    AlphaInv = multiplyMatrix(Pinv, Uinv, n, n);
    AlphaInv = multiplyMatrix(AlphaInv, Linv, n, n);

}

double** calculateLUWithMainElement(double** M, double** L, double** U, int n)
{
    double** Permutation = new double* [n];
    initializeMatrix(Permutation, n, n);
    identityMatrix(Permutation, n);

    identityMatrix(L, n);
    double tmp;
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
    /*printf("A=L*U? \n");
    printMatrix(multiplyMatrix(L, U, n, n),n,n);*/
    return Permutation;
}

double* solveWithLU(double** A, double* b, int n)
{
    double** ACopy = copyMatrix(A, n,n);
    double* y = new double[n];
    double** L = new double* [n];
    double** U = new double * [n];
    initializeMatrix(L, n, n);
    initializeMatrix(U, n, n);
    double** P = calculateLUWithMainElement(ACopy, L, U, n);

    y[0] = b[0] / L[0][0];
    double tmp;
    for (int i = 1; i < n; i++)
    {
        tmp = 0;
        for (int j = 0; j < i; j++)
        {
            tmp += L[i][j] * y[j];
        }
        y[i] = (b[i] - tmp) / L[i][i];
    }
    double* x = new double[n];
    for (int i = n - 1; i >= 0; i--)
    {
        tmp = 0;
        for (int j = i + 1; j < n; j++)
        {
            tmp += U[i][j] * x[j];
        }
        x[i] = (y[i] - tmp) / U[i][i];
    }
    double** Pt = transpose(P, n,n);
    double* res = multiplyByVector(Pt, x, n,n);
    //double* test = multiplyByVector(ACopy, res, n,n);
    InverseMatrix(L, U, P, A, x, b, n);
    return res;
}



int main()
{
    int m = 20;
    int n = 4;
    double* data = new double[m];
    for (int i = 0; i < m; i++)
    {
        data[i] = Function(i+1);
    }
    double noise[20] = { -0.282474,0.573539,0.061420,-1.888663,-0.411312,0.852528,0.841401,-1.072427,-0.586583,2.835644,-0.684718,-0.160293,-0.167703,-0.365004,-1.124204,0.189771,-0.750554,-0.725392,-1.591409,-0.971162 };
    double noise2[20];
    for (int i = 0; i < 20; i++)
    {
        noise2[i] = 1;
    }
    double* noisedData = new double[m];
    for (int i = 0; i < m; i++)
    {
        noisedData[i] = data[i] + noise[i];
    }
    printf("================= SZUM ==================\n");
    double** A1 = new double* [m];
    initializeMatrix(A1, m,n);
    fillMatrixA(A1, m, n, noise);
    printf("A:\n");
    printMatrix(A1, m, n);
    double** Atr1 = transpose(A1, m, n);
    printf("Atr: \n");
    printMatrix(Atr1, n, m);
    double** Alpha1 = multiplyMatrix(Atr1, A1, m,n);
    double* b1 = new double[m];
    for (int i = 0; i < m; i++)
    {
        b1[i] = noisedData[i] / fabs(noise[i]);
    }

    printf("b: \n");
    printVector(b1, m);

    double* Beta = multiplyByVector(Atr1, b1, n, m);

    printf("Alpha:\n");
    printMatrix(Alpha1, n, n);
    printf("Beta:\n");
    printVector(Beta, n);

    double** Alpha1copy = copyMatrix(Alpha1, n, n);
    double* res = solveWithLU(Alpha1, Beta, n);
    printf("a: \n");
    printVector(res, n);

    printf("Alpha*a: \n");
    printVector(multiplyByVector(Alpha1, res, n, n),n);

    printf("Cond(Alpha): %e\n\n", matrixNormZero(Alpha1copy,n)*matrixNormZero(AlphaInv,n));

    printf("Chi^2: %e\n\n", ChiSquare(noisedData,res,noise,m,n));

    printf("Macierz kowariancji: \n");
    printMatrix(AlphaInv, n, n);


    printf("================= BEZ SZUMU ==================\n");
    double** A2 = new double* [m];
    initializeMatrix(A2, m, n);
    fillMatrixA(A2, m, n, noise2);
    printf("A:\n");
    printMatrix(A2, m, n);
    double** Atr2 = transpose(A2, m, n);
    printf("Atr: \n");
    printMatrix(Atr2, n, m);


    double** Alpha2 = multiplyMatrix(Atr2, A2, m, n);
    double* b2 = new double[m];
    for (int i = 0; i < m; i++)
    {
        b2[i] = data[i]/noise2[i];
    }

    printf("b: \n");
    printVector(b2, m);

    double* Beta2 = multiplyByVector(Atr2, b2, n, m);

    printf("Alpha:\n");
    printMatrix(Alpha2, n, n);
    printf("Beta:\n");
    printVector(Beta2, n);

    double** Alpha2copy = copyMatrix(Alpha2, n, n);
    double* res2 = solveWithLU(Alpha2, Beta2, n);
    printf("a: \n");
    printVector(res2, n);

    printf("Alpha*a: \n");
    printVector(multiplyByVector(Alpha2, res2, n, n), n);

    printf("Cond(Alpha): %e\n\n", matrixNormZero(Alpha2copy, n) * matrixNormZero(AlphaInv, n));

    printf("Chi^2: %e\n\n", ChiSquare(data, res2, noise2, m, n));

    printf("Macierz kowariancji: \n");
    printMatrix(AlphaInv, n, n);
}

