#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>

float calculateX(int sign, float m, int e, int bias)
{
    return powf(-1, float(sign)) * m * powf(2.0, float(e - bias));
}
float calculateMantissa(int mantissa[23])
{
    float sum = 0;
    float tmp = 1.0 / 2.0;
    for (int i = 0; i < 23; i++)
    {
        sum += mantissa[i] * tmp;
        tmp /= 2.0;
    }
    return sum;
}
void getMantissa(int bytes[4][8], int mantissa[23])
{
    for (int i = 1; i < 8; i++)
    {
        mantissa[i - 1] = bytes[1][i];
    }
    for (int i = 0; i < 8; i++)
    {
        mantissa[i + 7] = bytes[2][i];
    }
    for (int i = 0; i < 8; i++)
    {
        mantissa[i + 15] = bytes[3][i];
    }
    printf("Mantissa: .");
    for (int i = 0; i < 23; i++)
    {
        printf("%d", mantissa[i]);
    }
    printf("\n");
}
int calculateExponent(int exponent[8])
{
    int sum = 0;
    for(int i = 7; i >= 0; i--)
    {
        sum += exponent[i] * pow(2, abs(i-7));
    }
    return sum;
}
void getExponent(int bytes[4][8], int exponent[8])
{
    for (int i = 1; i < 8; i++)
    {
        exponent[i-1] = bytes[0][i];
    }
    exponent[7] = bytes[1][0];

    printf("Exponent: ");
    for (int i = 0; i < 8; i++)
    {
        printf("%d", exponent[i]);
    }
    printf(".\n");
}
void printLittleEndian(int bytes[4][8])
{
    printf("Binary (little endian): ");
    for (int i = 3; i >= 0; i--)
    {
        for (int j = 0; j < 8; j++)
        {
            printf("%d", bytes[i][j]);
        }
        printf(" ");
    }
    printf("\n");
}
void printBigEndian(int bytes[4][8])
{
    printf("Binary (big endian): ");
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            printf("%d", bytes[i][j]);
        }
        printf(" ");
    }
    printf("\n");
}
void get_bytes(float f, int bytes[4][8])
{
    union {
        float f;
        uint32_t u;
    } fu = { f = f };
    int i = sizeof f * CHAR_BIT;


    while (i--)
    {
        bytes[3-(i / 8)][7-(i % 8)] = (fu.u >> i) & 0x1;
    }
}

int main()
{
    int bytes[4][8];
    float tmp;
    printf("Enter value: ");
    scanf("%f", &tmp);
    printf("Entered value: %e\n", tmp);
    get_bytes(tmp, bytes);
    printLittleEndian(bytes);
    printBigEndian(bytes);
    int sign = bytes[0][0];
    printf("Sign: %d\n", sign);
    int exponent[8];
    for (int i = 0; i < 8; i++)
    {
        exponent[i] = 0;
    }
    getExponent(bytes, exponent);
    int expVal = calculateExponent(exponent);
    printf("Exponent (e): %d\n", expVal);
    int bias = expVal == 0 ? 126 : 127;
    printf("Bias: %d\n", bias);
    int f0 = expVal == 0 ? 0 : 1;
    printf("f0: %d\n", f0);
    int mantissa[23];
    getMantissa(bytes, mantissa);
    float mantVal = calculateMantissa(mantissa);
    mantVal += float(f0);
    printf("Mantissa (m): %f\n", mantVal);
    printf("Decimal value: %.20f\n", calculateX(sign, mantVal, expVal, bias));
    return 0;
}


