#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cmath>

double NajwiekszaD()
{
    double wynik = 0;
    double tmp = 0.5;
    while (wynik + tmp < 1)
    {
        wynik += tmp;
        tmp /= 2.0;
    }
    double tmpMax = wynik;
    while ((tmpMax < INFINITY))
    {
        wynik = tmpMax;
        tmpMax *= 2.0;
    }
    return wynik;
}

float NajwiekszaF()
{
    float wynik = 0;
    float tmp = 0.5;
    while (wynik + tmp < 1)
    {
        wynik += tmp;
        tmp /= 2.0;
    }
    float tmpMax = wynik;
    while ((tmpMax < INFINITY))
    {
        wynik = tmpMax;
        tmpMax *= 2.0;
    }   
    return wynik;
}
double ZnormalizowanaD(double e)
{
    double wynik = (1.0 + e) / 2.0;
    double sprawdzenie = 1.0 + e;
    while (wynik * 2.0 == sprawdzenie)
    {
        sprawdzenie = wynik;
        wynik /= 2.0;
    }
    return sprawdzenie;
}

float ZnormalizowanaF(float e)
{
    float wynik = (1.0 + e)/2.0;
    float sprawdzenie = 1.0 + e;
    while (wynik * 2.0 == sprawdzenie)
    {
        sprawdzenie = wynik;
        wynik /= 2.0;      
    }
    return sprawdzenie;
}
double ZdenormalizowanaD()
{
    double wynik = 1;
    double dzielnik = 0.5;
    while (1.0*dzielnik > 0)
    {
        wynik = 1.0 * dzielnik;
        dzielnik /= 2.0;
    }
    return wynik;
}


float ZdenormalizowanaF()
{
    float wynik=1;
    float dzielnik = 0.5;
    while (1.0*dzielnik > 0)
    {
        wynik = 1.0*dzielnik;
        dzielnik /= 2.0;
    }
    return wynik;
}

float WartoscSprawdzajacaF()
{
    return (float)fabs(3 * ((float)4.0 / (float)3.0 - (float)1.0) - 1);
}

double WartoscSprawdzajacaD()
{
    return fabs(3 * (4.0 / 3.0 - 1.0) - 1);
}



void EpsilonF()
{
    printf("############# Float #############\n");
    printf("#\te\t\t1+e\n");
    float wartosc;
    float tmp=1;
    int n = 1;
    while ((float)1.0 + tmp > (float)1.0)
    {
        wartosc = tmp;
        printf("%d\t%e\t%.7f\n", n, tmp, (float)1.0 + tmp);
        tmp /= 2.0;
        n++;
    }
    printf("Ilosc bitow mantysy: %d\n", n - 1);
    printf("Obliczony e mach: %e,\t\t\t wg standardu IEEE 754: 1.19209e-07\n", wartosc);
    printf("Zaokraglenie przez obcinanie: %e\n", powf(2, 1 - 24));
    printf("Zaokraglenie do najblizszej: %e\n", (powf(2, 1 - 24)) / 2.0);
    printf("Wartosc ze wzoru |3(4/3-1)-1|: %e\n", WartoscSprawdzajacaF());
    printf("Najmniejsza zdenormalizowana: %e, \twg standardu IEEE 754: 1.4013e-45\n", ZdenormalizowanaF());
    printf("Najmniejsza znormalizowana: %e, \twg standardu IEEE 754: 1.17549e-38\n", ZnormalizowanaF(wartosc));
    printf("Najwieksza mozliwa: %e, \t\twg standardu IEEE 754: 3.40282e+38\n", NajwiekszaF());
}

void EpsilonD()
{
    printf("\n\n############# Double #############\n");
    printf("#\te\t\t1+e\n");
    double wartosc;
    double tmp = 1;
    int n = 1;
    while ((double)1.0 + (double)tmp > 1.0)
    {
        wartosc = tmp;
        printf("%d\t%e\t%.16f\n",n,tmp,(double)1.0+tmp);
        tmp /= 2.0;
        n++;
    }
    printf("Ilosc bitow mantysy: %d\n", n-1);
    printf("Obliczony e mach: %e, \t\twg standardu IEEE 754: 2.22045e-16\n", wartosc);
    printf("Zaokraglenie przez obcinanie: %e\n", pow(2, 1 - 53));
    printf("Zaokraglenie do najblizszej: %e\n", (pow(2, 1 - 53))/2.0);
    printf("Wartosc ze wzoru |3(4/3-1)-1|: %e\n", WartoscSprawdzajacaD());
    printf("Najmniejsza zdenormalizowana: %e, \twg standardu IEEE 754: 4.94066e-324\n", ZdenormalizowanaD());
    printf("Najmniejsza znormalizowana: %e, \twg standardu IEEE 754: 2.22507e-308\n", ZnormalizowanaD(wartosc));
    printf("Najwieksza mozliwa: %e, \t\twg standardu IEEE 754: 1.79769e+308\n", NajwiekszaD());
}





int main()
{
    EpsilonF();
    EpsilonD();
}


