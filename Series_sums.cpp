#define _CRT_SECURE_NO_WARNINGS
#include <iostream>

void LiczWartosci(int N, FILE* fptr)
{
    double wynikDoklWPrzod=0,wynikDoklWTyl=0;
    float wynikPrzyblWPrzod=0, wynikPrzyblWTyl = 0;
    double wartoscWPrzod, wartoscWTyl;


    for (int n = 1; n <= N; n++)
    {
        wartoscWPrzod = 1.0 / n;
        wartoscWTyl = 1.0 / (N - (n-1));
        wynikDoklWPrzod += wartoscWPrzod;
        wynikPrzyblWPrzod += wartoscWPrzod;
        wynikDoklWTyl += wartoscWTyl;
        wynikPrzyblWTyl += wartoscWTyl;
    }
    fprintf(fptr, "%d\t%.15lf\t%.15lf\t%.15lf\t%.15lf\n",N, wynikPrzyblWPrzod, wynikPrzyblWTyl, wynikDoklWPrzod, wynikDoklWTyl);
}

int main()
{
    FILE* fptr = fopen("data.txt", "w");
    fprintf(fptr, "#N\tS(a)F\t\t\tS(b)F\t\t\tS(a)D\t\t\tS(b)D\n");
    for (int i = 1; i <= 300; i++)
    {
        LiczWartosci(i * 10000, fptr);
    }
}


