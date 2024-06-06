#include <iostream>
#include <cmath>
#include <clocale>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;

double dY(int i, int j, double* y); //функция для подсчёта разностей (y[i+1]-y[i]), где j - порядок разности
int factorial(int a); //функция для факториала
double divsub(int i, int n, double* y, double* x); //разделённые разности (i - номер игрика, n - номер j(порядок))
double polresult(double* y, double* x, double a, int n, int pr); //подсчёт полинома по неравномерной сетке
double xmult(double* x, double a, int k); //подсчёт (x-x[i]) для полинома
double xmult1(double* x, double a, int k, int j); //подсчёт (x-x[i]) для производной полинома
double xmult2(double* x, double a, int k, int j, int l); //подсчёт (x-x[i]) для второй производной полинома
double ravnom_pol(double* x, double* y, double a, double h, int pr, int n); //подсчёт полинома по равномерной сетке
double SKOrav(double* res, double* x, double* y, double h, int n, int m, int pr); //подсчёт СКО для равном.полинома
double SKOner(double* res, double* x, double* y, double h, int n, int m, int pr);//СКО для неравном. полинома


int main()
{

    double a, b, * x, * y, * res;
    int i, n, m, pr, h;
    char setka; //вид сетки

    ifstream input;
    ofstream FileOut("Output.txt");
    input.open("Input.txt");
    input >> pr;//порядок производной
    input >> n;//порядок полинома
    input >> setka;//равномен/нет
    input >> a;//границы
    input >> b;
    x = new double[n + 1];//узлы сетки
    y = new double[n + 1];//значения в узлах
    for (i = 0; i <= n; i++) {
        input >> x[i];
    }
    for (i = 0; i <= n; i++) {
        input >> y[i];
    }
    input >> m;//кол-во интервалов
    res = new double[m + 1];//узлы резултатирующей сетки
    for (i = 0; i <= m; i++) {
        input >> res[i];
    }

    FileOut << "Y:\n";
    for (i = 0; i <= n; i++) {
        FileOut << y[i] << " ";
    }
    FileOut << "\nX:\n";
    for (i = 0; i <= n; i++) {
        FileOut << x[i] << " ";
    }
    if (setka == 'r') {
        h = (b - a) / n;
        FileOut << "\nResultat dlya proizvodnoy  = " << pr << "\n";
        for (i = 0; i <= m; i++) {
            FileOut << "\nResult f(x[" << res[i] << "] )= " << ravnom_pol(x, y, res[i], h, pr, n);
        }

        FileOut << setw(7) << "\nCKO = " << scientific << SKOrav(res, x, y, h, n, m, pr);
    }
    else {
        h = (b - a) / n;
        FileOut << "\nResultat dlya proizvodnoy  = " << pr << "\n";
        for (i = 0; i <= m; i++) {
            FileOut << "\nf(x[" << res[i] << "] )= " << polresult(y, x, res[i], n, pr);
        }
        FileOut << setw(7) << "\nCKO = " << scientific << SKOner(res, x, y, h, n, m, pr);
    }



    return 0;
}
int factorial(int a)
/**
Здесь происходит рекурсия, т.е. выдаваемый результат зависит от предыдущего
**/
{
    if (a > 0) {
        return a * factorial(a - 1);
    }
    return 1;
}

double divsub(int i, int j, double* y, double* x)
/**
Здесь происходит рекурсия, т.е. выдаваемый результат зависит от предыдущего
**/
{
    if (j == i) {
        return y[i];
    }

    return (divsub(i + 1, j, y, x) - divsub(i, j - 1, y, x)) / (x[j] - x[i]);

}

double xmult(double* x, double a, int k)
{
    /**
    Скобки (а-x[i]) умножаются друг на друга k раз (здесь "x" заменён на "а" для удобства.
    **/
    int i;
    double s = 1;

    for (i = 0; i <= k; i++) {
        s *= (a - x[i]);
    }
    return s;

}
double xmult1(double* x, double a, int k, int j)
{
    /**
    То же самое,  что в предыдущей, но пропускается скобка, где i == j;
    **/
    int i;
    double s = 1;

    for (i = 0; i <= k; i++) {
        if (i == j)
            i++;
        if (i <= k)
            s *= (a - x[i]);
    }
    return s;
}

double xmult2(double* x, double a, int k, int j, int l)
{
    /**
    То же самое,  что в предыдущей, но пропускаются скобки, где i == j или i == l;
    **/

    int i;
    double s = 1;

    for (i = 0; i <= k; i++) {
        if (i != j && i != l)
            s *= (a - x[i]);
    }
    return s;
}

double dY(int i, int j, double* y)
/**
Здесь происходит рекурсия, т.е. выдаваемый результат зависит от предыдущего
Меняется порядок разности - j (уменьшается с каждой рекурсией)
**/
{
    if (j != 0) { return dY(i + 1, j - 1, y) - dY(i, j - 1, y); }
    return y[i];

}
double ravnom_pol(double* x, double* y, double a, double h, int pr, int n)
{
    int i, j, k, l;
    double s1 = 1, s2 = 0, s3 = 0, result = 0;
    if (pr == 0)
    {
        for (i = 0; i <= n; i++)
        {
            s1 = 1;
            for (j = 0; j < i; j++) {
                s1 *= (((a - x[0]) / h) - j); //П((q-x[i]))
            }
            result += s1 * (dY(0, i, y) / factorial(i));
        }
    }

    if (pr == 1)
    {
        for (i = 1; i <= n; i++) {
            for (j = 0; j < i; j++) {
                s1 = 1;
                for (k = 0; k < i; k++) {
                    if (k != j)
                        s1 *= (((a - x[0]) / h) - k); //П((q-x[i]))
                }
                s2 += s1;
            }
            result += s2 * (dY(0, i, y) / factorial(i));

            s2 = 0;
        }
        result /= h;
    }

    if (pr == 2)
    {
        for (i = 2; i <= n; i++) {

            for (j = 0; j < i; j++) {
                for (k = 0; k < i; k++) {
                    if (k != j) {
                        s1 = 1;
                        for (l = 0; l < i; l++) {
                            if (l != k && l != j)
                                s1 *= (((a - x[0]) / h) - l); // П((q-x[i]))
                        }
                    }
                    s2 += s1;
                }
                s3 += s2;
                s2 = 0;
            }
            result += (dY(0, i, y) / factorial(i)) * s3;
            s3 = 0;
        }
        result /= pow(h, 2);
    }
    return result;
}
double polresult(double* p, double* x, double a, int n, int pr)
{
    int j, i, k;
    double result = 0, s = 0, s1 = 0;

    if (pr == 0) {
        result += p[0];
        for (i = 1; i <= n; i++) {
            s1 = 0;
            s1 = xmult(x, a, i - 1); //то же самое, что П((x-x[i])); (здесь цикл происходит внутри функции xmult)

            result += divsub(0, i, p, x) * s1;
        }

    }
    else if (pr == 2) {
        for (i = 2; i <= n; i++) {
            for (j = 0; j < i; j++) {
                for (k = 0; k < i; k++) {
                    if (k != j)
                        s1 += xmult2(x, a, i - 1, k, j); //то же самое, что П((x-x[i])); (здесь цикл происходит внутри функции xmult)
                }
                s += s1;
                s1 = 0;
            }
            result += divsub(0, i, p, x) * s;
            s = 0;
        }
    }
    else if (pr == 1) {
        for (i = 1; i <= n; i++) {
            s = 0;
            for (j = 0; j < i; j++) {
                s += xmult1(x, a, i - 1, j); //то же самое, что П((x-x[i])); (здесь цикл происходит внутри функции xmult)
            }
            result += divsub(0, i, p, x) * s;
        }
    }
    return result;

}

double SKOrav(double* res, double* x, double* y, double h, int n, int m, int pr)
{

    double s = 0, t = 0, E;
    int i;
    if (pr == 0) {
        for (i = 0; i <= m; i++) {
            s = ((1 - pow(res[i], 4)) - ravnom_pol(x, y, res[i], h, pr, n)); // 1-x^4

            t += pow(s, 2);
        }
    }
    if (pr == 1) {
        for (i = 0; i <= m; i++) {
            s = ((-4) * pow(res[i], 3) - ravnom_pol(x, y, res[i], h, pr, n)); // -4x^3

            t += pow(s, 2);
        }
    }
    if (pr == 2) {
        for (i = 0; i <= m; i++) {
            s = ((-12 * pow(res[i], 2)) - ravnom_pol(x, y, res[i], h, pr, n)); // -12x^2

            t += pow(s, 2);
        }
    }

    E = sqrt(t);
    E = sqrt(E);
    return E / n;
}
double SKOner(double* res, double* x, double* y, double h, int n, int m, int pr)
{

    double s = 0, t = 0, E;
    int i;
    if (pr == 0) {
        for (i = 0; i <= m; i++) {
            s = ((1 - pow(res[i], 4)) - polresult(y, x, res[i], n, pr)); //f(x) = 1 - x^4

            t += pow(s, 2);
        }
    }
    if (pr == 1) {
        for (i = 0; i <= m; i++) {
            s = (((-4) * pow(res[i], 3)) - polresult(y, x, res[i], n, pr)); //f'(x) = -4x^3

            t += pow(s, 2);
        }
    }
    if (pr == 2) {
        for (i = 0; i <= m; i++) {
            s = ((-12 * pow(res[i], 2)) - polresult(y, x, res[i], n, pr)); //f''(x) = -12x^2

            t += pow(s, 2);
        }
    }

    E = sqrt(t);
    E = sqrt(E);
    return E / n;
}