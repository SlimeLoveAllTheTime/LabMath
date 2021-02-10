#include <iostream>
#include "lagrange.h"
#include "SPLINES.H"
#include "quanc8.h"

// массивы в которых будут храниться:
// точки по которым ищем значения функции и узлы для построения полиномов
double arrayX[21], *x = arrayX;

// значения функции в точках, которые хранятся в массиве выше
double arrayF[21], *f = arrayF;

// массив, в котором хранятся коэффициенты b (сплайн интерполяции)
double arrayB[21], *b = arrayB;

// массив, в котором хранятся коэффициенты c (сплайн интерполяции)
double arrayC[21], *c = arrayC;

// массив, в котором хранятся коэффициенты d (сплайн интерполяции)
double arrayD[21], *d = arrayD;


// функция, которая вычисляет k точек, по которым строятся полиномы
double X(int k) {
    return -1 + 0.1 * k;
}

// функция, которая считает значение заданной функции
double F(double t) {
    return (1/(1 + 25 * pow(t, 2)));
}

// функция-обертка для нахождения значения полинома Лагранжа с помощью функции quanc8
double lagrangeForQuanc(double t) {
    return lagrange(20, x, f, t);
}

// функция-обертка для нахождения значения сплайн-полинома с помощью функции quanc8
double splineForQuanc(double t) {
    double *p = &t;
    return seval(20, p, x, f, b, c, d);
}



int main() {

    // заполняем массивы arrayX и arrayF
    for (int i = 0; i < 21; i++) arrayX[i] = X(i);
    for (int i = 0; i < 21; i++) arrayF[i] = F(arrayX[i]);

    // интерполируем сплайнами
    spline(20, x, f, b, c, d);

    for (int i = 0; i < 21; i++) {

        // интерполируем полиномом Лагранжа и ищем значение функции в точках arrayX
        double resLagrange = lagrange(20, x, f, arrayX[i]);

        // рещаем систему, состоящую из полиномов третьей степени
        double spline = seval(20, &arrayX[i], x, f, b, c, d);

        // выводим результаты в консоль
        if (i == 0) printf("\n x                 F(x)              L(x)             S(x)\n");

        if (i < 9) printf("%f          %f          %f          %f\n", arrayX[i], arrayF[i], resLagrange, spline);
        if (i > 9) printf(" %f          %f          %f          %f\n", arrayX[i], arrayF[i], resLagrange, spline);
    }

    // абсолютная и относительная погрешность
    double abserr = 1.0e-14, relerr = 0.0;

    // реальная погрешность
    double errest, *Errest = &errest;

    // индикатор надежности
    double flag, *Flag = &flag;

    // количество внутренних вычислений функции
    int nofun, *Nofan = &nofun;

    // результат нахождения интеграла
    double res, *result = &res;

    // использование функции, которая находит значения функции в точках x
    double (*FUN1) (double x); FUN1 = F;

    // использование функции-обертки интерполирования полиномом Лагранжа
    double (*FUN2) (double x); FUN2 = lagrangeForQuanc;

    // использование функции-обертки интерполирования сплайнами
    double (*FUN3) (double x); FUN3 = splineForQuanc;

    printf("\nFunction    Integration result    Errest      NoFan     Flag\n");

    // нахождение и вывод интеграла от функции
    quanc8(FUN1, -1.0, 1.0, abserr, relerr, result, Errest, Nofan, Flag);
    printf("F(x)        %f              %f    %d       %f\n", *result, errest, nofun, flag);

    // нахождение и вывод интеграла от полинома Лагранжа
    quanc8(FUN2, -1.0, 1.0, abserr, relerr, result, Errest, Nofan, Flag);
    printf("Lagrange    %f             %f    %d       %f\n",  *result, errest, nofun, flag);

    // нахождение и вывод интеграла от сплайнов
    quanc8(FUN3, -1.0, 1.0, abserr, relerr, result, Errest, Nofan, Flag);
    printf("Spline      %f              %f    %d      %f\n", *result, errest, nofun, flag);

    return 0;
    
}
