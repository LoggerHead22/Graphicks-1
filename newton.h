#ifndef NEWTON_H
#define NEWTON_H

#include <QVector>
#include <functional>
using namespace std;

QVector<double> newtonXs(double a, double b, int n) {
    QVector<double> xs(n + 1);
    double step = (b - a) / n;
    for(int i = 0; i <= n; i++){
        xs[i] = a + step*i;
    }
    return xs;
}

QVector<double> newtonYs(function<double(double)> f, const QVector<double>& xs, int n) {
    QVector<double> ys(n + 1);
    for(int i = 0; i <= n; i++){
        ys[i] = f(xs[i]);
    }
    return ys;
}

QVector<double> Newton222(int n, const QVector<double>& x, const QVector<double>& f){

    QVector<double> coef = f;

    double step = x[1] - x[0];

    //cout<<endl;
    //print_vector(coef, n + 1);
    for(int i = 1; i<= n ; i++){
        for( int k = i;k <= n; k++){
            coef[n + i - k] = (coef[n + i - k] - coef[n + i - 1 - k])/(step*i);

        }
        //print_vector(coef, n + 1);
    }

    return coef;
}

double f_aprox222(double x, double a, double b, int n, const QVector<double>& coef) {
    double step = (b - a) / n;
    double res = 0;
    for(int i= n; i >= 0 ;i--){
        res = res*(x - (a + step*(i))) + coef[i];
    }
    return res;
}

#endif // NEWTON_H
