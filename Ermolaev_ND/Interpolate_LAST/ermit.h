#ifndef ERMIT_H
#define ERMIT_H

#include <QVector>
#include <tuple>
#include <QDebug>
#include <functional>
#include "newton.h"
using namespace std;





double diff_Newton(double y, int k, const QVector<double>& xs , const QVector<double> coef ){
    double Ly = f_aprox222(y,xs[0] , xs[k] , k , coef);
    auto Ln = coef;

    if(abs(xs[0] - y) > 1e-16){
        Ln[0] = 0;
    }else{
    Ln[0] = (Ln[0] - Ly)/( xs[0] - y);
    }



    for (int i = 1 ; i < k; i++){
        if(abs(xs[i] - y) > 1e-16){
            Ln[i] = 0;
        }else{
          Ln[i] = (Ln[i] - Ln[i-1])/(xs[i] - y);
        }

    }

    double diff = f_aprox222(y,xs[0] , xs[k] , k, Ln);
    return diff;
}

double diff_func(function <double(double)> func ,double x){
    return (func(x + 1e-10) - func(x - 1e-10))/(2e-10);
}


std::tuple<QVector<double>, QVector<double>, QVector<double>, QVector<double>> Ermit(int n, const QVector<double>& x, const QVector<double>& f, function<double(double)> ddx){

    QVector<double> c1 = f;
    QVector<double> c2(n + 1);
    QVector<double> c3(n);
    QVector<double> c4(n);

    double step = x[1] - x[0];

    for(int i = 1; i < n; i++){
        c2[i] = ddx(x[i]);
    }

    const int k = std::min(f.size() -1 , 4);

    auto newtonCoef = Newton222(k , x, f);
    //c2[0] = diff_func(std::bind(f_aprox222, std::placeholders::_1, x[0], x[k + 1] , k, newtonCoef), x[0]);
    c2[0] = (f_aprox222(x[0] + 1e-10,x[0] , x[k] ,k,newtonCoef) - f_aprox222(x[0] - 1e-10 ,x[0] , x[k] , k , newtonCoef))/(2e-10);
    //c2[0] = diff_Newton(x[0] , k, x, newtonCoef);
	QVector<double> copy_x(k + 1); QVector<double> copy_f(k + 1);
	
	for (int i = 0; i < k + 1; i++) {
		copy_x[i] = x[x.size() - (k + 1 - i)];
		copy_f[i] = f[f.size() - (k  + 1 - i)];
	}
	
    newtonCoef = Newton222(k , copy_x, copy_f);
    c2[n] = (f_aprox222(x.back() + 1e-10 ,x[n - k] , x.back() ,k,newtonCoef) - f_aprox222(x.back() - 1e-10,x[n-k] , x.back() , k , newtonCoef))/(2e-10);

	c2[0] = ddx(x[0]);
	c2[n] = ddx(x[n]);
    for(int i = 0; i < n; i++){
        c3[i] = (3*(f[i + 1] - f[i])/(step) - 2*c2[i] - c2[i+1])/(step);
        c4[i] = (-2*(f[i + 1] - f[i])/(step) + c2[i] + c2[i+1])/(step*step);
    }


    return {c1, c2, c3, c4};
}

double f_aprox_polin(double x, double a, double b, int n, const std::tuple<QVector<double>, QVector<double>, QVector<double>, QVector<double>>& coef){
    const QVector<double>& c1 = std::get<0>(coef);
    const QVector<double>& c2 = std::get<1>(coef);
    const QVector<double>& c3 = std::get<2>(coef);
    const QVector<double>& c4 = std::get<3>(coef);

    double step = (b - a)/n;

    int i = (x - a)/step;
    if (i == n) {
        i--;
    }
    double x_i = a + i*step;
    double res = c1[i] + c2[i]*(x - x_i) + c3[i]*(x - x_i)*(x - x_i) + c4[i]*(x - x_i)*(x - x_i)*(x - x_i);

    return res;
}






#endif // ERMIT_H
