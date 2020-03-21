#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <functional>
using namespace std;

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(double aa, double bb, int nn, QWidget *parent = nullptr);
    ~MainWindow();

    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void compute_Newton();
    function<double(double)> currentFunction();
    void computeResidual();


    double a = -3, b = 3;
    int n = 5 , k = 4;
    double stepp = (b - a) / n;
    double newton_a = 0 , newton_b = 0;
    double* coef = new double[k+1];
    void Newton(double a, double b , int n);
    double f_aprox(double x);
    void Ermit(double a, double b , int n);
    double f_aprox_polin(double x);

    vector<double> c1;
    vector<double> c2;
    vector<double> c3;
    vector<double> c4;


    bool drawFunction = true, drawApproximated = true;

    QVector<QPair<function<double(double)>, QString>> functions;
    int currentFunctionIndex = 0;

    QPair<double, double> view;

    double zoomStart;
    bool isZooming = false;

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H

void print_vector(double *coef, int n);
double diff_func(function <double(double)> func,double x);
