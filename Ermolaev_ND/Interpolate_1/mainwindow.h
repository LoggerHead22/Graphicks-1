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
    int n = 5;

    double* coef = nullptr;
    double* Newton(double a, double b , int n);
    double f_aprox(double x);

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
