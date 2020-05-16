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
    MainWindow(double aa, double bb, int nn, int kk, QWidget *parent = nullptr);
    ~MainWindow();

    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void keyReleaseEvent(QKeyEvent* event) override;
    void compute_Newton();
    function<double(double)> currentFunction();
    void computeResidual();


    double a = -3, b = 3;
    int n = 5 , k = 4;
    int kit = 0;
    double stepp = (b - a) / n;
    double newton_a = 0 , newton_b = 0;
    QPair<double, double> minMax;
    double noize = 0;

    QPair<double, double> fnMinMax;
    QString prevText;


    QVector<double> coef;
    QVector <double> x;
    QVector <double> f;
    QVector<double> c1;
    QVector<double> c2;
    QVector<double> c3;
    QVector<double> c4;

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

