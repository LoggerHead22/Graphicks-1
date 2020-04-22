#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QDebug>
#include <QPainter>
#include <QWheelEvent>
#include <cmath>
#include "newton.h"
#include "ermit.h"
#include <iostream>


QPair<double, double> functionMinMax(function<double(double)> f, QPair<double, double> xRange) {
    double start = xRange.first, end = xRange.second;
    double min = f(start), max = f(start);

    for (double x = start; x <= end; x += (xRange.second - xRange.first) / 1000) {
        if (f(x) < min) {
            min = f(x);
        }
        if (f(x) > max) {
            max = f(x);
        }
    }

    return {min, max};
}

QVector<QLineF> functionLines(function<double(double)> f, QPair<double, double> xRange) {
    QVector<QLineF> ret;
    double start = xRange.first, end = xRange.second;
    QPointF lastPoint = {start, f(start)};

    for (double x = start; x <= end; x += (xRange.second - xRange.first) / 1000) {
        ret.push_back(QLineF(lastPoint, {x, f(x)}));
        lastPoint = {x, f(x)};
    }
    ret.push_back(QLineF(lastPoint, {end, f(end)}));

    return ret;
}



void print_vector(double *coef, int n){
    for(int i = 0; i < n; i++){
        std::cout<<coef[i]<<" ";
    }
    std::cout<<endl;
}

double ddx_f0(double)   { return 0; }
double ddx_f1(double)   { return 1; }
double ddx_f2(double x) { return 2 * x; }
double ddx_f3(double x) { return 3 * x*x; }
double ddx_f4(double x) { return 4 * x*x*x; }
double ddx_f5(double x) { return exp(x); }
double ddx_f6(double x) { return -(50 * x) / pow(25 * x*x + 1, 2); }

QVector<double(*)(double)> ddxs = {ddx_f0, ddx_f1, ddx_f2, ddx_f3, ddx_f4, ddx_f5, ddx_f6};



MainWindow::MainWindow(double aa, double bb, int nn, int kk, QWidget *parent): QMainWindow(parent), ui(new Ui::MainWindow)
{
    a=aa;
    b=bb;
    n=nn;
    currentFunctionIndex = kk;


    ui->setupUi(this);
    menuBar()->addAction("Exit", this, SLOT(close()))->setShortcut(QString("Ctrl+X"));
//    QMenuBar *tool_bar = ui->menubar;
//    QAction *action;

//    action = tool_bar->addAction ("&Exit", this, SLOT (close ()));
//    action->setShortcut (QString ("Ctrl+X"));


    ui->aLineEdit->setText(QString::number(a));
    ui->bLineEdit->setText(QString::number(b));
    ui->nLineEdit->setText(QString::number(n));

    connect(ui->aLineEdit, &QLineEdit::returnPressed, this, [&]() { a = ui->aLineEdit->text().toDouble(); view = {a, b}; compute_Newton(); });
    connect(ui->bLineEdit, &QLineEdit::returnPressed, this, [&]() { b = ui->bLineEdit->text().toDouble(); view = {a, b}; compute_Newton(); });
    connect(ui->nLineEdit, &QLineEdit::returnPressed, this, [&]() { n = ui->nLineEdit->text().toInt(); compute_Newton(); });

    connect(ui->mlt2PushButton, &QPushButton::released, this, [&]() { n *= 2; ui->nLineEdit->setText(QString::number(n));compute_Newton(); });
    connect(ui->div2PushButton, &QPushButton::released, this, [&]() { n = (n>1?n/2:n); ui->nLineEdit->setText(QString::number(n));compute_Newton(); });

    functions.push_back({[](double)   { return 1; }, "f(x) = 1"});
    functions.push_back({[](double x)   { return x; }, "f(x) = x"});
    functions.push_back({[](double x) { return x * x; }, "f(x) = x ^ 2"});
    functions.push_back({[](double x) { return x * x * x; }, "f(x) = x ^ 3"});
    functions.push_back({[](double x) { return x * x * x * x; }, "f(x) = x ^ 4"});
    functions.push_back({[](double x) { return std::exp(x); }, "f(x) = e ^ x"});
    functions.push_back({[](double x) { return 1.0 / (25 * x * x + 1); }, "f(x) = 1 / (25 * x ^ 2 + 1)"});

    ui->functionNameLabel->setText(functions[currentFunctionIndex].second);

    connect(ui->changeFunctionPushButton, &QPushButton::released, this, [&]() {
        currentFunctionIndex = (currentFunctionIndex + 1) % functions.size();
        ui->functionNameLabel->setText(functions[currentFunctionIndex].second);
        view = {a, b} ; compute_Newton();
        noize = 0;
    });

    connect(ui->resetZoomPushButton, &QPushButton::released, this, [&]() { view = {a, b}; });

    view = {a, b};
    x = newtonXs(a,b,n);
    f = newtonYs(currentFunction(),x,n);
    if(n < 50){
        coef = Newton222(n,x,f);
    }
    std::tie(c1,c2,c3,c4) = Ermit(n,x,f,ddxs[currentFunctionIndex]);

    computeResidual();


}

MainWindow::~MainWindow() {
    delete ui;

}

void MainWindow::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event);

    Q_UNUSED(event);

    QVector<QLineF> fnLines, apprxLines, apprx222Lines, fnMinusApprxLines, fnMinusApprx222Lines;

    std::tuple<QVector<double>,QVector<double>,QVector<double>,QVector<double>> cs  (c1,c2,c3,c4);

    fnLines = functionLines(currentFunction(), view);
    if (kit == 0 || kit == 2) apprxLines = functionLines(std::bind(f_aprox_polin,std::placeholders::_1,a,b,n,cs), view);

    if (kit == 1 || kit == 2) if (n < 50) apprx222Lines = functionLines(std::bind(f_aprox222, std::placeholders::_1, a, b, n,coef), view);
    if (kit == 3) fnMinusApprxLines = functionLines([&](double x) { return std::abs(currentFunction()(x) - f_aprox_polin(x,a,b,n,cs)); }, view);
    if (kit == 3) if( n < 50) fnMinusApprx222Lines = functionLines([&](double x) { return std::abs(currentFunction()(x) - f_aprox222(x, a, b, n, coef)); }, view);

    QPainter painter(this);
    QRectF plotRect = ui->plotWidget->geometry(); //окно с графиком
    fnMinMax = functionMinMax(currentFunction(), view);

    QPair<double, double> apprxMinMax, apprx222MinMax, fnMinusApprxMinMax, fnMinusApprx222MinMax;
    if (kit == 0 || kit == 2) apprxMinMax = functionMinMax(std::bind(f_aprox_polin,std::placeholders::_1,a,b,n,cs), view);
    if (kit == 1 || kit == 2) if (n < 50) apprx222MinMax = functionMinMax(std::bind(f_aprox222, std::placeholders::_1, a, b, n, coef), view);
    if (kit == 3) fnMinusApprxMinMax = functionMinMax([&](double x) { return std::abs(currentFunction()(x) - f_aprox_polin(x,a,b,n,cs)); }, view);
    if (kit == 3) if (n < 50) fnMinusApprx222MinMax = functionMinMax([&](double x) { return std::abs(currentFunction()(x) - f_aprox222(x, a, b, n, coef)); }, view);

    if (kit == 0) {
        minMax = {qMin(fnMinMax.first, apprxMinMax.first),
                  qMax(fnMinMax.second, apprxMinMax.second)};
    }
    if (kit == 1) {
        if (n < 50) minMax = {qMin(fnMinMax.first, apprx222MinMax.first),
                  qMax(fnMinMax.second, apprx222MinMax.second)};
        else minMax = fnMinMax;
    }
    if (kit == 2) {
        minMax = {qMin(fnMinMax.first, apprxMinMax.first),
                  qMax(fnMinMax.second, apprxMinMax.second)};
        if (n < 50) minMax = {qMin(minMax.first, apprx222MinMax.first),
                              qMax(minMax.second, apprx222MinMax.second)};
    }
    if (kit == 3) {
        if( n < 50){
            minMax = {qMin(fnMinusApprxMinMax.first, fnMinusApprx222MinMax.first),
                      qMax(fnMinusApprxMinMax.second, fnMinusApprx222MinMax.second)};
        }else {
            minMax = fnMinusApprxMinMax;
        }
    }

    if (std::abs(minMax.first - minMax.second) < 1e-15) {
        minMax.first = minMax.first - 1e-15;
        minMax.second = minMax.second + 1e-15;
    }

//    if (kit == 0) painter.drawText(1200, 40, "Набор = ф-ия + ermit");
//    if (kit == 1) painter.drawText(1200, 40, "Набор = ф-ия + newton");
//    if (kit == 2) painter.drawText(1200, 40, "Набор = погрешность newton");

    QRectF viewRect = {QPointF(view.first, minMax.first), QPointF(view.second, minMax.second)}; //границы кординат

    QPen axisPen(Qt::red, 0);
    QPen plotPen(Qt::darkMagenta, 0);
    QPen AprPen(Qt::green, 0);
    QPen numPen(Qt::black, 0);
    QPen gridPen(Qt::lightGray, 0);

    painter.fillRect(plotRect, Qt::white);
    painter.setPen(QPen(Qt::black, 2));
    painter.drawRect(plotRect);

    int numCount = 10;
    for (double i = plotRect.left(), x = view.first;
        i <= plotRect.right() + 1;
        i += plotRect.width() / numCount, x += (view.second - view.first) / numCount)
    {
        painter.setPen(numPen);
        painter.drawText(i, plotRect.bottom() + painter.fontMetrics().height(), QString::number(x, 'e', 4));
        painter.setPen(gridPen);
        painter.drawLine(i, plotRect.top(), i, plotRect.bottom());
    }
    for (double i = plotRect.bottom(), y = minMax.first;
        i >= plotRect.top() - 1;
        i -= plotRect.height() / numCount, y += (minMax.second - minMax.first) / numCount)
    {
        painter.setPen(numPen);
        painter.drawText(plotRect.left(), i, QString::number(y, 'e', 4));
        painter.setPen(gridPen);
        painter.drawLine(plotRect.left(), i, plotRect.right(), i);
    }

    painter.translate(plotRect.center());
    painter.scale(plotRect.width() / (view.second - view.first), plotRect.height() / -(minMax.second - minMax.first));
    painter.translate(-viewRect.center());
    if (view.first * view.second < 0) {
        painter.setPen(axisPen);
        painter.drawLine(QLineF(0, minMax.first, 0, minMax.second));
    }
    if (minMax.first * minMax.second < 0) {
        painter.setPen(axisPen);
        painter.drawLine(QLineF(view.first, 0, view.second, 0));
    }

    if (kit == 0) {
        painter.setPen(plotPen);
        painter.drawLines(fnLines);
        painter.setPen(AprPen);
        painter.drawLines(apprxLines);
    }
    if (kit == 1) {
        painter.setPen(plotPen);
        painter.drawLines(fnLines);
        if (n < 50) {
            painter.setPen(QPen(Qt::blue, 0));
            painter.drawLines(apprx222Lines);
        }
    }
    if (kit == 2) {
        painter.setPen(plotPen);
        painter.drawLines(fnLines);
        painter.setPen(AprPen);
        painter.drawLines(apprxLines);
        if (n < 50) {
            painter.setPen(QPen(Qt::blue, 0));
            painter.drawLines(apprx222Lines);
        }
    }
    if (kit == 3) {
        painter.setPen(AprPen);
        painter.drawLines(fnMinusApprxLines);
        if( n < 50){
        painter.setPen(QPen(Qt::blue, 0));
        painter.drawLines(fnMinusApprx222Lines);
        }
    }

    if (isZooming) {
        double viewWidth = view.second - view.first;
        double zoomEnd = view.first + double(ui->plotWidget->mapFromGlobal(QCursor::pos()).x()) / ui->plotWidget->geometry().width() * viewWidth;
        zoomEnd = qMin(zoomEnd, view.second);
        zoomEnd = qMax(zoomEnd, view.first);
        painter.fillRect(QRectF(QPointF(zoomStart, minMax.first), QPointF(zoomEnd, minMax.second)), QBrush(QColor(255, 128, 8, 128)));
    }
    //painter.restore ();
    ui->residualLabel->setText(prevText +
                               "\nfn max: " + QString::number(fnMinMax.second) +
                               "\nk = " + QString::number(k) +
                               "\nmax{|Fmin|, |Fmax|} = " + QString::number(std::max(std::abs(fnMinMax.first), std::abs(fnMinMax.second))) +
                               "\nНабор: " + QVector<QString>{"ф-ия + Эрмит", "ф-ия + Ньютон", "ф-ия + Эрмит + Ньютон", "ошибки Эрмита и Ньютона"}[kit] +
                               "\nМасшбтаб: " + QString::number(std::log2((b - a)/(view.second - view.first))) +
                               "\np = " + QString::number(noize));
    update();
}

void MainWindow::keyReleaseEvent(QKeyEvent *event) {
    if (event->key() == Qt::Key_0) {
        ui->changeFunctionPushButton->animateClick();
    }
    if (event->key() == Qt::Key_1) {
        kit = (kit + 1) % 4;
        update();
    }
    if (event->key() == Qt::Key_2) {
        double center = (view.first + view.second) / 2.0;
        double r = (view.second - view.first) / 2.0;       
        view = {center - r / 2, center + r / 2};
        update();
    }
    if (event->key() == Qt::Key_3) {
        double center = (view.first + view.second) / 2.0;
        double r = (view.second - view.first) / 2.0;
        if (center - r * 2 >= a) {
            view = {center - r * 2, center + r * 2};
        }
        update();
    }
    if (event->key() == Qt::Key_4) {
        ui->mlt2PushButton->animateClick();
    }
    if (event->key() == Qt::Key_5) {
        ui->div2PushButton->animateClick();
    }
    if (event->key() == Qt::Key_6 && kit != 3) {
        noize += fnMinMax.second / 10;
        compute_Newton();
    }
    if (event->key() == Qt::Key_7 && kit != 3) {
        noize -= fnMinMax.second / 10;
        compute_Newton();
    }
}

void MainWindow::compute_Newton(){

    x = newtonXs(a,b,n);
    f = newtonYs(currentFunction(),x,n);

    if(n < 50){
        coef = Newton222(n,x,f);
    }
    std::tie(c1,c2,c3,c4) = Ermit(n,x,f,ddxs[currentFunctionIndex]);

    computeResidual();
}

void MainWindow::mousePressEvent(QMouseEvent *event) {
    double viewWidth = view.second - view.first;
    zoomStart = view.first + double(ui->plotWidget->mapFromGlobal(event->globalPos()).x()) / ui->plotWidget->geometry().width() * viewWidth;
    zoomStart = qMin(zoomStart, view.second);
    zoomStart = qMax(zoomStart, view.first);
    isZooming = true;
}

void MainWindow::mouseReleaseEvent(QMouseEvent *event) {
    Q_UNUSED(event);
    double viewWidth = view.second - view.first;
    double zoomEnd = view.first + double(ui->plotWidget->mapFromGlobal(QCursor::pos()).x()) / ui->plotWidget->geometry().width() * viewWidth;
    zoomEnd = qMin(zoomEnd, view.second);
    zoomEnd = qMax(zoomEnd, view.first);

    if (zoomStart > zoomEnd) {
        std::swap(zoomStart, zoomEnd);
    }

    if (std::abs(zoomStart - zoomEnd) > 1e-4) {
        view.first = zoomStart;
        view.second = zoomEnd;
    }

    isZooming = false;
}

function<double(double)> MainWindow::currentFunction() {
    auto ret = [&](double x) -> double {
        if (std::abs(x - (view.first + view.second) / 2) < 1e-6) {
            return functions[currentFunctionIndex].first(x) + noize;
        }
        return functions[currentFunctionIndex].first(x);
    };
    return ret;
}

void MainWindow::computeResidual() {
    QPair<double, double> fnMinusApprxMinMax, fnMinusApprx222MinMax;
    auto xs = x;
    auto ys = f;
    std::tuple<QVector<double>,QVector<double>,QVector<double>,QVector<double>> cs  (c1,c2,c3,c4);
    fnMinusApprxMinMax = functionMinMax([&](double x) { return std::abs(currentFunction()(x) - f_aprox_polin(x,a,b,n,cs)); }, view);
    if(n < 50 ) fnMinusApprx222MinMax = functionMinMax([&](double x) { return std::abs(currentFunction()(x) - f_aprox222(x, a, b, n, coef)); }, view);
    std::cout << "ermit error:" << qMax(std::abs(fnMinusApprxMinMax.first), std::abs(fnMinusApprxMinMax.second));
    if (n < 50) std::cout << "newton error:" << qMax(std::abs(fnMinusApprx222MinMax.first), std::abs(fnMinusApprx222MinMax.second));
    QString newtonResidualAsStr = "n >= 50";
    if (n < 50) {
        double newtonResidual = qMax(std::abs(fnMinusApprx222MinMax.first), std::abs(fnMinusApprx222MinMax.second));
        newtonResidualAsStr = QString::number(newtonResidual);
    }
    ui->residualLabel->setText(QString("ermit error: ") + QString::number(qMax(std::abs(fnMinusApprxMinMax.first), std::abs(fnMinusApprxMinMax.second))) + "\nnewton error: " + newtonResidualAsStr);
    prevText = ui->residualLabel->text();
    std::cout <<  "max{|Fmin|, |Fmax|} = " << std::max(std::abs(fnMinMax.first), std::abs(fnMinMax.second));
}



