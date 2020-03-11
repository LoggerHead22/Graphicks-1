#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QDebug>
#include <QPainter>
#include <QWheelEvent>
#include <cmath>



QPair<double, double> functionMinMax(function<double(double)> f, QPair<double, double> xRange) {
    double start = xRange.first, end = xRange.second;
    double min = f(start), max = f(start);

    for (double x = start; x <= end; x += (xRange.second - xRange.first) / 5000) {
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

    for (double x = start; x <= end; x += (xRange.second - xRange.first) / 5000) {
        ret.push_back(QLineF(lastPoint, {x, f(x)}));
        lastPoint = {x, f(x)};
    }
    ret.push_back(QLineF(lastPoint, {end, f(end)}));

    return ret;
}

double* MainWindow::Newton(double a, double b , int n){

    double *coef;
    coef = new double [n + 1];

    double step = (b - a) / n;

    for(int i = 0; i <= n; i++){
        coef[i] = currentFunction()(a + step*i);
        //cout<<coef[i]<<" ";
    }

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

void print_vector(double *coef, int n){
    for(int i = 0; i < n; i++){
        qDebug()<<coef[i]<<" ";
    }
    qDebug()<<endl;
}

double MainWindow::f_aprox(double x) {
    double step = (b - a) / n;
    double res = 0;
    for(int i= n; i >= 0 ;i--){
        res = res*(x - (a + step*(i))) + coef[i];
    }


    return res;
}



MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->aLineEdit->setText(QString::number(a));
    ui->bLineEdit->setText(QString::number(b));
    ui->nLineEdit->setText(QString::number(n));

    connect(ui->aLineEdit, &QLineEdit::returnPressed, this, [&]() { a = ui->aLineEdit->text().toDouble(); view = {a, b}; compute_Newton(); });
    connect(ui->bLineEdit, &QLineEdit::returnPressed, this, [&]() { b = ui->bLineEdit->text().toDouble(); view = {a, b}; compute_Newton(); });
    connect(ui->nLineEdit, &QLineEdit::returnPressed, this, [&]() { n = ui->nLineEdit->text().toInt(); compute_Newton(); });

    connect(ui->mlt2PushButton, &QPushButton::released, this, [&] { n *= 2; ui->nLineEdit->setText(QString::number(n));compute_Newton(); });
    connect(ui->div2PushButton, &QPushButton::released, this, [&] { n = ceil(double(n)/2); ui->nLineEdit->setText(QString::number(n));compute_Newton(); });

    ui->functionCheckBox->setChecked(drawFunction);
    ui->approximatedCheckBox->setChecked(drawApproximated);

    connect(ui->functionCheckBox, &QCheckBox::stateChanged, this, [&](int state) { drawFunction = state == Qt::Checked; });
    connect(ui->approximatedCheckBox, &QCheckBox::stateChanged, this, [&](int state) { drawApproximated = state == Qt::Checked; });

    functions.push_back({[](double x) { return x*(x - 1)*(x+1); }, "f(x) = x*(x-1)*(x+1)"});
    functions.push_back({[](double x) { return cos(x)*x*x; }, "f(x) = cos(x) * x * x"});
    functions.push_back({[](double x) { return sin(x); }, "f(x) = sin(x)"});

    ui->functionNameLabel->setText(functions[0].second);

    connect(ui->changeFunctionPushButton, &QPushButton::released, this, [&]() {
        currentFunctionIndex = (currentFunctionIndex + 1) % functions.size();
        ui->functionNameLabel->setText(functions[currentFunctionIndex].second);
        view = {a, b} ; compute_Newton();
    });

    connect(ui->resetZoomPushButton, &QPushButton::released, this, [&]() { view = {a, b}; });

    view = {a, b};
    coef = Newton(a, b, n);
    computeResidual();
}

MainWindow::~MainWindow() {
    delete ui;
    delete[] coef;
}

void MainWindow::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event);
    QPainter painter(this);
    QRectF plotRect = ui->plotWidget->geometry();
    QPair<double, double> fnMinMax = functionMinMax(currentFunction(), view);
    QPair<double, double> apprxMinMax = functionMinMax(std::bind(&MainWindow::f_aprox, this, std::placeholders::_1), view);
    QPair<double, double> minMax = {qMin(fnMinMax.first, apprxMinMax.first),
                                    qMax(fnMinMax.second, apprxMinMax.second)};
    QRectF viewRect = {QPointF(view.first, minMax.first), QPointF(view.second, minMax.second)};

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
        painter.drawText(i, plotRect.bottom() + painter.fontMetrics().height(), QString::number(x, 'f', 4));
        painter.setPen(gridPen);
        painter.drawLine(i, plotRect.top(), i, plotRect.bottom());
    }
    for (double i = plotRect.bottom(), y = minMax.first;
        i >= plotRect.top() - 1;
        i -= plotRect.height() / numCount, y += (minMax.second - minMax.first) / numCount)
    {
        painter.setPen(numPen);
        painter.drawText(plotRect.left(), i, QString::number(y, 'f', 4));
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
    if (drawFunction) {
        painter.setPen(plotPen);
        painter.drawLines(functionLines(currentFunction(), view));
    }
    if(drawApproximated){
        painter.setPen(AprPen);
        painter.drawLines(functionLines(std::bind(&MainWindow::f_aprox, this, std::placeholders::_1), view));

    }

    if (isZooming) {
        double viewWidth = view.second - view.first;
        double zoomEnd = view.first + double(ui->plotWidget->mapFromGlobal(QCursor::pos()).x()) / ui->plotWidget->geometry().width() * viewWidth;
        zoomEnd = qMin(zoomEnd, view.second);
        zoomEnd = qMax(zoomEnd, view.first);
        painter.fillRect(QRectF(QPointF(zoomStart, minMax.first), QPointF(zoomEnd, minMax.second)), QBrush(QColor(255, 128, 8, 128)));
    }

    update();
}

void MainWindow::compute_Newton(){
    delete[] coef;
    coef = Newton(a, b , n);
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
    return functions[currentFunctionIndex].first;
}


void MainWindow::computeResidual() {
    double residual = 0;
    auto lines1 = functionLines(currentFunction(), {a, b});
    auto lines2 = functionLines(bind(&MainWindow::f_aprox, this, placeholders::_1), {a, b});

    int XXX = 0;
    for (int i = 0; i < lines1.size(); i++) {
        double diff = qMax(std::abs(lines1[i].p1().y() - lines2[i].p1().y()),
                           std::abs(lines1[i].p2().y() - lines2[i].p2().y()));
        if (diff > residual) {
            XXX = i;
            residual = diff;
        }
    }
    //qDebug() << XXX;

    ui->residualLabel->setText(QString("Residual: ") + QString::number(residual));
}


