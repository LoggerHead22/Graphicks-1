#include "mainwindow.h"

#include <QApplication>
#include <QMenuBar>
#include <QAction>
#include <QMessageBox>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    if(argc!=5){
       QMessageBox::warning (0, "Wrong input arguments!","Usage a.exe a b n k");
       return -1;
    }
    double a, b;
    int n , k;
    if (   sscanf (argv[1], "%lf", &a) != 1
        || sscanf (argv[2], "%lf", &b) != 1
        || b - a < 1.e-2
        || (sscanf (argv[3], "%d", &n) != 1)
        || n <= 0
        || (sscanf (argv[4], "%d", &k) != 1)
          || (k <0) || (k > 6)){

       QMessageBox::warning (0, "Error!","Usage a.exe a b n k");
      return -2;
    }


    MainWindow w(a ,b , n, k);
    w.setWindowTitle ("Interpolate 1");
    w.show();
    return app.exec();
}
