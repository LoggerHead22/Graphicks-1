#include "mainwindow.h"

#include <QApplication>
#include <QMenuBar>
#include <QAction>
#include <QMessageBox>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    if(argc!=4){
       QMessageBox::warning (0, "Wrong input arguments!", "usage: a.out a b n");
       return -1;
    }
    double a, b;
    int n;
    if (   sscanf (argv[1], "%lf", &a) != 1
        || sscanf (argv[2], "%lf", &b) != 1
        || b - a < 1.e-2
        || (sscanf (argv[3], "%d", &n) != 1)
        || n <= 0){
      return -2;
    }


    MainWindow w(a,b,n);
    //QMenuBar tool_bar(&w);


    //tool_bar.setMaximumHeight (30);

    //w.setMenuBar (tool_bar);
    w.setWindowTitle ("Newton Approximate");
    w.show();
    return app.exec();
}
