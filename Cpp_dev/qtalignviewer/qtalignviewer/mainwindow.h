#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include "mainwindow.h"

QT_BEGIN_NAMESPACE // QT_BEGIN_NAMESPACE / QT_END_NAMESPACE are not needed in Qt user code
class QTableView; //forward declaration
QT_END_NAMESPACE


class MainWindow : public QMainWindow
{
    Q_OBJECT
private:
    QTableView *tableView;
public:
    MainWindow(QWidget *parent = 0);
public slots:
    void showWindowTitle(const QString & title);
};
#endif // MAINWINDOW_H
