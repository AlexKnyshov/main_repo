#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <QVBoxLayout>

QT_BEGIN_NAMESPACE // QT_BEGIN_NAMESPACE / QT_END_NAMESPACE are not needed in Qt user code
class QTableView; //forward declaration
QT_END_NAMESPACE

// native code

//class MainWindow : public QMainWindow
//{
//    Q_OBJECT
//private:
//    QTableView *tableView;
//public:
//    MainWindow(QWidget *parent = 0);
//public slots:
//    void showWindowTitle(const QString & title);
//};


class MainWindow : public QMainWindow
{
    Q_OBJECT
    QWidget *centralWidget;
    QVBoxLayout *layout;
public:
    MainWindow(QWidget *parent = 0);
    //~MainWindow();
public slots:
    void testfunc(); //button func
    void savefunc(); //save button func
private:

};

#endif // MAINWINDOW_H
