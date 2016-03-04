#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <QVBoxLayout>


#include <QAbstractTableModel>

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
    void translate(); //save button func
    void translate_option1(bool); //set transl option
    void translate_option2(bool); //set transl option
    void translate_option3(bool); //set transl option
    void deleteRow(); //set transl option
private:

};

#endif // MAINWINDOW_H
