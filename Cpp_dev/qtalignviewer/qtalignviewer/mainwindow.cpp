#include <QTableView>
#include <QHeaderView>
#include <QPushButton>
#include "mainwindow.h"
#include "table.h"
#include "functions.h"
#include "extern.h"

std::vector< std::vector<std::string> > vector1;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    QWidget *centralWidget = new QWidget(this); //add mainwidget
    this->setCentralWidget( centralWidget ); //cet it central
    QVBoxLayout *layout = new QVBoxLayout(centralWidget); //add box layout

    //add table
    QTableView *tableView;
    tableView = new QTableView(centralWidget);
    //add model
    QAbstractTableModel *myModel = new MyModel(this);
    tableView->setModel(myModel);
    //add header
    QHeaderView* header = tableView->horizontalHeader();
    header->setSectionResizeMode(QHeaderView::Stretch);
    //connect model to window title
    connect(myModel, SIGNAL(editCompleted(const QString &)), this, SLOT(setWindowTitle(const QString &)));

    //add and connect button
    QPushButton *train_button = new QPushButton();
    connect(train_button, SIGNAL(released()), this, SLOT(testfunc()));
    train_button->setText(tr("something"));

    //set layouts
    layout->addWidget(train_button);
    layout->addWidget(tableView);

    setWindowTitle(tr("Basic Layouts"));

    std::map <std::string, std::string> test = readfastafile("winglesscat.fas");
    std::cout << "call func";
    std::cout << "fillvec called" << test.size() << std::endl;
    fillvector (vector1, test);

}

//void MainWindow::showWindowTitle(const QString & title)
//{
//setWindowTitle(title);
//}
void MainWindow::testfunc()
{
//    int size;
    std::map <std::string, std::string> test = readfastafile("winglesscat.fas");
    //std::vector< std::vector<std::string> > vector1;
//    size = fastalen(test);
//    int *pArray = new int[size];
    //return;
    //MyModel::MyModel->m_gridData[0][2] = QString::number(size);
    updatevector(vector1, test);
}
