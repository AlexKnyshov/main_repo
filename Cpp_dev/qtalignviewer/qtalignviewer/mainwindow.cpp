#include <QTableView>
#include <QHeaderView>
#include <QPushButton>
#include "mainwindow.h"
#include "table.h"
#include "functions.h"
#include "extern.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    QWidget *centralWidget = new QWidget(this);
    this->setCentralWidget( centralWidget );
    QVBoxLayout *layout = new QVBoxLayout(centralWidget);
    QTableView *tableView;
    //QVBoxLayout *layout = new QVBoxLayout();
    tableView = new QTableView(centralWidget);
    //setCentralWidget(tableView);
    QAbstractTableModel *myModel = new MyModel(this);
    tableView->setModel(myModel);
    QHeaderView* header = tableView->horizontalHeader();
    header->setSectionResizeMode(QHeaderView::Stretch);
    //readfastafile("winglesscat.fas");
    //transfer changes to the model to the window title
    connect(myModel, SIGNAL(editCompleted(const QString &)), this, SLOT(setWindowTitle(const QString &)));
    //QVBoxLayout *layout = new QVBoxLayout(tableView);
    QPushButton *train_button = new QPushButton();
    connect(train_button, SIGNAL(released()), this, SLOT(testfunc()));
    //connect(train_button, SIGNAL(released()), tableView, SLOT(timerHit()));
    train_button->setText(tr("something"));
    layout->addWidget(train_button);
    layout->addWidget(tableView);
    //setLayout(layout);
    setWindowTitle(tr("Basic Layouts"));
}

//void MainWindow::showWindowTitle(const QString & title)
//{
//setWindowTitle(title);
//}
void MainWindow::testfunc()
{
//    int size;
    std::map <std::string, std::string> test = readfastafile("winglesscat.fas");
    std::vector< std::vector<std::string> > vector1;
//    size = fastalen(test);
//    int *pArray = new int[size];
    //return;
    //MyModel::MyModel->m_gridData[0][2] = QString::number(size);
    fillvector(vector1, test);
}
