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

    //std::map <std::string, std::string> test = readfastafile("winglesscat.fas");
    std::map <std::string, std::string> test = readfastafile("COI_aligned");
    std::cout << "call func";
    std::cout << "fillvec called" << test.size() << std::endl;
    //fillvector (vector1, test);
    int size = fastalen(test);
    COLS = size;
    std::cout << "columns " << COLS << std::endl;
    int rowsize = test.size();
    ROWS = rowsize;
    std::cout << "rows " << ROWS << std::endl;
    fillvector (vector1, test);



    //add table
    QTableView *tableView;
    tableView = new QTableView(centralWidget);
    //add model
    QAbstractTableModel *myModel = new MyModel(this);
    tableView->setModel(myModel);
    //add header
    QHeaderView* header = tableView->horizontalHeader();
    //header->setSectionResizeMode(QHeaderView::Stretch);
    for (int col=1; col<COLS; col++)
    {
       tableView->setColumnWidth(col,20);
    }
    //header->setSectionResizeMode(QHeaderView::Stretch);
    //connect model to window title
    connect(myModel, SIGNAL(editCompleted(const QString &)), this, SLOT(setWindowTitle(const QString &)));

    //add and connect button
    QPushButton *train_button = new QPushButton();
    connect(train_button, SIGNAL(released()), this, SLOT(testfunc()));
    train_button->setText(tr("Load the file"));

    //add and connect button
    QPushButton *save_button = new QPushButton();
    connect(save_button, SIGNAL(released()), this, SLOT(savefunc()));
    save_button->setText(tr("Save the file"));

    //set layouts
    layout->addWidget(train_button);
    layout->addWidget(save_button);
    layout->addWidget(tableView);

    setWindowTitle(tr("Basic Layouts"));
}

//void MainWindow::showWindowTitle(const QString & title)
//{
//setWindowTitle(title);
//}
void MainWindow::testfunc()
{
//    int size;
    //std::map <std::string, std::string> test = readfastafile("winglesscat.fas");
    std::map <std::string, std::string> test = readfastafile("COI_aligned");
    //std::vector< std::vector<std::string> > vector1;
//    size = fastalen(test);
//    int *pArray = new int[size];
    //return;
    //MyModel::MyModel->m_gridData[0][2] = QString::number(size);
    updatevector(vector1, test);
}

void MainWindow::savefunc()
{
//    int size;
    //std::map <std::string, std::string> test = readfastafile("winglesscat.fas");
    //std::map <std::string, std::string> test = readfastafile("COI_aligned");
    //std::vector< std::vector<std::string> > vector1;
//    size = fastalen(test);
//    int *pArray = new int[size];
    //return;
    //MyModel::MyModel->m_gridData[0][2] = QString::number(size);
    //updatevector(vector1, test);
    std::cout << "savefunc called" << std::endl;
    savefile(vector1, "testicula");
}
