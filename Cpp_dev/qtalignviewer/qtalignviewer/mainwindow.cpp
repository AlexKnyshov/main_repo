#include <QTableView>
#include <QHeaderView>
#include <QPushButton>
#include "mainwindow.h"
#include "table.h"
#include "functions.h"
#include "extern.h"
#include <QHBoxLayout>
#include <QRadioButton>

std::vector< std::vector<std::string> > vector1;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    QWidget *centralWidget = new QWidget(this); //add mainwidget
    this->setCentralWidget( centralWidget ); //cet it central
    QVBoxLayout *layout = new QVBoxLayout(centralWidget); //add box layout

    //std::map <std::string, std::string> test = readfastafile("winglesscat.fas");
    std::map <std::string, std::string> test = readfastafile("winconcat.fas.fas");
    std::cout << "call func";
    std::cout << "fillvec called" << test.size() << std::endl;

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

    QHBoxLayout *horiz_menu = new QHBoxLayout;

    //add and connect button
    QPushButton *train_button = new QPushButton();
    connect(train_button, SIGNAL(released()), this, SLOT(testfunc()));
    train_button->setText(tr("Load the file"));

    //add and connect button
    QPushButton *save_button = new QPushButton();
    connect(save_button, SIGNAL(released()), this, SLOT(savefunc()));
    save_button->setText(tr("Save the file"));

    //add and connect button
    QPushButton *translate_button = new QPushButton();
    connect(translate_button, SIGNAL(released()), this, SLOT(translate()));
    translate_button->setText(tr("Translate"));

    //add and connect button
    QPushButton *translate_button = new QPushButton();
    connect(translate_button, SIGNAL(released()), this, SLOT(translate()));
    translate_button->setText(tr("Translate"));

    //set layouts
    //layout->addWidget(train_button);
    //layout->addWidget(save_button);
    layout->addLayout(horiz_menu);
    horiz_menu->addWidget(train_button);
    horiz_menu->addWidget(save_button);
    horiz_menu->addWidget(translate_button);

    layout->addWidget(tableView);

    setWindowTitle(tr("Basic Layouts"));
}

void MainWindow::testfunc()
{

    //std::map <std::string, std::string> test = readfastafile("winglesscat.fas");
    std::map <std::string, std::string> test = readfastafile("winconcat.fas.fas");

    updatevector(vector1, test);
}

void MainWindow::savefunc()
{

    std::cout << "savefunc called" << std::endl;
    savefile(vector1, "testicula");
}

void MainWindow::translate()
{

    std::cout << "translate() called" << std::endl;
    translatevector(vector1);
}
