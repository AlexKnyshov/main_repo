#include <QTableView>
#include <QHeaderView>
#include <QPushButton>
#include "mainwindow.h"
#include "table.h"
#include "functions.h"
#include "extern.h"
#include <QHBoxLayout>
#include <QRadioButton>
#include <QGroupBox>
#include <QShortcut>
std::vector< std::vector<std::string> > vector1;

int opt = 1;
    QTableView *tableView;
    QAbstractTableModel *myModel;

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


    std::cout << "opt value - " << opt << std::endl;


    //add table

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
    //header->setSectionResizeMode(QHeaderView::Fixed);
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

    //group box
    QGroupBox *groupBox = new QGroupBox(tr("Translate options"));
    QRadioButton *radio1 = new QRadioButton(tr("Option 1"));
    QObject::connect(radio1,SIGNAL(clicked(bool)),this,SLOT(translate_option1(bool)));
    //radio1->setAutoExclusive(false);
    QRadioButton *radio2 = new QRadioButton(tr("Option 2"));
    QObject::connect(radio2,SIGNAL(clicked(bool)),this,SLOT(translate_option2(bool)));
    //radio2->setAutoExclusive(false);
    QRadioButton *radio3 = new QRadioButton(tr("Option 3"));
    //radio3->setAutoExclusive(false);
    QObject::connect(radio3,SIGNAL(clicked(bool)),this,SLOT(translate_option3(bool)));
    radio1->setChecked(true);
    radio2->setChecked(false);
    radio3->setChecked(false);
    QVBoxLayout *vbox = new QVBoxLayout;
    vbox->addWidget(radio1);
    vbox->addWidget(radio2);
    vbox->addWidget(radio3);
    //vbox->addStretch(1);
    groupBox->setLayout(vbox);
    //add and connect button
    //QRadioButton *translate_option1 = new QRadioButton("option1");
    //connect(translate_option1, SIGNAL(toggled()), this, SLOT(translate()));
    //translate_button->setText(tr("Translate"));

    //set layouts
    //layout->addWidget(train_button);
    //layout->addWidget(save_button);
    layout->addLayout(horiz_menu);
    horiz_menu->addWidget(train_button);
    horiz_menu->addWidget(save_button);
    horiz_menu->addWidget(translate_button);
    //horiz_menu->addWidget(translate_option1);
    horiz_menu->addWidget(groupBox);

    layout->addWidget(tableView);

    setWindowTitle(tr("Basic Layouts"));



    QShortcut* shortcut = new QShortcut(QKeySequence(Qt::Key_A), tableView);
    connect(shortcut, SIGNAL(activated()), this, SLOT(deleteRow()));








}

void MainWindow::testfunc()
{

    //std::map <std::string, std::string> test = readfastafile("winglesscat.fas");
    std::map <std::string, std::string> test = readfastafile("winconcat.fas.fas");

    updatevector(vector1, test);
    for (int col=1; col<COLS; col++)
    {
        tableView->setColumnWidth(col,20);
        std::cout << "test funct goes " << col << " " << tableView->columnWidth(col) << std::endl;
    }
//    QModelIndex topLeft = tableView->model()->QAbstractItemModel::createIndex(0,0);
//    QModelIndex bottomRight = createIndex(ROWS,COLS);
//    //emit a signal to make the view reread identified data
//    tableView->model()->dataChanged(topLeft, bottomRight);
    tableView->reset();

    //myModel->select();
}

void MainWindow::savefunc()
{

    std::cout << "savefunc called" << std::endl;
    savefile(vector1, "testicula");
}

void MainWindow::translate()
{

    std::cout << "translate() called" << std::endl;
    switch(opt)
    {
        case 1:
            std::cout << "opt 1 will be executed";
            break;
        case 2:
            std::cout << "opt 2 will be executed";
            break;
        case 3:
            std::cout << "opt 3 will be executed";
            break;
    }
    int oldcols = translatevector(vector1, opt);
    tableView->reset();
    tableView->model()->headerDataChanged(Qt::Horizontal, 0 , oldcols);
//    QModelIndex topLeft = createIndex(0,0);
//    QModelIndex bottomRight = createIndex(ROWS,COLS);
//    //emit a signal to make the view reread identified data
//    tableView->model()->dataChanged(topLeft, bottomRight);
}
void MainWindow::translate_option1(bool)
{
    opt = 1;
    std::cout << "translate(1) called " << opt << std::endl;

    //translatevector(vector1);
}
void MainWindow::translate_option2(bool)
{

    opt = 2;
    std::cout << "translate(2) called " << opt << std::endl;
    //translatevector(vector1);
}
void MainWindow::translate_option3(bool)
{
    opt = 3;
    std::cout << "translate(3) called " << opt << std::endl;
    //translatevector(vector1);
}
void MainWindow::deleteRow()
{
    std::cout << "deleterow called" << std::endl;
    //QModelIndex idx = tableView->currentIndex();
    //if (idx.isValid())
    //   tableView->model()->removeRow(idx.row(), idx.parent());
        std::cout << tableView->currentIndex().row() << std::endl;
        vector1[tableView->currentIndex().row()][tableView->currentIndex().column()]="A";
        tableView->model()->dataChanged(tableView->currentIndex(), tableView->currentIndex());
}
