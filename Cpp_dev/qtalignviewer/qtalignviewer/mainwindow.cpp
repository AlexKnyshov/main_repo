#include <QTableView>
#include <QHeaderView>
#include "mainwindow.h"
#include "table.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    tableView = new QTableView(this);
    setCentralWidget(tableView);
    QAbstractTableModel *myModel = new MyModel(this);
    tableView->setModel(myModel);
    QHeaderView* header = tableView->horizontalHeader();
    header->setSectionResizeMode(QHeaderView::Stretch);

    //transfer changes to the model to the window title
    connect(myModel, SIGNAL(editCompleted(const QString &)), this, SLOT(setWindowTitle(const QString &)));
}

void MainWindow::showWindowTitle(const QString & title)
{
setWindowTitle(title);
}
