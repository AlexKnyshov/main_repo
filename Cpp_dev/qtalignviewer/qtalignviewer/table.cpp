#include "extern.h"
#include "table.h"
#include "functions.h"
#include <QTime>

int ROWS; //for some reason, first mentioned here
int COLS;

MyModel::MyModel(QObject *parent)
    :QAbstractTableModel(parent)
{
    timer = new QTimer(this);
    timer->setInterval(2000);
    connect(timer, SIGNAL(timeout()) , this, SLOT(timerHit()));
    timer->start();
}

int MyModel::rowCount(const QModelIndex & /*parent*/) const
{
   return ROWS;
}

int MyModel::columnCount(const QModelIndex & /*parent*/) const
{
    return COLS;
}

QVariant MyModel::data(const QModelIndex &index, int role) const
{
    int row = index.row();
    int col = index.column();
    switch(role){
    case Qt::DisplayRole:
        return  QString::fromStdString(vector1[index.row()][index.column()]);
        break;
    case Qt::BackgroundRole:
        std::map<std::string, int> m = {{"A", 1}, {"C", 2}, {"T", 3}, {"G", 4},{"a", 1}, {"c", 2}, {"t", 3}, {"g", 4}, {"*", 5}};
        int nt = m[vector1[index.row()][index.column()]];
        QBrush redBackground(Qt::red);
        QBrush greenBackground(Qt::green);
        QBrush blueBackground(Qt::blue);
        QBrush yellowBackground(Qt::yellow);
        QBrush blackBackground(Qt::black);
        switch(nt){
        case 1:
            return redBackground;
            break;
        case 2:
            return greenBackground;
            break;
        case 3:
            return blueBackground;
            break;
        case 4:
            return yellowBackground;
            break;
        case 5:
            return blackBackground;
            break;
        }
    }
    return QVariant();
}
bool MyModel::setData(const QModelIndex & index, const QVariant & value, int role)
{
    if (role == Qt::EditRole)
    {
        //save value from editor to member m_gridData
        QString result = value.toString();
        vector1[index.row()][index.column()] = result.toUtf8().constData();
    }
    return true;
}
Qt::ItemFlags MyModel::flags(const QModelIndex & /*index*/) const
{
    return Qt::ItemIsSelectable |  Qt::ItemIsEditable | Qt::ItemIsEnabled ;
}
void MyModel::timerHit()
{
    //we identify the top left cell
    QModelIndex topLeft = createIndex(0,0);
    QModelIndex bottomRight = createIndex(ROWS,COLS);
    //emit a signal to make the view reread identified data
    emit dataChanged(topLeft, bottomRight);
}
