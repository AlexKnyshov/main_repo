#include "extern.h"
#include "table.h"
#include "functions.h"
#include <QTime>

std::vector< std::vector<std::string> > vector1;

MyModel::MyModel(QObject *parent)
    :QAbstractTableModel(parent)
{
    std::map <std::string, std::string> test = readfastafile("winglesscat.fas");
    std::cout << "call func";
    //fastalen(test);
    //for ( auto x : test){
      //cout << x.first << string(maxsize - x.first.length(), ' ') << " " << x.second.substr(0, 30) << endl;
      //m_gridData[0][0] = QString::fromStdString(x.first);
      //m_gridData[0][1] = QString::fromStdString(x.second.substr(0, 1));
      //m_gridData[0][2] = QString::fromStdString(x.second.substr(1, 1));
    //}

    //std::vector< std::vector<std::string> > vector1;
    std::cout << "fillvec called" << test.size() << std::endl;
    fillvector (vector1, test);
    std::cout << "filling..." << vector1.size() << std::endl;
//    for (int x = 0; x < ROWS; x++)
//    {
//        for (int y = 0; y < COLS; y++)
//        {
//            //std::cout << "test" << x << " " << y << " " << vector1[x][y] << std::endl;
//            m_gridData[x][y] = QString::fromStdString(vector1[x][y]);
//        }
//    }
    //m_gridData[0][1] = QString::fromStdString(valuetogrid(0, 1));
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
        //return m_gridData[index.row()][index.column()];
        return  QString::fromStdString(vector1[index.row()][index.column()]);
        break;
    case Qt::BackgroundRole:
        if (row % 10 == 0 && col % 10 == 0)
        {
            QBrush redBackground(Qt::red);
            return redBackground;
        }
        break;
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
        //m_gridData[index.row()][index.column()] = value.toString();
        //for presentation purposes only: build and emit a joined string
//        QString result;
//        for(int row= 0; row < ROWS; row++)
//        {
//            for(int col= 0; col < COLS; col++)
//            {
//                result += m_gridData[row][col] + " ";
//            }
//        }
        //MyModel->setData( index.row(), index.column(), m_gridData[index.row()][index.column()] );
        //emit editCompleted( result );
        //emit dataChanged(index.row(), index.column());
        //emit layoutAboutToBeChanged();
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
    //emit a signal to make the view reread identified data
    std::map <std::string, std::string> test = readfastafile("winglesscat.fas");
    fillvector (vector1, test);
    emit dataChanged(topLeft, topLeft);
}
