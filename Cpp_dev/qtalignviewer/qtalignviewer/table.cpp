#include "extern.h"
#include "table.h"
#include "functions.h"

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
    for (int x = 0; x < ROWS; x++)
    {
        for (int y = 0; y < COLS; y++)
        {
            std::cout << "test" << x << " " << y << " " << valuetogrid(x, y, test) << std::endl;
            m_gridData[x][y] = QString::fromStdString(valuetogrid(x, y, test));
        }
    }
    //m_gridData[0][1] = QString::fromStdString(valuetogrid(0, 1));
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
        return m_gridData[index.row()][index.column()];
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
        m_gridData[index.row()][index.column()] = value.toString();
        //for presentation purposes only: build and emit a joined string
        QString result;
        for(int row= 0; row < ROWS; row++)
        {
            for(int col= 0; col < COLS; col++)
            {
                result += m_gridData[row][col] + " ";
            }
        }
        //MyModel->setData( index.row(), index.column(), m_gridData[index.row()][index.column()] );
        emit editCompleted( result );
    }
    return true;
}
Qt::ItemFlags MyModel::flags(const QModelIndex & /*index*/) const
{
    return Qt::ItemIsSelectable |  Qt::ItemIsEditable | Qt::ItemIsEnabled ;
}
