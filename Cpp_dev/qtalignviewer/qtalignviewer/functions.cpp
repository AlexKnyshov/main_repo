#include "extern.h"
#include "functions.h"

using namespace std;

//int COLS;
//int ROWS;

std::map <std::string, std::string> readfastafile(string filename)
{
  std::map <std::string, std::string> dict;
  cout << "intializing .... \n";

  cout <<"infile given as follows:" << filename << endl;

  string current = "0";
  string line;
  string seq;
  int maxsize = 0;
  ifstream myfile;
  cout << "opening the file..." << endl;
  myfile.open(filename);
  if (myfile.is_open()){
    cout << "file opened, reading..." << endl;
    while (getline (myfile, line))
    {
      line.erase(line.find_last_not_of(" \n\r\t")+1);
      if (line[0] == '>')
      {
        if (current != "0")
        {
          dict[current.substr(1, current.length())]=seq;
        }
        seq = "";
        current = line;
        if (current.length() > maxsize)
        {
          maxsize = current.length();
        }
      }
      else
      {
        seq += line;
      }
    }
    dict[current.substr(1, current.length())]=seq;
    myfile.close();
  }
  else
  {
    cout << "cannot open the file" << endl;
  }
//  for ( auto x : dict){
//    cout << x.first << string(maxsize - x.first.length(), ' ') << " " << x.second.substr(0, 30) << endl;
//  }
  cout << "total sequences read: " << dict.size() << endl;
  cout << "done" << endl;

  return dict;
}

int fastalen(std::map<string, string> dict2)
{
    //cout << "fastalen called";
    size_t maxsize = 0;
    for ( auto x : dict2){
        //cout << "maxs " << maxsize << " x.sec " << x.second.length() << endl;
        if (maxsize < x.second.length())
        {
          maxsize = x.second.length();
        }
    }
    return int(maxsize);
}

//string valuetogrid (int row, int col, map <string, string> d)
//{

//    //cout << vec[0][0];
//    return vec[row][col];
//    //int *pArray = new int[rowsize];
//}
void fillvector (std::vector< std::vector<std::string> >& vec, std::map <std::string, std::string> dict)
{
    //int ROWS, COLS;
    cout << "test" << endl;
    cout << ROWS << " " << COLS << endl;
    int size;
    //map <string, string> test = readfastafile("winglesscat.fas");
    size = fastalen(dict);
    int rowsize = dict.size();
    //vector< vector<string> > vec;
    for (int i = 0; i < ROWS; i++) {
        vector<string> row; // Create an empty row
        //cout << "row" << i << endl;
        for (int j = 0; j < COLS; j++) {
            row.push_back("Temp"); // Add an element (column) to the row
            //cout << "line" << j << endl;
        }
        vec.push_back(row); // Add the row to the main vector
    }
//    map<string, string>::iterator it;
//    int rowcount = 0;
//    for ( it = dict.begin(); it != dict.end(); it++ )
//    {
//        vec[rowcount][0] = it->first;
//        //cout << "row" << rowcount << endl;
//        for (int z = 1; z < COLS; z++)//it->second.length(); z++)
//        {
//            vec[rowcount][z] = it->second.substr(z, 1);
//            //cout << "line" << z << endl;
//        }
//        rowcount += 1;
//    }
}

void updatevector (vector< vector<string> >& vec, map <string, string> dict)
{
    int size = fastalen(dict);
    int rowsize = dict.size();
    COLS = size;
    cout << COLS << endl;
    ROWS = rowsize;
    cout << ROWS << endl;

    for (int i = 0; i < ROWS; i++) {
        vector<string> row; // Create an empty row
        //cout << "row" << i << endl;
        for (int j = 0; j < COLS; j++) {
            row.push_back("Temp"); // Add an element (column) to the row
            //cout << "line" << j << endl;
        }
        vec.push_back(row); // Add the row to the main vector
    }



    map<string, string>::iterator it;
    int rowcount = 0;
    for ( it = dict.begin(); it != dict.end(); it++ )
    {
        vec[rowcount][0] = it->first;
        cout << "row" << rowcount << " " << it->second.substr(0, it->second.length()) << endl;
        for (int z = 0; z < (it->second.length())-1; z++)
        {
            vec[rowcount][z+1] = it->second.substr(z, 1);
            //cout << "line" << z << " " << it->second.substr(z, 1) << endl;
        }
        rowcount += 1;
    }
}

void savefile (vector< vector<string> >& vec, string outfile)
{
    cout << "outfile " << outfile << endl;
    ofstream myoutfile;
    myoutfile.open (outfile);
    for (int x = 0; x < ROWS; x++)
        {
            string temp = "";
            for (int y = 0; y < COLS; y++)
            {
                cout << "test" << x << " " << y << " " << vector1[x][y] << std::endl;
                if (y == 0)
                {
                    myoutfile << ">"+vector1[x][y] << endl;
                }
                else
                {
                    temp  += vector1[x][y];
                }
            }
            myoutfile << temp << endl;
        }
}

