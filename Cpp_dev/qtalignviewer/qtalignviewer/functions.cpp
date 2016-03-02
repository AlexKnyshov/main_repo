#include "extern.h"
#include "functions.h"

using namespace std;

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
  for ( auto x : dict){
    cout << x.first << string(maxsize - x.first.length(), ' ') << " " << x.second.substr(0, 30) << endl;
  }
  cout << "total sequences read: " << dict.size() << endl;
  cout << "done" << endl;

  return dict;
}

int fastalen(std::map<string, string> dict2)
{
    cout << "fastalen called";
    size_t maxsize = 0;
    for ( auto x : dict2){
        cout << maxsize << " " << x.second.length();
        if (maxsize > x.second.length())
        {
          maxsize = x.second.length();
        }
    }
    return maxsize;
}
