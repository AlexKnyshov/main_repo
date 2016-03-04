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

void fillvector (std::vector< std::vector<std::string> >& vec, std::map <std::string, std::string> dict)
{

    cout << "test" << endl;
    cout << ROWS << " " << COLS << endl;
    int size;

    size = fastalen(dict);
    int rowsize = dict.size();

    for (int i = 0; i < ROWS; i++) {
        vector<string> row; // Create an empty row
        //cout << "row" << i << endl;
        for (int j = 0; j < COLS; j++) {
            row.push_back("?"); // Add an element (column) to the row
            //cout << "line" << j << endl;
        }
        vec.push_back(row); // Add the row to the main vector
    }

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
            row.push_back("?"); // Add an element (column) to the row
            //cout << "line" << j << endl;
        }
        vec[i].swap(row);
        //vec.push_back(row); // Add the row to the main vector
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
                cout << "test" << x << " " << y << " " << vec[x][y] << std::endl;
                if (y == 0)
                {
                    myoutfile << ">"+vec[x][y] << endl;
                }
                else
                {
                    temp  += vec[x][y];
                }
            }
            myoutfile << temp << endl;
        }
}

int translatevector (std::vector< std::vector<std::string> >& vec, int option)
{
    cout << "translatevector() called" << endl;
    int out = COLS;
    map<string, string> transtable = {{"TTT", "F"}, {"TTC", "F"}, {"TTA", "L"}, {"TTG", "i"},
                                                     {"CTT", "L"}, {"CTC", "L"}, {"CTA", "L"}, {"CTG", "i"},
                                                     {"ATT", "I"}, {"ATC", "I"}, {"ATA", "I"}, {"ATG", "i"},
                                                     {"GTT", "V"}, {"GTC", "V"}, {"GTA", "V"}, {"GTG", "V"},
                                                     {"TCT", "S"}, {"TCC", "S"}, {"TCA", "S"}, {"TCG", "S"},
                                                     {"CCT", "P"}, {"CCC", "P"}, {"CCA", "P"}, {"CCG", "P"},
                                                     {"ACT", "T"}, {"ACC", "T"}, {"ACA", "T"}, {"ACG", "T"},
                                                     {"GCT", "A"}, {"GCC", "A"}, {"GCA", "A"}, {"GCG", "A"},
                                                     {"TAT", "Y"}, {"TAC", "Y"}, {"TAA", "*"}, {"TAG", "*"},
                                                     {"CAT", "H"}, {"CAC", "H"}, {"CAA", "Q"}, {"CAG", "Q"},
                                                     {"AAT", "N"}, {"AAC", "N"}, {"AAA", "K"}, {"AAG", "K"},
                                                     {"GAT", "D"}, {"GAC", "D"}, {"GAA", "E"}, {"GAG", "E"},
                                                     {"TGT", "C"}, {"TGC", "C"}, {"TGA", "*"}, {"TGG", "W"},
                                                     {"CGT", "R"}, {"CGC", "R"}, {"CGA", "R"}, {"CGG", "R"},
                                                     {"AGT", "S"}, {"AGC", "S"}, {"AGA", "R"}, {"AGG", "R"},
                                                     {"GGT", "G"}, {"GGC", "G"}, {"GGA", "G"}, {"GGG", "G"}};
    vector< vector<string> > tempvec;
    int prot_col = 0;
    for (int x = 0; x < ROWS; x++)
        {
            int temp_col = 0;
            string temp = "";
            string codontemp = "";
            vector<string> row;
            for (int y = 0; y < COLS-option; y++)
            {
                //cout << "test" << x << " " << y << " " << vec[x][y] << std::endl;
                if (y == 0)
                {
                    //cout << ">"+vec[x][y] << endl; //title
                    row.push_back(vec[x][y]);
                }
                else
                {
                    if (y % 3 == 0)
                    {
                        codontemp  += vec[x][y+option];//creating seq
                        for (auto & c: codontemp) c = toupper(c);
                        //cout << "yes " << y << transtable[codontemp] << endl;
                        if ( transtable.find(codontemp) == transtable.end() ) {
                          row.push_back("-");
                        } else {
                          row.push_back(transtable[codontemp]);
                        }

                        temp_col += 1;
                        codontemp = "";

                    }
                    else
                    {

                        codontemp  += vec[x][y+option];//creating seq
                        //cout << "else " << y << codontemp << endl;
                    }
                }
            }
            //cout << temp << endl;
            tempvec.push_back(row);
            if (temp_col > prot_col)
            {
                prot_col = temp_col;
            }

            //myoutfile << temp << endl;
        }
        //cout << prot_col;
    //vector output
        for (int i = 0; i < ROWS; i++) {
            for (int j = 0; j < prot_col; j++) {
                vec[i][j] = tempvec[i][j];
            }
            vec[i].erase(vec[i].begin()+prot_col, vec[i].end());
        }
        cout << COLS << " cols" << endl;
        COLS = prot_col;
        cout << COLS << " cols" << endl;
        return out;
}
