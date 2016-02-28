#include <iostream>
#include <fstream>
#include <string>
#include <map>


using namespace std;

string colour_out (string seq)
{
	string new_seq = "";
	for (auto i : seq)
				{
					if (i == 'A')
					{
						new_seq.append("\e[39;41m");
						new_seq.push_back(i);
						new_seq.append("\e[39;49m");
					}
					else if (i == 'G')
					{
						new_seq += "\e[30;43m" + string(1, i) + "\e[39;49m";
					}
					else if (i == 'C')
					{
						new_seq += "\e[30;102m" + string(1, i) + "\e[39;49m";
					}
					else if (i == 'T')
					{
						new_seq += "\e[30;104m" + string(1, i) + "\e[39;49m";
					}
					else
					{
						new_seq += i;
					}
				}
	new_seq += "\e[39;49m";
	return new_seq;
}


int main(int argc, char* argv[]){
	//cout << "\e[31m" << "This is colourful text" << "\e[39m" << endl ;
	cout << "intializing .... \n";
	string infile;
	if (argc == 2) {
		infile = argv[1];
	}
	cout <<"infile given as follows:" << infile << endl;
	map <string, string> dict;
	string current = "0";
	string line;
	string seq;
	int maxsize = 0;
	ifstream myfile;
	cout << "opening the file..." << endl;
	myfile.open(infile);
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
		cout << x.first << string(maxsize - x.first.length(), ' ') << " " << colour_out(x.second.substr(0, 30)) << endl; 
	}
	cout << "total sequences read: " << dict.size() << endl;
	cout << "done" << endl;
	return 0;

}

