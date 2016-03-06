#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <cctype>

std::map <std::string, std::string> readfastafile(std::string filename); //reading fasta into dict

int fastalen(std::map <std::string, std::string> dict2); //determing the length of the longest sequence

void fillvector (std::vector< std::vector<std::string> >& vec, std::map <std::string, std::string> dict); //add value to global vector

void updatevector (std::vector< std::vector<std::string> >& vec, std::map <std::string, std::string> dict); //reload file

void savefile (std::vector< std::vector<std::string> >& vec, std::string outfile); // save file

int translatevector(std::vector< std::vector<std::string> >& vec, int option); //translate and update seq

void shiftnucl (std::vector< std::vector<std::string> >& vec, int xpos, int ypos);

void deletenucl (std::vector< std::vector<std::string> >& vec, int xpos, int ypos);

void deletecol (std::vector< std::vector<std::string> >& vec, int ypos, int count);

void dnd_row (std::vector< std::vector<std::string> >& vec, int row, int count);

void insert_dnd_row (int row, int count);

#endif // FUNCTIONS_H
