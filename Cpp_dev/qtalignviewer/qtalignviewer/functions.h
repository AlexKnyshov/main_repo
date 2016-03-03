#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

std::map <std::string, std::string> readfastafile(std::string filename); //reading fasta into dict

int fastalen(std::map <std::string, std::string> dict2); //determing the length of the longest sequence

std::string valuetogrid (int row, int col, std::map <std::string, std::string> d); //extract value - deprecated

void fillvector (std::vector< std::vector<std::string> >& vec, std::map <std::string, std::string> dict); //add value to global vector

void updatevector (std::vector< std::vector<std::string> >& vec, std::map <std::string, std::string> dict);

void savefile (std::vector< std::vector<std::string> >& vec, std::string outfile);

#endif // FUNCTIONS_H
