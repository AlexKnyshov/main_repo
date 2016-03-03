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

void translatevector (std::vector< std::vector<std::string> >& vec); //translate and update seq

#endif // FUNCTIONS_H
