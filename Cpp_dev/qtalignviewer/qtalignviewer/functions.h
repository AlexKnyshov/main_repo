#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

std::map <std::string, std::string> readfastafile(std::string filename);

int fastalen(std::map <std::string, std::string> dict2);

std::string valuetogrid (int row, int col, std::map <std::string, std::string> d);

//std::map <std::string, std::string> dict;


#endif // FUNCTIONS_H
