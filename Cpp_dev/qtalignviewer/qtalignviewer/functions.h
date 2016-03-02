#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <iostream>
#include <fstream>
#include <string>
#include <map>


std::map <std::string, std::string> readfastafile(std::string filename);

int fastalen(std::map <std::string, std::string> dict2);

//std::map <std::string, std::string> dict;


#endif // FUNCTIONS_H
