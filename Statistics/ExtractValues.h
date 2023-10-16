#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>

std::string extractSecondColumn(const std::string& filePath, int col);
std::vector<double> getValues(int energy, int num);