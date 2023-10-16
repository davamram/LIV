#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include "ExtractValues.h"

std::string extractSecondColumn(const std::string& filePath, int col) {
    std::vector<std::string> secondColumnNumbers;
    std::ifstream file(filePath);

    if (!file.is_open()) {
        std::cerr << "Erreur lors de l'ouverture du fichier." << std::endl;
        return "";
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> columns;
        std::istringstream iss(line);
        std::string token;

        while (std::getline(iss, token, '\t')) {
            columns.push_back(token);
        }

        if (columns.size() == 5) {
            secondColumnNumbers.push_back(columns[col]);
        }
    }

    file.close();

    if (secondColumnNumbers.size() >= 16) {
        std::string result;
        for (int i = secondColumnNumbers.size() - 16; i < secondColumnNumbers.size(); ++i) {
            result += secondColumnNumbers[i];
            if (i < secondColumnNumbers.size() - 1) {
                result += ", ";
            }
        }
        return result;
    } else {
        return "";
    }
}

std::vector<double> getValues(int energy, int num) {

    std::string filePath = "/home/amram/Documents/LorentzPhotons/Rivet/LIV/Plots/MadGraph/Reweight/" + std::to_string(energy) + "GeV/TEST_ANALYSIS/d01-x01-y01.dat";
    // num : 2 for values, 3 for errors
    std::string result = extractSecondColumn(filePath, num);

    std::vector<double> numbers;
    std::istringstream ss(result);
    std::string token;

    while (std::getline(ss, token, ',')) {
        double number;
        std::istringstream(token) >> number;
        numbers.push_back(number);
    }

    return numbers;
}