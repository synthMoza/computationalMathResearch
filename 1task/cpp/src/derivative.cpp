#include <iostream>
#include <fstream>
#include <string>

#include "derivative.hpp"

using namespace se;

enum order {
    NONE = -1, FIRST, SECOND, COUNT
};

static double secondNumericalDerivative(mathFunc f, mathFunc f2d, double step, std::pair<double, double> section, bool output = false) {
    double r = 0;

    std::vector<double> x;
    std::vector<double> derivative;

    // Fill table function
    std::size_t size = (section.second - section.first) / step + 1;

    x.resize(size);
    for (std::size_t i = 0; i < size; ++i)
        x.at(i) = section.first + step * i;

    derivative.resize(size);
    
    derivative.at(0) = (2 * f(x.at(0)) - 5 * f(x.at(1)) + 4 * f(x.at(2)) - f(x.at(3))) / step / step;
    for (std::size_t i = 1; i < size - 1; ++i)
        derivative.at(i) = (f(x.at(i-1)) - 2 * f(x.at(i)) + f(x.at(i+1))) / step / step;
    derivative.at(size-1) = (2 * f(x.at(size-1))  - 5 * f(x.at(size-2)) + 4 * f(x.at(size-3)) - f(x.at(size - 4))) / step / step;

    double tmp = 0;
    for (std::size_t i = 0; i < size; ++i) {
        tmp = abs(derivative.at(i) - f2d(x.at(i)));
        if (tmp > r)
            r = tmp;
    }

    if (output) {
        // Put derivative into file
        std::ofstream file;
        file.open(outputSecondDerivative);
            
        for (std::size_t i = 0; i < size; ++i)
            file << x.at(i) << "\t" << derivative.at(i) << std::endl; 
    }

    return r;
}

static double firstNumericalDerivative(mathFunc f, mathFunc f1d, double step, std::pair<double, double> section, order o, bool output = false) {
    double r = 0;

    std::vector<double> x;
    std::vector<double> derivative;
    
    // Fill table function
    std::size_t size = (section.second - section.first) / step + 1;
    
    x.resize(size);
    for (std::size_t i = 0; i < size; ++i)
        x.at(i) = section.first + step * i;

    switch (o) {
        case FIRST:
            // Calculate derivative
            derivative.resize(size);
            for (std::size_t i = 0; i < size - 1; ++i)
                derivative.at(i) = (f(x.at(i + 1)) - f(x.at(i))) / step;

            derivative.at(size - 1) = (f(x.at(size - 1)) - f(x.at(size - 2))) / step;

            break;
        case SECOND:
            // Calculate deviation
            derivative.resize(size);
            for (std::size_t i = 0; i < size - 2; ++i)
                derivative.at(i) = (-3 * f(x.at(i)) + 4 * f(x.at(i + 1)) - f(x.at(i + 2))) / (2 * step);

            for (std::size_t i = size - 2; i < size; ++i)
                derivative.at(i) = (3 * f(x.at(i)) - 4 * f(x.at(i - 1)) + f(x.at(i - 2))) / (2 * step);

            break;
        default:
            throw std::runtime_error("Unknown first derivative deviation order!");
    }

    // Calculate deviation
    double tmp = 0;
    for (std::size_t i = 0; i < size; ++i) {
        tmp = abs(derivative.at(i) - f1d(x.at(i)));
        if (tmp > r)
            r = tmp;
    }

    if (output) {
        // Put derivative into file
        std::ofstream file;
        file.open(outputFirstDerivative);

        for (std::size_t i = 0; i < size; ++i)
            file << x.at(i) << "\t" << derivative.at(i) << std::endl; 
    }

    return r;
}

void se::numericalDerivative(mathFunc f, mathFunc f1d, mathFunc f2d, std::pair<double, double> section, const std::vector<double> steps) {
    std::ofstream file;
    file.open(outputDeviation);
    
    // Iterate through all given steps
    for (std::size_t i = 0; i < steps.size() - 1; ++i) {
        file << steps.at(i) << "\t";
        file << firstNumericalDerivative(f, f1d, steps.at(i), section, FIRST) << "\t";
        file << firstNumericalDerivative(f, f1d, steps.at(i), section, SECOND) << "\t";
        file << secondNumericalDerivative(f, f2d, steps.at(i), section) << std::endl;
    }

    // Output the last one to the file
    file << steps.back() << "\t";
    file << firstNumericalDerivative(f, f1d, steps.back(), section, FIRST) << "\t";
    file << firstNumericalDerivative(f, f1d, steps.back(), section, SECOND, true) << "\t";
    file << secondNumericalDerivative(f, f2d, steps.back(), section, true) << std::endl;

    file.close();
}