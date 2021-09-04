#pragma once

#include <vector>

namespace se {
    // Function type to calculate numerical derivative from
    using mathFunc = double (*)(double);

    // Output file name
    const std::string outputDeviation = "../deviation.csv";
    const std::string outputFirstDerivative = "../first_derivative.csv";
    const std::string outputSecondDerivative = "../second_derivative.csv";

    // Calculate the numerical 1/2 derivative of the given function f, compare it with 
    // the given right oner by its maximum deviation from the original, write the table
    // with step-deviation dependence into file
    // @param f Original function
    // @param f1d First derivative of the original function
    // @param f2d Second derivative of the original function
    // @param section The section to calcualte derivative on
    // @param steps The set of steps to build table from
    void numericalDerivative(mathFunc f, mathFunc f1d, mathFunc f2d, std::pair<double, double> section, const std::vector<double> steps);
}