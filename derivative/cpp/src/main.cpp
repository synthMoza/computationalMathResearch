#include <iostream>
#include <vector>

#include "derivative.hpp"

// Epsilon for comparing double values
const double eps = 1e-7;

// Get the section from the standart input with help info written to standart output
// @return Pair of doubles that represents a correct section [a;b]
std::pair<double, double> getSection() {
    std::pair<double, double> result;

    std::cout << "Enter the section to take derivate on: ";
    std::cin >> result.first >> result.second;
    while (std::cin.fail() || result.first >= result.second) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "Invalid section! Please try again: ";
        std::cin >> result.first >> result.second;
    }

    return result;
}

// Get the step from the standart input with help info written to standart output
// @return Steps for building table function
std::vector<double> getSteps() {
    std::vector<double> result;

    double start = 0, end = 0;
    int amount = 0;

    // Start step
    std::cout << "Enter the start step: ";
    std::cin >> start;
    while (std::cin.fail() || start < eps) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "Invalid step! Please try again: ";
        std::cin >> start;
    }

    // End step
    std::cout << "Enter the end step: ";
    std::cin >> end;
    while (std::cin.fail() || end < eps) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "Invalid step! Please try again: ";
        std::cin >> end;
    }

    // Amount
    std::cout << "Enter the amount of steps: ";
    std::cin >> amount;
    while (std::cin.fail() || amount < eps) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "Invalid step! Please try again: ";
        std::cin >> amount;
    }

    // Fil steps
    result.resize(amount);

    double step = (end - start) / (amount - 1);
    for (int i = 0; i < amount; ++i)
        result.at(i) = start + i * step;

    return result;
}

// Target function and its first and second derivative for calculating numerical derivative
// and compare results

double targetFunction(double x) {
    return sin(x) * x;
}

double targetFirstDerivative(double x) {
    return x * cos(x) + sin(x);
}

double targetSecondDerivative(double x) {
    return 2 * cos(x) - x * sin(x); 
}

using namespace se;

int main() {
    // Get the section to take derivative on
    auto section = getSection();
    // Get the step
    auto steps = getSteps();

    numericalDerivative(targetFunction, targetFirstDerivative, targetSecondDerivative, section, steps);
}