//
// Created by michal on 15.3.18.
#ifndef GRAVITATIONAL_POTENCIAL_LOCALORTHOGONALSYSTEM_H
#define GRAVITATIONAL_POTENCIAL_LOCALORTHOGONALSYSTEM_H

#include "matrix.h"
#include "sortByLatitude.h"
#include <array>
#include "algorithms.h"

class localOrthogonalSystem {

public:

    localOrthogonalSystem(const int, const std::array<double, 4> &);
    localOrthogonalSystem(const int, const std::array<double, 4> &, const std::string&);
    double getOrthoFunc(const int, const int, const Point &, bool = true);
    std::vector<double>getVector(const int, const int, bool);
    ~localOrthogonalSystem(){delete M;}

private:

    void coeffToTxt();
    void writeCoeff(std::ofstream&, const std::vector<double>&, const double, const double);
    void computeValues(int);
    double computeCoefficients(int, int, int, int, bool, bool);
    double getOrthoFunc(int, int, double, double, bool = true);
    double gaussIntegration(const std::array<double , 4> &, const int, const int, const int, const int, bool, bool);
    double gaussIntegration(const std::array<double , 4> &, const int, const int, bool);
    double gaussINT(const std::array<double , 4> &, const int, const int, const int, const int, bool, bool);
    double gaussINT(const std::array<double , 4> &, const int, const int, bool);

    matrix* M = nullptr;
    const int _max;
    std::array<double , 4> boundary;
    std::vector<double> duration_sys;
};


#endif //GRAVITATIONAL_POTENCIAL_LOCALORTHOGONALSYSTEM_H
