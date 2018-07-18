//
// Created by michal on 15.3.18.
//
#ifndef GRAVITATIONAL_POTENCIAL_MATRIX_H
#define GRAVITATIONAL_POTENCIAL_MATRIX_H

#include <algorithm>
#include <vector>

class matrix {

public:
    matrix(){};
    matrix(const int, bool state = true);
    ~matrix() {}
    void setPointers(const int, const int, const std::pair<double*, double*> , bool);
    std::pair<double*, double*> getPointers(int, int, bool = true);
    std::vector<double> getCoefficients(const std::pair<double*, double*>&);
    void setCoefficients(std::vector<double>&, const int, const int, bool);
    std::vector<double> validateCoefficients(std::vector<double>&);
    const std::vector<double>& getVector(const int);
    void addCoefficient(const int, const double);
    void setNorm(const int, const int, const double, bool);
    double getNorm(const int, const int, bool);
    void coeffFromVec(const std::vector<double>&);
    void makePointers(int, const int, const int, bool);

private:

    std::vector<std::pair<double*, double*>> Matrix;
    std::vector<double> coefficients;
    std::vector<std::vector<double>> Coefficients;
    std::vector<double> norms;
    const double* _iter;
    int _width;
    int _size;
    bool _state;
};


#endif //GRAVITATIONAL_POTENCIAL_MATRIX_H
