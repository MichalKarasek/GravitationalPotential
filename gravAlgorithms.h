//
// Created by michal on 17.7.18.
//
#ifndef GRAVITATIONAL_POTENCIAL_GRAVALGORITHMS_H
#define GRAVITATIONAL_POTENCIAL_GRAVALGORITHMS_H

#include "localOrthogonalSystem.h"
#include <iostream>
#include <mutex>

class gravAlgorithms
{
public:
    static void gravitationalPotential(model&, localOrthogonalSystem&, const Point&, bool, int, const std::array<double, 4>&);

private:
    static double gaussINT(const int, const int, int, model &, localOrthogonalSystem&, bool, bool, const int,
                           const std::array<double, 4>& arr={-PI/2, PI/2, 0, 2*PI}, const int idwpar = 2);
    static double gauss(const std::array<double, 4> &, model &, int, int, int, bool, bool, const int, const int, localOrthogonalSystem& );
    static inline double function(const Point &, const int, const int, const int, bool choice);
    static inline double functionLocal(const Point &, const int, const int, const int, bool, localOrthogonalSystem&);
    static void compute(const int, const int, model&, localOrthogonalSystem& ,const Point&, bool, int, const std::array<double, 4>&);
};

#endif //GRAVITATIONAL_POTENCIAL_GRAVALGORITHMS_H
