//
// Created by michal on 18.3.18.
//
#ifndef GRAVITATIONAL_POTENCIAL_ALGORITHMS_H
#define GRAVITATIONAL_POTENCIAL_ALGORITHMS_H

#include <array>
#include <atomic>
#include "model.h"
#include "sortByLatitude.h"

#define PI 4*atan(1)

class algorithms {

public:
    static double legenderNorm(const double, const int, const int);
    static double legendre(const double, const int, const int, const int);
    static double sphereFunction(const double, const double, const int, const int, bool = true);
    static double sphereFunctionNorm(const Point&, const int, const int, bool = true);
    static void bilinearInterpolation(Point&, const model&);
    static std::array<std::pair<double , double>, 4> findPoints(Point& gauss, const model& M, const int dist_choice = 0);
    static double earthDistance(const Point&, const Point&);
    static double harvesineDist(const Point&, const Point&);
    static void IDW(Point &, const model&, const int);
    static void nearestNeighbour(Point&, const model&);
    static double euclideanDist(const Point&, const Point&);
};

#endif //GRAVITATIONAL_POTENCIAL_ALGORITHMS_H
