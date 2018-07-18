//
// Created by michal on 1.3.18.
//
#ifndef GRAVITATIONAL_POTENCIAL_MODEL_H
#define GRAVITATIONAL_POTENCIAL_MODEL_H

#include <vector>
#include <string>
#include <array>
#include <mutex>
#include "sortByLatitude.h"

class model
{
public:

    model(const std::string &, const int, const int);
    model(const std::string &, const int , const int, const int, const double latitude = -(4*atan(1))/2, const double longitude = 0);
    ~model(){}
    std::array<Point,4> findPoints(const Point &)const;
    std::pair<double , int> getPair(int pos)const{ return  X_coordinates.at(pos);}
    void setCoefficient(const int&, const int&, const double&, bool);
    double getCoefficients(const int&, const int&, bool)const;
    int getMax()const{ return  _width-1;}
    void addPotential(const double&, const int&, const int&);
    double getPotential()const{ return gravity_potential;}
    double getPotentialDiff(int pos)const{ return potential_diff.at(pos);}
    void writeToTxt(std::ofstream&);

private:

    void processDiff();
    std::vector<Point> getVector(const int)const;
    Point* getPoints(std::vector<Point> &, const double)const;
    std::vector<Point> Model;
    std::vector<std::pair<double, int >> X_coordinates;
    std::vector<double> coefficients;
    std::vector<double> potential_diff;
    int _width;
    double gravity_potential;
    std::mutex lock;
};
#endif //GRAVITATIONAL_POTENCIAL_MODEL_H
