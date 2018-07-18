//
// Created by michal on 2.3.18.
//
#ifndef GRAVITATIONAL_POTENCIAL_SORTBYLATITUDE_H
#define GRAVITATIONAL_POTENCIAL_SORTBYLATITUDE_H

#include <cmath>
#include <atomic>
#include <mutex>

struct Point {
    double x;
    double y;
    double potencial;
};

class sortByLatitude
{
public:
    sortByLatitude();
    bool operator ()(Point &a, Point &b)
    {
        const double diff = 1e-7;

        if(fabs(a.x - b.x) < diff)
        {
            return (a.y < b.y);
        }
        else
        {
            return (a.x < b.x);
        }
    }

};
#endif //GRAVITATIONAL_POTENCIAL_SORTBYLATITUDE_H
