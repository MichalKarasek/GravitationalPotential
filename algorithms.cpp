//
// Created by michal on 18.3.18.
//

#include <thread>
#include <iostream>
#include "algorithms.h"

/**
 * @param x - latitude
 * @param a - order of Legendre associated polynomials
 * @param b - degree of Legendre associated polynomials
 * @return value of normed Legendre associated polynomial of order a and degree b for latitude x
 */
double algorithms::legenderNorm(const double x, const int a, const int b)
{
    // Compute value of Legendre associated functions in point x
    int size = 41;
    int position = (a*(a+1))/2 + b;

    std::vector<double> legender = std::vector<double>(static_cast<size_t>((size*(size+1)))/2);

    legender.at(0) = 1;
    legender.at(1) = std::sqrt(3)*sin(x);
    legender.at(2) = std::sqrt(3)*cos(x);
    legender.at(4) = std::sqrt(15)*sin(x)*cos(x);

    // Cycle filling main diagonal
    for(int i=2; i<= size-1; ++i)
    {
        int pos = (i*(i+1))/2 + i;
        int pos_reverse = (i*(i-1))/2 + i-1;
        legender.at((u_int)pos)=sqrt((2.0*i+1)/(2.0*i))*cos(x)*legender.at((u_int)pos_reverse);
        if(pos == position)
        {
            return legender.at((u_int)pos);
        }
    }

    // Cycle filling first minor diagonal
    for(int i=3; i<= size-1; ++i)
    {
        int pos = (i*(i+1))/2 + i-1;
        int pos_reverse = (i*(i-1))/2 + i-2;
        legender.at((u_int)pos) = sqrt((2.0*i+1)/(2.0*i-2))*cos(x)*legender.at((u_int)pos_reverse);
        if(pos == position)
        {
            return legender.at((u_int)pos);
        }
    }

    // Cycle for remaining elements of lower triangular matrix
    for(int i=2; i<=size-1; ++i)
    {
        for(int j=0; j<i-1; ++j)
        {
            double w=sqrt(((2.0*i+1)*(2.0*i-1))/((i+j)*(i-j)));
            double w1=sqrt(((2.0*i-1)*(2.0*i-3))/((i+j-1)*(i-j-1)));
            int pos = (i*(i+1))/2 + j;
            int pos_reverse = (i*(i-1))/2 + j;
            int pos_reserve2 = ((i-2)*(i-1))/2 + j;
            legender.at((u_int)pos)=w*(sin(x)*legender.at((u_int)pos_reverse)-(1.0/w1)*legender.at((u_int)pos_reserve2));
            if(pos == position)
            {
                return legender.at((u_int)pos);
            }
        }
    }
    return legender.at((u_int)position);
}
/**
* @brief Compute value of spheric function
* @param x - latitude
* @param y - longitude
* @param m - degree of LAF
* @param n - order of LAF (n <= m)
* @param choice - compute first or second Spheric function
* @return value in defined longitude and latitude of the spheric function
*/
double algorithms::sphereFunction(const double x, const double y, const int m, const int n, bool choice)
{
    return legenderNorm(x,m,n) * (choice ? cos : sin)(n * y);
}
/**
 * @brief compute value of gravitational potential in gauss point using Bilinear interpolation method
 * @param p - point in which is gravitational potential interpolated
 * @param M - model of gravitational potential
 */
void algorithms::bilinearInterpolation(Point& p, const model& M)
{
    double a0{};
    double a1{};
    double a2{};
    double a3{};

    std::array<Point, 4> points = M.findPoints(p);

    double x1 = points.front().x;
    double y1 = points.front().y;
    double x2 = points.back().x;
    double y2 = points.back().y;

    a0 = (points.front().potencial*x2*y2)/((x1-x2)*(y1-y2)) + (points.at(1).potencial*x2*y1)/((x1-x2)*(y2-y1)) + (points.at(2).potencial*x1*y2)/((x1-x2)*(y2-y1)) + (points.back().potencial*x1*y1)/((x1-x2)*(y1-y2));
    a1 = (points.front().potencial*y2)/((x1-x2)*(y2-y1)) + (points.at(1).potencial*y1)/((x1-x2)*(y1-y2)) + (points.at(2).potencial*y2)/((x1-x2)*(y1-y2)) + (points.back().potencial*y1)/((x1-x2)*(y2-y1));
    a2 = (points.front().potencial*x2)/((x1-x2)*(y2-y1)) + (points.at(1).potencial*x2)/((x1-x2)*(y1-y2)) + (points.at(2).potencial*x1)/((x1-x2)*(y1-y2)) + (points.back().potencial*x1)/((x1-x2)*(y2-y1));
    a3 = (points.front().potencial)/((x1-x2)*(y1-y2)) + (points.at(1).potencial)/((x1-x2)*(y2-y1)) + (points.at(2).potencial)/((x1-x2)*(y2-y1)) + (points.back().potencial)/((x1-x2)*(y1-y2));

    p.potencial = a0 + a1*p.x + a2*p.y + a3*p.x*p.y;
}

/**
 * @param gauss - point for which are nearest points seek
 * @param M - model of gravitational potential
 * @param dist_choice - algorithm for distance computation
 * @return array of four nearest points
 */
std::array<std::pair<double, double>, 4> algorithms::findPoints(Point& gauss, const model& Model, const int dist_choice)
{
    std::array<std::pair<double, double>, 4> result;
    std::array<Point, 4> points;

    points = Model.findPoints(gauss);

    for (size_t i=0; i < points.size(); ++i)
    {
        switch(dist_choice)
        {
            case 1:
            {
                result.at(i).first = earthDistance(gauss, points.at(i));
                break;
            }
            case 2:
            {
                result.at(i).first = euclideanDist(gauss, points.at(i));
                break;
            }
            default:
            {
                result.at(i).first = harvesineDist(gauss, points.at(i));
                break;
            }
        }
        result.at(i).second = points.at(i).potencial;
    }
    return result;
}
/**
 * @param a - Point structure
 * @param b - Point structure
 * @return sphere distance of points a and b in meters
 */
double algorithms::earthDistance(const Point& a, const Point& b)
{
    double radius = 6371000.79;
    double leng = std::abs(a.y - b.y);
    double angle = acos(sin(a.x) * sin(b.x) + cos(a.x) * cos(b.x) * cos(leng));

    return  radius * angle;
}
/**
 * @brief compute value of gravitational potential in gauss point using IDW interpolation method
 * @param gauss - point with unknown gravitational potential
 * @param M - gravitational potential model
 * @param parameter - power parameter, by increasing its value the nearest neighbour multiplied impact
 *        on final gravitational potential value
 */
void algorithms::IDW(Point &gauss, const model &M, const int parameter)
{
    double dividend{};
    double aliquot{};
    double weight{};
    std::array<std::pair<double, double>, 4> arr;
    arr = algorithms::findPoints(gauss, M);

    for(auto &v:arr)
    {
        weight = 1/pow(v.first, parameter);
        dividend += weight*v.second;
        aliquot += weight;
    }
    gauss.potencial = dividend/aliquot;
}
/**
 * @brief compute - value of gravitational potential in point gauss using nearest
 *        neighbour interpolation method
 * @param gauss - point with unknown gravitational potential
 * @param M - gravitational potential model
 */
void algorithms::nearestNeighbour(Point &gauss, const model &M)
{
    std::array<std::pair<double, double>, 4> arr;
    arr = algorithms::findPoints(gauss, M);
    double distance = 1e6;

    for(auto &v:arr)
    {
        if(v.first < distance)
        {
            distance = v.first;
            gauss.potencial = v.second;
        }
    }
}
/**
 * @brief on short distance in meters and smaller it is more accurate approach of
 *        calculation then earthDistance algorithm
 * @param a - Point structure
 * @param b - Point structure
 * @return sphere distance of points a and b in meters
 */
double algorithms::harvesineDist(const Point& a, const Point& b)
{
    double radius = 637100.79;
    double deltaFi = std::abs(a.x - b.x);
    double deltaLa = std::abs(a.y - b.y);

    double temp = std::sin(deltaFi/2)*std::sin(deltaFi/2) + std::cos(a.x) * std::cos(b.x) * std::sin(deltaLa/2)*std::sin(deltaLa/2);
    temp = 2 * std::atan2(std::sqrt(temp), std::sqrt(1-temp));

    return temp * radius;
}
/**
 * @param a - Point
 * @param b - Point
 * @return Euclidean distance between a and b
 */
double algorithms::euclideanDist(const Point& a, const Point& b)
{
    return std::sqrt((a.y-b.y)*(a.y-b.y)+(a.x-b.x)*(a.x-b.x));
}
/**
 * @param x - latitude
 * @param a - order of Legendre associated polynomials more info (https://en.wikipedia.org/wiki/Associated_Legendre_polynomials)
 * @param b - degree of Legendre associated polynomials
 * @return value of Legendre associated polynomial of order a and degree b for latitude x
 */
double algorithms::legendre(const double x, const int a, const int b, const int size)
{
    int position = (a*(a+1))/2 + b;

    std::vector<double> legender = std::vector<double>(static_cast<size_t>((size*(size+1)))/2);

    legender.at(0) = 1;
    legender.at(1) = sin(x);

    // Cycle filling main diagonal
    for(int i=1; i<= size-1; ++i)
    {
        int pos = (i*(i+1))/2 + i;
        int pos_reverse = (i*(i-1))/2 + i-1;
        legender.at((u_int)pos)=legender.at((u_int)pos_reverse)*(2.0*i-1.0)*cos(x);
        if(pos == position)
        {
            return legender.at((u_int)pos);
        }
    }
    // Cycle filling first minor diagonal
    for(int i=2; i<= size-1; i++)
    {
        int pos = (i*(i+1))/2 + i-1;
        int pos_reverse = (i*(i-1))/2 + i-2;
        legender.at((u_int)pos) = (2.0*i-1)*cos(x)*legender.at((u_int)pos_reverse);
        if(pos == position)
        {
            return legender.at((u_int)pos);
        }
    }

    // Cycle for remaining elements of lower triangular matrix
    for(int i=2; i<=size-1; i++)
    {
        for(int j=0; j<i-1; j++)
        {
            int pos = (i*(i+1))/2 + j;
            int pos_reverse = (i*(i-1))/2 + j;
            int pos_reserve2 = ((i-2)*(i-1))/2 + j;
            legender.at((u_int)pos)=((2.0*i-1)/(i-j))*sin(x)*legender.at((u_int)pos_reverse)-((i+j-1.0)/(i-j))*legender.at((u_int)pos_reserve2);
            if(pos == position)
            {
                return legender.at((u_int)pos);
            }
        }
    }
    return legender.at((u_int)position);
}
/**
 * @param p - Point in which is sphere function value being evaluated
 * @param m - degree of Legendre associated polynomials
 * @param n - order of Legendre associated polynomials
 * @param choice - specify shape of sphere function
 * @return value of sphere function in point p
 */
double algorithms::sphereFunctionNorm(const Point& p, const int m, const int n, bool choice)
{
    return legenderNorm(p.x,m,n) * (choice ? cos : sin)(n * p.y);
}

