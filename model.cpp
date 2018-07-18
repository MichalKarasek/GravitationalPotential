//
// Created by michal on 1.3.18.
//
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include "model.h"

/**
 * @param file_path - path to txt file containing data
 * @param size - number of points in model
 * @param max - maximum partial sum of Fourier series
 */
model::model(const std::string &file_path, const int size, const int max):_width(max+1)
{
    std::ifstream file;
    std::string line;
    std::string::size_type sz;
    std::pair<double, int> value;
    coefficients = std::vector<double>(static_cast<size_t>(_width*_width));
    potential_diff = std::vector<double>(static_cast<size_t>(_width*(_width+1)));
    double temp_x{};
    int counter{};
    Point p;

    file.open(file_path);

    if(!file.is_open())
    {
        std::cerr<<"Fail to open data model file"<<'\n';
        exit(-1);
    }
    else
    {
        std::cout<<"Data model file open successfully"<<std::endl;
    }

    Model.reserve(static_cast<size_t>(size));

    while(getline(file, line))
    {
        p.x = std::stod(line, &sz);

        if(temp_x != p.x)
        {
            value.first = p.x;
            value.second = counter;
            X_coordinates.push_back(value);
        }
        line = line.substr(sz);
        p.y = std::stod(line, &sz);
        p.potencial = std::stod(line.substr(sz));

        Model.push_back(p);
        temp_x = p.x;
        ++counter;
    }
    file.close();

    //Sort Model by latitude. Value with same latitude sort by longitude. Ascending
    std::sort(Model.begin(), Model.end(), sortByLatitude());
}
/**
 * @brief Constructor that accept gravitational potential model as matrix (size_1 x size_2) of
 *        gravitational potential values. Constructor treat matrix as regular grid where distance between
 *        points in latitude and longitude is constant. It is necessary to put in coordinates of point at position [0,0].
 *        Usually less used way of data file.
 * @param file_path - path to txt file containing data
 * @param size_1 - number of intervals in latitude
 * @param size_2 - number of intervals in longitude
 * @param latitude - latitude of the first point in grid
 * @param longitude - longitude of the first point in grid
 */
model::model(const std::string &file_path, const int size_1, const int size_2, int max, const double latitude, const double longitude):_width(max+1)
{
    std::ifstream file;
    std::string line;
    std::string::size_type sz;
    double PI = 4*atan(1);
    coefficients = std::vector<double>(static_cast<size_t>(_width*_width));
    potential_diff = std::vector<double>(static_cast<size_t>(_width*(_width+1)));
    Point p;

    file.open(file_path);

    if(!file.is_open())
    {
        std::cerr<<"Fail to open data model file"<<'\n';
        exit(-1);
    }
    else
    {
        std::cout<<"Data model file open successfully"<<std::endl;
    }

    Model.reserve(size_1*size_2);

    double latitude_step = PI/(size_1-1);
    double longitude_step = 2*PI/size_2;
    double longt = longitude;
    double latit = latitude;
    int line_number = 0;
    std::pair<double, int> values;

    while(getline(file, line))
    {
        for(size_t i = 0; i < size_2; ++i)
        {
            p.x = latit;
            p.y = longt;
            p.potencial = std::stod(line, &sz);
            line = line.substr(sz);

            Model.push_back(p);
            longt += longitude_step;
        }

        values.first = p.x;
        values.second = line_number;
        X_coordinates.push_back(values);
        latit += latitude_step;
        longt = longitude;
        line_number += size_2;
    }
    file.close();
}
//TODO add option to find N - nearest points to next interpolation ???
/**
 * @brief function that will find four closest points in model from which the new value of
 *        gravitational potential will be interpolated
 * @param a - Point of interest with unknown gravitational potential
 * @return four nearest points from point a
 */
std::array<Point, 4> model::findPoints(const Point &a)const
{
    std::array<Point, 4> result;
    std::vector<Point> upper_latitude;
    std::vector<Point> lower_latitude;
    size_t size = X_coordinates.size();
    bool split = true;
    size_t half = size/2;
    size_t interval = size - half;
    const Point *point = nullptr;

    if(std::abs(a.x) > 4*atan(1)/2 || a.y > 8*atan(1))
        std::cerr<<"Out of bounds"<<'\n';

    // Create two vectors of points such as point A latitude lies between these two latitudes
    while(split)
    {
        if(a.x < X_coordinates.at(half).first)
        {
            if(a.x > X_coordinates.at(half-1).first)
            {
                upper_latitude = getVector(X_coordinates.at(half).second);
                lower_latitude = getVector(X_coordinates.at(half-1).second);

                // Search for proper points and add them to result array
                point = getPoints(upper_latitude, a.y);
                result.at(0) = *point;
                result.at(1) = *++point;

                point = getPoints(lower_latitude, a.y);
                result.at(2) = *point;
                result.at(3) = *++point;

                split = false;
            }
            // In edge of intervals for latitude and longitude check against
            // infinity loop needs to be done.
            if(interval/2 == 0 && (half-2) >= 0)
            {
                --half;
            }
            else
            {
                half -= interval/2;
                interval /= 2;
            }
        }
        else
        {
            if(a.x < X_coordinates.at(half+1).first)
            {
                upper_latitude = getVector(X_coordinates.at(half).second);
                lower_latitude = getVector(X_coordinates.at(half+1).second);

                // Search for proper points and add them to result array
                point = getPoints(upper_latitude, a.y);
                result.at(0) = *point;
                result.at(1) = *++point;

                point = getPoints(lower_latitude, a.y);
                result.at(2) = *point;
                result.at(3) = *++point;

                split = false;
            }
            if(interval/2 == 0 && (half+2) <= size)
            {
                ++half;
            }
            else
            {
                half += interval/2;
                interval /= 2;
            }
        }
    }
    return result;
}
/**
 * @brief
 * @param a
 * @return
 */
std::vector<Point> model::getVector(const int a)const
{
    std::vector<Point> points;
    Point p;
    int b = a + 1;

    p = Model.at(a);
    points.push_back(p);
    // add points to vector until they have same value or maximum size of
    // model is reached
    while(b < Model.size() && Model.at(a).x == Model.at(b).x)
    {
        p = Model.at(b);
        points.push_back(p);
        ++b;
    }
    // Longitude starting and ending point are identical. Due to next calculation first point needs to
    // be added to vector as last one with longitude 2*PI==(0)
    p = {
            points.front().x,
            8*atan(1),
            points.front().potencial
    };
    points.push_back(p);

    return  points;
}
/**
 *
 * @param vec
 * @param y
 * @return
 */
Point* model::getPoints(std::vector<Point> &vec, const double y)const
{
    int size = vec.size();
    bool split = true;
    int half = size/2;
    int interval = size - half;

    // Create two vectors of points such as point A latitude lies between these two latitudes. First if
    // consider if longitude in vector of points are sorted ascending or descending
    while(split)
    {
        if(y < vec.at(half).y)
        {
            if(y >= vec.at(half-1).y)
            {
                return &vec.at(half-1);
            }

            if(interval/2 == 0 && (half-2) >= 0)
            {
                --half;
            }
            else
            {
                half -= interval/2;
                interval /= 2;
            }
        }
        else
        {
            if(y < vec.at(half+1).y)
            {
                return &vec.at(half);
            }
            if(interval/2 == 0 && (half+2) <= size)
            {
                ++half;
            }
            else
            {
                half += interval/2;
                interval /= 2;
            }
        }
    }
}
/**
 * @param m - degree of spherical harmonics
 * @param n - order of spherical harmonics
 * @param value - value of Fourier coefficient
 * @param choice - LAP*cos(n*longitude) - true An coefficient, LAP*sin(n*longitude) - false Bn coefficient
 */
void model::setCoefficient(const int& m, const int& n, const double& value, bool choice)
{
    if(m*n > _width*_width)
    {
        std::cerr<<"Out of range"<<'\n';
    }
    coefficients.at(static_cast<size_t>((choice ? m * _width + n : (_width*_width) - (m * _width + n)))) = value;
}
/**
 * @param m - degree of spherical harmonics
 * @param n - order of spherical harmonics
 * @param choice - LAP*cos(n*longitude) - true, LAP*sin(n*longitude) - false
 * @return value of Fourier coefficient
 */
double model::getCoefficients(const int& m, const int& n, bool choice)const
{
    return coefficients.at(static_cast<size_t>((choice ? m * _width + n : (_width*_width) - (m * _width + n))));
}
/**
 * @brief Add value of gravitational potential to partial sum
 * @param temp - gravitational potential growth
 * @param m - degree of spherical harmonics
 * @param n - order of spherical harmonics
 */
void model::addPotential(const double& temp, const int& m, const int& n)
{
    std::lock_guard<std::mutex> guard(lock);
    int pos = (int) ((m * (m + 1)) / 2 + n);
    std::cout << "Pos: " << pos << '\n';
    gravity_potential += temp;
    potential_diff.at(pos) = temp;
}
/**
 * @brief Write values of gravitational potential in partial Fourier series to TXT file
 * @param output - output file where values of gravitational potential will be written
 */
void model::writeToTxt(std::ofstream& output)
{
    int max = getMax();
    int temp_m{};
    int temp_n{};
    int pos{};

    processDiff();

    while(true)
    {
        if(temp_m == max && temp_n == max)
        {
            output<<std::scientific<<std::setprecision(12)<<"Order: "<<temp_m<<", degree: "
                  <<temp_n<<" Coefficient An: "<<getCoefficients(temp_m, temp_n, true)<<
                  " Coefficient Bn: "<<getCoefficients(temp_m, temp_n, false)<<" Value of gravity potential: "<<
                  getPotentialDiff(pos)<<'\n';
            break;
        }
        else if(temp_m != temp_n)
        {
            if(!temp_n)
            {
                output<<std::scientific<<std::setprecision(12)<<"Order: "<<temp_m<<", degree: "
                      <<temp_n<<" Coefficient An: "<<getCoefficients(temp_m, temp_n, true)<<
                      " Coefficient Bn: "<<0.00<<" Value of gravity potential: "<<getPotentialDiff(pos)<<'\n';
            }
            else
            {
                output<<std::scientific<<std::setprecision(12)<<"Order: "<<temp_m<<", degree: "
                      <<temp_n<<" Coefficient An: "<<getCoefficients(temp_m, temp_n, true)<<
                      " Coefficient Bn: "<<getCoefficients(temp_m, temp_n, false)<<" Value of gravity potential: "<<
                      getPotentialDiff(pos)<<'\n';
            }
            ++temp_n;
        }
        else
        {
            output<<std::scientific<<std::setprecision(12)<<"Order: "<<temp_m<<", degree: "
                  <<temp_n<<" Coefficient An: "<<getCoefficients(temp_m, temp_n, true)<<
                  " Coefficient Bn: "<<(temp_m ? getCoefficients(temp_m, temp_n, false): 0.00)<<" Value of gravity potential: "<<
                  getPotentialDiff(pos)<<'\n';
            ++temp_m;
            temp_n = 0;
        }
        ++pos;
        output.flush();
    }
}
/**
 * @brief Sum up value of gravitational potential for each partial sum of Fourier series
 */
void model::processDiff()
{
    size_t temp = {1};

    while(temp != potential_diff.size())
    {
        potential_diff.at(temp) += potential_diff.at(temp - 1);
        ++temp;
    }
}
