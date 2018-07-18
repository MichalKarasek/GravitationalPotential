//
// Created by michal on 15.3.18.
//
#include <iostream>
#include <vector>
#include <algorithm>
#include "matrix.h"

matrix::matrix(const int m, bool state): _state(state)
{
    int coeff_size = ((m+1)*(1+(2*m+1)))/2;
    int coeff_num = ((coeff_size-1)*coeff_size)/2;
    _width = m+1;
    _size = _width*_width;
    Coefficients = std::vector<std::vector<double>>(static_cast<size_t>(coeff_size), std::vector<double>());
    coefficients = std::vector<double>(static_cast<size_t>(coeff_num));
    norms = std::vector<double>(_size, -1);
    _iter = &coefficients.at(0);

    for(size_t i = 0; i < _size; ++i)
    {
        Matrix.push_back(std::pair<double *, double *>());
    }
}
/**
 * @brief set pair of iterators on its right position in private data structure Matrix based
 * on the degree m and order n of OSF function and its shape defined by choice. These pairs are
 * accessed to obtain a vector of coefficients for each function creating OSF on defined interval
 * of longitude and latitude on sphere. See( Associated Legendre polynomials, Spherical harmonics)
 * for info on computation that is being made.
 *
 * @param m degree of function
 * @param n order of function
 * @param pair pair of start-end pointers pointing on coefficients
 * @param choice shape of function. true - LAF * cos(m * longitude), false - LAF*sin(m*latitude)
 */
void matrix::setPointers(const int m, const int n, const std::pair<double*, double*> pair, bool choice)
{
    if (m * n > _size)
    {
        std::cerr << "Setting matrix value " << m << ", " << n << "out of scope" << '\n';
    }

    Matrix.at(static_cast<size_t>((choice ? m * _width + n : _size - (m * _width + n)))) = pair;
}
/**
 * @brief return vector of coefficients based on start and end iterator
 *
 * @param pair of iterators on coefficients vector containing all coefficients.
 * @return vector of coefficients obtain from private coefficient vector
 */
std::vector<double> matrix::getCoefficients(const std::pair<double*, double*>& pair)
{
    std::vector<double> vec;
    double* temp = &(*pair.first);

    if(!coefficients.empty())
    {
        if(temp == pair.second)
        {
            vec.push_back(*temp);
        }
        else
        {
            while (temp != pair.second)
            {
                vec.push_back(*temp);
                ++temp;
            }
            vec.push_back(*temp);
        }
    }
    return vec;
}
 /**
  * @brief - return pair of pointers pointing on start and end coefficients of OSF function
  * @param m degree of function
  * @param n order of function
  * @param choice depend on OSF function true - Legendre associated function multiplied by cos(m*longitude), false - LAF multiplied by sin(m*longitude)
  * @return a pair of pointers to vector of coefficients
  */
std::pair<double*, double*> matrix::getPointers(int m, int n, bool choice)
{
    return Matrix.at(choice ? static_cast<size_t>(m * _width + n) : static_cast<size_t>(_size - (m * _width + n)));
}
/**
 * @brief - set vector of function coefficients to a class vector that contain all the
 * coefficients. First and last coefficients of each vector are always pointed by two iterators
 * which are stored in Matrix private vector.
 *
 * @parameters - vec - vector of coefficients. m - degree of OSF (based on Legendre associated
 * functions) n - order of OSF. choice - based on type of spheric function which coefficients
 * are set. true - spheric function of first type Legendre associated polynomial (m,n) * cos(m*longitude)
 * false - Legendre associated polynomial (m,n) * sin(m*longitude)
 */
void matrix::setCoefficients(std::vector<double>& vec, const int m, const int n, bool choice)
{
    // recalculate each coefficient so it can be used in gauss integration
    if(_state)
        vec = validateCoefficients(vec);

    int size = static_cast<int>(vec.size());

    // create a pair of iterators to store the position of function coefficients in vector of coefficients
    // of all functions
    std::pair<double*, double*> pair;
    pair.first = const_cast<double*>(_iter);
    double* temp = const_cast<double*>(_iter);
    _iter += size;
    pair.second = const_cast<double*>(_iter-1);

    for (auto& value:vec)
    {
        *temp = value;
        ++temp;
    }

    setPointers(m, n, pair, choice);
}
/**
 * @brief - coefficients for each function of OSF(orthogonal system functions) are computed
 * for the previous functions (see Schmidt orthogonalization). Each of this function can be divided
 * into spheric functions (m, n) multiplied by computed coefficients. Function guarantee that each
 * coefficient stored in coefficient private vector can be later multiplied directly by spheric function
 * without any additional recursion needed.
 *
 * @parameters - vec - vector of coefficients which real value need to be calculated and returned
 *  for storing in coefficient private vector.
 */
std::vector<double> matrix::validateCoefficients(std::vector<double> &vec)
{
    std::vector<double> temp;
    std::vector<double> valid;
    std::vector<double>::iterator temp_it;
    double  coeff{};

    for (int i = 0; i < vec.size(); ++i)
    {
        temp = getVector(i);
        temp_it = temp.begin();
        for (int j = i; j < vec.size(); ++j)
        {
            if(j == i)
            {
                coeff += vec.at(static_cast<size_t>(j));
            }
            else
            {
                coeff += *temp_it * vec.at(static_cast<size_t>(j));
                ++temp_it;
            }
        }
        addCoefficient(i, coeff);
        valid.push_back(coeff);
        coeff = 0;
    }
    return  valid;
}
/**
 * @param m position on which the vector of coefficients will be returned
 *
 * @return vector of coefficients
 */
const std::vector<double>& matrix::getVector(const int m)
{
    return Coefficients.at(static_cast<size_t>(m));
}
/**
 * @brief adding coefficient in private data structure Coefficients. Coefficients is used for
 * convenience and quick access of coefficients at desired position. Used in validateCoefficients
 *
 * @param m position on which the new coefficient x will be added
 * @param x value of added coefficients
 */
void matrix::addCoefficient(const int m, const double x)
{
    Coefficients.at(static_cast<size_t>(m)).push_back(x);
}
/**
 * @brief set norm in norms matrix to avoid multiple calculation
 * @param m degree of function
 * @param n order of function
 * @param x value to be set
 * @param choice
 */
void matrix::setNorm(const int m, const int n, const double x, bool choice)
{
    norms.at(choice ? static_cast<size_t>(m * _width + n) : static_cast<size_t>(_size - (m * _width + n))) = x;
}
/**
 * @brief Check whether the norm of spheric function has been
 * @param m degree of function
 * @param n order of function
 * @return norm at specified position
 */
double matrix::getNorm(const int m, const int n, bool choice)
{
    return norms.at(choice ? static_cast<size_t>(m * _width + n) : static_cast<size_t>(_size - (m * _width + n)));
}
/**
 * @brief
 * @param vec
 */
void matrix::coeffFromVec(const std::vector<double> &vec)
{
    coefficients = vec;
    _iter = &coefficients.at(0);

    std::pair<double*, double*> pointers;

    u_int n{};
    u_int m{};
    u_int coeff{};

    while(m != _width && n != _width)
    {
        if(m == 0)
        {
            pointers.first = (double*)_iter;
            pointers.second = (double*)_iter;
            setPointers(m, n, pointers, true);
            ++m;
        }
        else if( m != n)
        {
            if(n)
            {
                makePointers(++coeff, m, n, true);
                makePointers(++coeff, m, n, false);
            }
            else
            {
                makePointers(++coeff, m, n, true);
            }
            ++n;
        }
        else
        {
            makePointers(++coeff, m, n, true);
            makePointers(++coeff, m, n, false);

            ++m;
            n = 0;
        }
    }
}
/**
 * @brief Set pointers on coefficients
 * @param coeff -
 * @param m -
 * @param n -
 * @param choice -
 */
void matrix::makePointers(int coeff, const int m, const int n, bool choice)
{
    std::pair<double*, double*> pointers;
    pointers.first = (double*)(_iter+1);
    _iter += coeff;
    pointers.second = (double*)_iter;
    setPointers(m, n, pointers, choice);
}

