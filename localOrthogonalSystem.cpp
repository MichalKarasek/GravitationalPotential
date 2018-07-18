//
// Created by michal on 15.3.18.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include "localOrthogonalSystem.h"

/**
 * @param m - maxim degree of local orthogonal system that will be created
 * @param bound - boundary in which will be system orthogonal
 */
localOrthogonalSystem::localOrthogonalSystem(const int m, const std::array<double, 4> &bound):
boundary(bound), _max(m)
{
    M = new matrix(m);
    computeValues(m);
    coeffToTxt();
}
/**
 * @param m - maximum degree of local orthogonal sphere function
 * @param bound - boundary in which are functions orthogonal
 * @param path - path to the source file containing local ortohogonal system
 */
localOrthogonalSystem::localOrthogonalSystem(const int m, const std::array<double, 4> &bound, const std::string &path):
boundary(bound), _max(m)
{
    double temp{};
    int counter{};
    int i{};
    std::string line;
    std::ifstream file;
    std::string::size_type sz;
    std::vector<double> coefficients;
    coefficients.reserve((u_int)(m*(m+1))/2);
    M = new matrix(m, false);

    file.open(path);

    if(!file.is_open())
    {
        std::cerr<<"Fail to open local orthogonal system file"<<'\n';
        exit(-1);
    }
    else
    {
        std::cout<<"Local orthogonal system file open successfully"<<'\n';
    }

    while(getline(file,line))
    {
        ++counter;

        if(counter == 1)
            coefficients.push_back(1);

        while ( i != counter+1)
        {
            temp = std::stod(line, &sz);
            line = line.substr(sz);

            if(i > 1)
            {
                coefficients.push_back(temp);
            }
            ++i;
        }
        i = 0;
    }
    M->coeffFromVec(coefficients);
}
/**
 * @brief compute local orthogonal system of spherical function. Calculation is done by Gran-Schmidt process
 *        where Spherical harmonics (https://en.wikipedia.org/wiki/Spherical_harmonics) are used as base system
 * @param max - maxim degree of local orthogonal system
 */
void localOrthogonalSystem::computeValues(const int max)
{
    int m = 1;
    int n{};
    int temp_n{};
    int temp_m{};
    double coefficient{};
    bool arr[2] = {true, false};
    std::vector<double> vec;

    while(m <= max && n <= max)
    {
        std::cout<<"Compute "<<m<<" "<<n<<'\n';
        auto start = std::chrono::high_resolution_clock::now();

        if(n == 0)
        {
            do
            {
                if(temp_n == 0)
                {
                    coefficient = computeCoefficients(temp_m, temp_n, m, n, true, true);
                    vec.push_back(coefficient);
                    if(temp_m == 0)
                    {
                        temp_m += 1;
                    }
                    else
                    {
                        temp_n += 1;
                    }
                }
                else
                {
                    coefficient = computeCoefficients(temp_m, temp_n, m, n, true, true);
                    vec.push_back(coefficient);
                    coefficient = computeCoefficients(temp_m, temp_n, m, n, true, false);
                    vec.push_back(coefficient);

                    if(temp_m == temp_n)
                    {
                        temp_m +=1;
                        temp_n = 0;
                    }
                    else
                    {
                        temp_n += 1;
                    }
                }
            }while(temp_m < m && temp_n <= temp_m);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time = end - start;
            duration_sys.push_back(time.count());

            M->setCoefficients(vec, m, n, true);
            vec.clear();
            vec.resize(0);
            temp_m = 0;
            temp_n = 0;
        }
        else
        {
            for (size_t i = 0; i < 2; ++i)
            {
                do
                {
                    if(temp_n == 0)
                    {
                        coefficient = computeCoefficients(temp_m, temp_n, m, n, arr[i], true);
                        vec.push_back(coefficient);

                        temp_m == 0 ? temp_m += 1 : temp_n += 1;
                    }
                    else
                    {
                        if(temp_m == m && temp_n == n)
                        {
                            if(!arr[i])
                            {
                                coefficient = computeCoefficients(temp_m, temp_n, m, n, arr[i], true);
                                vec.push_back(coefficient);
                            }
                        }
                        else
                        {
                            coefficient = computeCoefficients(temp_m, temp_n, m, n, arr[i], true);
                            vec.push_back(coefficient);
                            coefficient = computeCoefficients(temp_m, temp_n, m, n, arr[i], false);
                            vec.push_back(coefficient);
                        }

                        if(temp_m == temp_n)
                        {
                            temp_m +=1;
                            temp_n = 0;
                        }
                        else
                        {
                            temp_n += 1;
                        }
                    }
                }while((temp_m < m && temp_n <= temp_m) || (temp_m == m && temp_n <= n));

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time = end - start;
                duration_sys.push_back(time.count());

                M->setCoefficients(vec, m, n, arr[i]);
                vec.clear();
                vec.resize(0);
                temp_m = 0;
                temp_n = 0;
            }
        }
        if (m == n)
        {
            m += 1;
            n = 0;
        }
        else
        {
            n += 1;
        }
    }
}
/**
 * @param m_max - degree of local spherical harmonics
 * @param n_max - order of local spherical harmonics
 * @param x - latitude in which is function evaluated
 * @param y - longitude in which is function evaluated
 * @param choice - are we interested in spherical harmonics in shape LAP*cos(order*longitude) - true, or LAP*sin(order*longitude) - false
 * @return return the value of local spherical harmonics of specified degree and order at specified longitude and latitude
 */
double localOrthogonalSystem::getOrthoFunc(int m_max, int n_max, double x, double y, bool choice)
{
    double Sum{0};
    int temp_m{};
    int temp_n{};
    auto pair = M->getPointers(m_max, n_max, choice);
    if(m_max == 0)
    {
        return algorithms::sphereFunction(x, y, m_max, n_max, true);
    }

    std::vector<double > vec = M->getCoefficients(pair);
    std::vector<double>::iterator it = vec.begin();

    Sum += algorithms::sphereFunction(x, y, m_max, n_max, choice);

    while(it != vec.end())
    {
        if(temp_n == 0)
        {
            Sum += *it*algorithms::sphereFunction(x, y, temp_m, temp_n, true);
            ++it;

            temp_m == 0 ? temp_m += 1 : temp_n += 1;
        }
        else
        {
            if(temp_m == m_max && temp_n == n_max)
            {
                if(!choice)
                {
                    Sum += *it*algorithms::sphereFunction(x, y, temp_m, temp_n, true);
                    ++it;
                }
            }
            else
            {
                Sum += *it*algorithms::sphereFunction(x, y, temp_m, temp_n, true);
                ++it;
                Sum += *it*algorithms::sphereFunction(x, y, temp_m, temp_n, false);
                ++it;
            }

            if (temp_m == temp_n)
            {
                temp_m+=1;
                temp_n = 0;
            }
            else
            {
                temp_n += 1;
            }
        }
    }
    return Sum;
}
/**
 * @param m
 * @param n
 * @param k
 * @param l
 * @param choice
 * @param choice2
 * @return
 */
double localOrthogonalSystem::computeCoefficients(int m, int n,int k, int l, bool choice, bool choice2)
{
    double  coeff{};
    double  temp = M->getNorm(m, n, choice2);

    coeff = gaussINT(boundary, m, n, k, l, choice2, choice);

    if(temp == -1)
    {
        temp = gaussINT(boundary, m, n, choice2);
        M->setNorm(m, n, temp, choice2);
    }
    return (-1)*(coeff/temp);
}
/**
 * @param interval
 * @param m
 * @param k
 * @param n
 * @param l
 * @param choice
 * @param choice1
 * @return
 */
double localOrthogonalSystem::gaussIntegration(const std::array<double , 4> &interval, const int m, const int k, const int n, const int l, bool choice, bool choice1)
{
    double result{};
    double x{};
    double y{};

    double weights [] = {0.46791393457269104738987034398, 0.46791393457269104738987034398, 0.36076157304813860756983351384, 0.36076157304813860756983351384, 0.171324492379170345040296142182, 0.171324492379170345040296142182};
    double roots [] = {-0.23861918608319690863050172168, 0.23861918608319690863050172168, -0.66120938646626451366139959501, 0.66120938646626451366139959501, -0.93246951420315202781230155449, 0.93246951420315202781230155449};

    for (size_t i = 0; i < sizeof(weights)/sizeof(double); ++i)
    {
        for (size_t j = 0; j < sizeof(weights)/sizeof(double); ++j)
        {
            x = ((interval.at(1)-interval.front()) / 2.0 * roots[i] + (interval.at(1)+interval.front())/2.0);
            y = ((interval.back()-interval.at(2)) / 2.0 * roots[j] + (interval.at(2)+interval.back())/2.0);
            result+= weights[i]*weights[j]*getOrthoFunc(m, k, x, y, choice)*algorithms::sphereFunction(x, y, n, l, choice1)*cos(x);
        }
    }
    result *= ((interval.back()-interval.at(2))*(interval.at(1)-interval.front()))/4.0;

    return result;
}
/**
 * @param interval
 * @param m
 * @param k
 * @param choice
 * @return
 */
double localOrthogonalSystem::gaussIntegration(const std::array<double , 4> &interval, const int m, const int k, bool choice)
{
    double result{};
    double x{};
    double y{};

    double weights [] = {0.46791393457269104738987034398, 0.46791393457269104738987034398, 0.36076157304813860756983351384, 0.36076157304813860756983351384, 0.171324492379170345040296142182, 0.171324492379170345040296142182};
    double roots [] = {-0.23861918608319690863050172168, 0.23861918608319690863050172168, -0.66120938646626451366139959501, 0.66120938646626451366139959501, -0.93246951420315202781230155449, 0.93246951420315202781230155449};

    for (size_t i = 0; i < sizeof(weights)/sizeof(double); ++i)
    {
        for (size_t j = 0; j < sizeof(weights)/sizeof(double); ++j)
        {
            x = ((interval.at(1)-interval.front()) / 2.0 * roots[i] + (interval.at(1)+interval.front())/2.0);
            y = ((interval.back()-interval.at(2)) / 2.0 * roots[j] + (interval.at(2)+interval.back())/2.0);
            result+= weights[i]*weights[j]*getOrthoFunc(m, k, x, y, choice)*getOrthoFunc(m, k, x, y, choice)*cos(x);
        }
    }
    result *= ((interval.back()-interval.at(2))*(interval.at(1)-interval.front()))/4.0;

    return result;
}
/**
 * @param arr
 * @param m
 * @param k
 * @param n
 * @param l
 * @param choice
 * @param choice2
 * @return
 */
double localOrthogonalSystem::gaussINT(const std::array<double, 4> &arr, const int m, const int k, const int n, const int l, bool choice, bool choice2)
{
    std::vector<std::array<double, 4>> IN;
    std::array<double, 4> arr1;
    std::array<double, 4> arr2;
    std::array<double, 4> arr3;
    std::array<double, 4> arr4;
    double result{};
    IN.push_back(arr);
    double gauss_int2{};
    double diff{};
    double gauss_int_last{};

    double gauss_int = gaussIntegration(arr, m, k, n, l, choice, choice2);

    do
    {
        arr1 = {IN.back().front(), (IN.back().at(1) + IN.back().front()) / 2.0, IN.back().at(2), (IN.back().at(2)+IN.back().back()) / 2.0};
        arr2 = {(IN.back().at(1) + IN.back().front()) / 2.0, IN.back().at(1), IN.back().at(2), (IN.back().at(2)+IN.back().back()) / 2.0};
        arr3 = {IN.back().front(), (IN.back().at(1) + IN.back().front()) / 2.0, (IN.back().at(2)+IN.back().back()) / 2.0, IN.back().back()};
        arr4 = {(IN.back().at(1) + IN.back().front()) / 2.0, IN.back().at(1), (IN.back().at(2)+IN.back().back()) / 2.0, IN.back().back()};

        gauss_int2 = gaussIntegration(arr1, m, k, n, l, choice, choice2) + gaussIntegration(arr2, m, k, n, l, choice, choice2)
                     + gaussIntegration(arr3, m, k, n, l, choice, choice2);
        /* Compute last sub-interval separately due to the fact that last interval
         * of computation is in each cycle divide in four smaller sub-intervals and then compare
         * with result reach by the sum of all four sub-intervals. No need to compute it again */
        gauss_int_last = gaussIntegration(arr4, m, k, n, l, choice, choice2);
        gauss_int2 += gauss_int_last;
        diff = std::abs(gauss_int - gauss_int2);

        /*Set the precision of calculation. Compare value of last stored interval
         * computed with and without dividing into smaller sub-intervals */
        if(diff < 1e-12)
        {
            IN.pop_back();
            result+= gauss_int2;
            if(IN.size())
                gauss_int = gaussIntegration(IN.back(), m, k, n, l, choice, choice2);
        }
        else
        {
            IN.pop_back();
            IN.push_back(arr1);
            IN.push_back(arr2);
            IN.push_back(arr3);
            IN.push_back(arr4);
            gauss_int = gauss_int_last;
        }
    }while(!IN.empty());

    return result;
}
/**
 *
 * @param arr
 * @param m
 * @param k
 * @param choice
 * @return
 */
double localOrthogonalSystem::gaussINT(const std::array<double, 4> &arr, const int m, const int k, bool choice)
{
    std::vector<std::array<double, 4>> IN;
    std::array<double, 4> arr1;
    std::array<double, 4> arr2;
    std::array<double, 4> arr3;
    std::array<double, 4> arr4;
    double result{};
    IN.push_back(arr);
    double gauss_int2{};
    double diff{};
    double gauss_int_last{};

    double gauss_int = gaussIntegration(arr, m, k, choice);

    do
    {
        arr1 = {IN.back().front(), (IN.back().at(1) + IN.back().front()) / 2.0, IN.back().at(2), (IN.back().at(2)+IN.back().back()) / 2.0};
        arr2 = {(IN.back().at(1) + IN.back().front()) / 2.0, IN.back().at(1), IN.back().at(2), (IN.back().at(2)+IN.back().back()) / 2.0};
        arr3 = {IN.back().front(), (IN.back().at(1) + IN.back().front()) / 2.0, (IN.back().at(2)+IN.back().back()) / 2.0, IN.back().back()};
        arr4 = {(IN.back().at(1) + IN.back().front()) / 2.0, IN.back().at(1), (IN.back().at(2)+IN.back().back()) / 2.0, IN.back().back()};

        gauss_int2 = gaussIntegration(arr1, m, k, choice) + gaussIntegration(arr2, m, k, choice)
                     + gaussIntegration(arr3, m, k, choice);
        /* Compute last sub-interval separately due to the fact that last interval
         * of computation is in each cycle divide in four smaller sub-intervals and then compare
         * with result reach by the sum of all four sub-intervals. No need to compute it again */
        gauss_int_last = gaussIntegration(arr4, m, k, choice);
        gauss_int2 += gauss_int_last;
        diff = std::abs(gauss_int - gauss_int2);

        /*Set the precision of calculation. Compare value of last stored interval
         * computed with and without dividing into smaller sub-intervals */
        if(diff < 1e-12)
        {
            IN.pop_back();
            result+= gauss_int2;
            if(IN.size())
                gauss_int = gaussIntegration(IN.back(), m, k, choice);
        }
        else
        {
            IN.pop_back();
            IN.push_back(arr1);
            IN.push_back(arr2);
            IN.push_back(arr3);
            IN.push_back(arr4);
            gauss_int = gauss_int_last;
        }
    }while(!IN.empty());

    return result;
}
/**
 *
 */
void localOrthogonalSystem::coeffToTxt()
{
    //TODO upravit tvorbu systemu aby v coefficients byla na zacatku jednicka pro P00 ted musim delat navic jeden if
    std::ofstream output;
    output.open( "/home/michal/Desktop/Koefcienty.txt", std::ofstream::out | std::ofstream::app);

    int temp_m{};
    int temp_n{};
    std::pair<double*, double*> pointers;
    std::vector<double> coefficients;

    while (temp_m <= _max && temp_n <= _max)
    {
        if (temp_m == 0 || temp_n == 0)
        {
            if(temp_m == 0 && temp_n == 0)
            {
                coefficients = {1};
                writeCoeff(output, coefficients, temp_m, temp_n);
            }
            else
            {
                pointers = M->getPointers(temp_m, temp_n, true);
                coefficients = M->getCoefficients(pointers);
                writeCoeff(output, coefficients, temp_m, temp_n);
            }
            temp_m ? ++temp_n : ++temp_m;
        }
        else
        {
            pointers = M->getPointers(temp_m, temp_n, true);
            coefficients = M->getCoefficients(pointers);
            writeCoeff(output, coefficients, temp_m, temp_n);

            pointers = M->getPointers(temp_m, temp_n, false);
            coefficients = M->getCoefficients(pointers);
            writeCoeff(output, coefficients, temp_m, temp_n);

            if(temp_m == temp_n)
            {
                ++temp_m;
                temp_n = 0;
            }
            else
            {
                ++temp_n;
            }
        }
    }
}
/**
 * @param output
 * @param coeff
 * @param m
 * @param n
 */
void localOrthogonalSystem::writeCoeff(std::ofstream& output, const std::vector<double>& coeff, const double m, const double n)
{
    output<<std::setprecision(3)<<m<<" "<<n;

    for(auto& coefficient:coeff)
    {
        output<<std::scientific<<std::setprecision(15)<<" "<<coefficient;
    }
    //TODO upravit coefficients pro 00
    if(m)
    {
        output<<"Time to compute: "<<duration_sys.at(0)<<'\n';
        duration_sys.erase(duration_sys.begin());
        output.flush();
    }
}
/**
 * @param m
 * @param n
 * @param choice
 * @return
 */
std::vector<double> localOrthogonalSystem::getVector(const int m, const int n, bool choice)
{
    std::pair<double*, double*> pair;

    pair = M->getPointers(m,n,choice);
    return M->getCoefficients(pair);
}

//TEST FUNCTION
double localOrthogonalSystem::getOrthoFunc(const int m_max, const int n_max, const Point &p, bool choice)
{
    double Sum{0};
    double temp_m = 0;
    double temp_n = 0;
    auto pair = M->getPointers(m_max, n_max, choice);
    if(m_max == 0)
    {
        return algorithms::sphereFunction(p.x, p.y, m_max, n_max, true);
    }

    std::vector<double > vec = M->getCoefficients(pair);
    std::vector<double>::iterator it = vec.begin();

    Sum += algorithms::sphereFunction(p.x, p.y, m_max, n_max, choice);

    while(it != vec.end())
    {
        if(temp_n == 0)
        {
            Sum += *it*algorithms::sphereFunction(p.x, p.y, (int) temp_m, (int) temp_n, true);
            ++it;

            temp_m == 0 ? temp_m += 1 : temp_n += 1;
        }
        else
        {
            if(temp_m == m_max && temp_n == n_max)
            {
                if(!choice)
                {
                    Sum += *it*algorithms::sphereFunction(p.x, p.y, (int) temp_m, (int) temp_n, true);
                    ++it;
                }
            }
            else
            {
                Sum += *it*algorithms::sphereFunction(p.x, p.y, (int) temp_m, (int) temp_n, true);
                ++it;
                Sum += *it*algorithms::sphereFunction(p.x, p.y, (int) temp_m, (int) temp_n, false);
                ++it;
            }

            if (temp_m == temp_n)
            {
                temp_m+=1;
                temp_n = 0;
            }
            else
            {
                temp_n += 1;
            }
        }
    }
    return Sum;
}
