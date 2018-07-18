//
// Created by michal on 17.7.18.
//
#include <thread>
#include "gravAlgorithms.h"

namespace safety
{
    std::atomic<int> m{};
    std::atomic<int> n{};
    std::mutex stop;
}

/**
* @brief Only callable function in class. It is basic function from which all other functions need for calculation of
*        gravitational potential on Earth surface are being called.
* @param Model - gravitational potential model
* @param sys - local orthogonal system
* @param p - Point in which is the gravitational model compute
* @param global_model - true if global gravitational model is available, false - if local gravitational model
* @param alg_choice - interpolation method algorithm
* @param arr - interval on which is model defined. Needs to be set in case of lacal gravitational model
*/
void gravAlgorithms::gravitationalPotential(model& model, localOrthogonalSystem& sys, const Point& p, bool global_model, int alg_choice,
                                            const std::array<double, 4>& arr={-PI/2, PI/2, 0, 2*PI})
{
    int temp_m{};
    int temp_n{};

    while(true)
    {
        if(safety::m == model.getMax() && safety::n == model.getMax())
        {
            temp_m = safety::m.load();
            temp_n = safety::n.load();
            if(!model.getCoefficients(static_cast<int>(safety::m), static_cast<int>(safety::n), false))
            {
                model.setCoefficient(static_cast<int>(safety::m), static_cast<int>(safety::n), 1, false);
                compute(temp_m, temp_n, model, sys, p, global_model, alg_choice, arr);
            }
            break;
        }
        else if(safety::m.load() != safety::n.load())
        {
            safety::n.fetch_add(1);
            temp_m = safety::m.load();
            temp_n = safety::n.load()-1;

        }
        else
        {
            safety::m.fetch_add(1);
            temp_m = safety::m.load()-1;
            temp_n = safety::n.load();
            safety::n.store(0);
        }
        compute(temp_m, temp_n, model, sys, p, global_model, alg_choice, arr);
    }
}
/**
 * @param m - degree of Legendre associated polynomials
 * @param k - order of Legendre associated polynomials
 * @param moz - parameter for Function and FunctionLocal
 * @param Model - gravitational potential model
 * @param choice - parameter for Function and FunctionLocal
 * @param global_model - true if global gravitational model is available, false - if local gravitational model
 * @param alg_choice - interpolation method algorithm
 * @param arr - interval on which is model defined. Needs to be set in case of lacal gravitational model
 * @param idwpar - power parameter, by increasing its value the nearest neighbour multiplied impact
 *                 on final gravitational potential value
 * @return integration value
 */
double gravAlgorithms::gaussINT(const int m, const int k, const int moz, model& Model, localOrthogonalSystem& system, bool choice, bool global_model,
                            int alg_choice, const std::array<double, 4>& arr, const int idwpar)
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
    double gauss_int = gauss(arr, Model, m, k, moz, choice, global_model, alg_choice, idwpar, system);

    do
    {
        arr1 = {IN.back().front(), (IN.back().at(1) + IN.back().front()) / 2.0, IN.back().at(2), (IN.back().at(2)+IN.back().back()) / 2.0};
        arr2 = {(IN.back().at(1) + IN.back().front()) / 2.0, IN.back().at(1), IN.back().at(2), (IN.back().at(2)+IN.back().back()) / 2.0};
        arr3 = {IN.back().front(), (IN.back().at(1) + IN.back().front()) / 2.0, (IN.back().at(2)+IN.back().back()) / 2.0, IN.back().back()};
        arr4 = {(IN.back().at(1) + IN.back().front()) / 2.0, IN.back().at(1), (IN.back().at(2)+IN.back().back()) / 2.0, IN.back().back()};

        gauss_int2 = gauss(arr1, Model, m, k, moz, choice, global_model, alg_choice, idwpar, system) + gauss(arr2, Model, m, k, moz, choice, global_model, alg_choice, idwpar, system)
                     + gauss(arr3, Model, m, k, moz, choice, global_model, alg_choice, idwpar, system);
        /* Compute last sub-interval separately due to the fact that last interval
         * of computation is in each cycle divide in four smaller sub-intervals and then compare
         * with result reach by the sum of all four sub-intervals. No need to compute it again */
        gauss_int_last = gauss(arr4, Model, m, k, moz, choice, global_model, alg_choice, idwpar, system);
        gauss_int2 += gauss_int_last;
        diff = std::abs(gauss_int - gauss_int2);

        /*Set the precision of calculation. Compare value of last stored interval
         * computed with and without dividing into smaller sub-intervals */
        if(diff < 1e-4)
        {
            IN.pop_back();
            result+= gauss_int2;
            if(IN.size())
                gauss_int = gauss(IN.back(), Model, m, k, moz, choice, global_model, alg_choice, idwpar, system);
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
 * @param interval - limits of interval in which is integration evaluated
 * @param Model - Model of gravitational potential
 * @param m - degree of Legendre associated polynomials
 * @param k - order of Legendre associated polynomials
 * @param moz - parameter for Function and FunctionLocal functions
 * @return value of integration on specified interval
 */
double gravAlgorithms::gauss(const std::array<double, 4> &interval, model& Model, int m, int k, int moz, bool choice,
                         bool global_model, const int alg_choice, const int idw_par, localOrthogonalSystem& system)
{
    Point gauss_point;
    double result{};

    //double weights [] = {0.46791393457269104738987034398, 0.46791393457269104738987034398, 0.36076157304813860756983351384, 0.36076157304813860756983351384, 0.171324492379170345040296142182, 0.171324492379170345040296142182};
    //double roots [] = {-0.23861918608319690863050172168, 0.23861918608319690863050172168, -0.66120938646626451366139959501, 0.66120938646626451366139959501, -0.93246951420315202781230155449, 0.93246951420315202781230155449};
    double weights [] = {0.2729250867779006, 0.2628045445102467, 0.2628045445102467, 0.2331937645919905, 0.2331937645919905, 0.1862902109277343, 0.1862902109277343, 0.1255803694649046, 0.1255803694649046, 0.0556685671161737, 0.0556685671161737};
    double roots [] = {0, -0.2695431559523450, 0.2695431559523450, -0.5190961292068118, 0.5190961292068118, -0.7301520055740494, 0.7301520055740494, -0.8870625997680953, 0.8870625997680953, -0.9782286581460570, 0.9782286581460570};
    //double weights [] = {0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538};
    //double roots [] = {-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526};
    //double weights [] = {0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891};
    //double roots [] = {0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640};

    for (size_t i = 0; i < sizeof(weights)/sizeof(double); ++i)
    {
        for (size_t j = 0; j < sizeof(weights)/sizeof(double); ++j)
        {
            gauss_point.x = ((interval.at(1)-interval.front()) / 2.0 * roots[i] + (interval.at(1)+interval.front())/2.0);
            gauss_point.y = ((interval.back()-interval.at(2)) / 2.0 * roots[j] + (interval.at(2)+interval.back())/2.0);

            switch (alg_choice)
            {
                case 1:
                {
                    algorithms::nearestNeighbour(gauss_point, Model);
                    break;
                }
                case 2:
                {
                    algorithms::bilinearInterpolation(gauss_point, Model);
                    break;
                }
                default:
                {
                    algorithms::IDW(gauss_point, Model, idw_par);
                }
            }

            if (global_model)
            {
                result+= weights[i]*weights[j]*function(gauss_point, moz, m, k, choice);
            }
            else
            {
                result+= weights[i]*weights[j]*functionLocal(gauss_point, moz, m, k, choice, system);
            }
        }
    }
    result *= ((interval.back()-interval.at(2))*(interval.at(1)-interval.front()))/4.0;

    return result;
}
/**
 * @brief Compute value of Fourier coefficients An, Bn at point p for global gravitational potential model defined
 *        for the whole Globe
 * @param p - point of Gauss-Legendre numerical integration
 * @param moz - 1 - An, 2 - Bn, 3 - Fourier coefficient denominator
 * @param m - degree of Legendre associated polynomials
 * @param k - order of Legendre associated polynomials
 * @param choice - true - An coefficient denominator, false - Bn coefficient denominator
 * @return growth of integration in Point p
 */
double gravAlgorithms::function(const Point &p, const int moz, const int m, const int k, bool choice)
{
    double Int = 0;

    switch (moz)
    {
        case 1:
        {
            Int = p.potencial * algorithms::legenderNorm(p.x, m, k) * cos(k * p.y) * cos(p.x);
            break;
        }
        case 2:
        {
            Int = p.potencial * algorithms::legenderNorm(p.x, m, k) * sin(k * p.y) * cos(p.x);
            break;
        }
        case 3:
        {
            Int = pow((choice ? cos : sin)(k * p.y), 2.0);
            Int *= algorithms::legenderNorm(p.x, m, k) * algorithms::legenderNorm(p.x, m, k) * cos(p.x);
            break;
        }
    }
    return Int;
}
/**
 * @param temp_m - degree of Fourier coefficient
 * @param temp_n - order of Fourier coefficient
 * @param M1 - model of gravitational potential
 * @param p - point in which we calculate gravitational potential
 * @param choice - true - global model, false - local model
 */
void gravAlgorithms::compute(const int temp_m, const int temp_n, model& M1, localOrthogonalSystem& system, const Point& p, bool choice, int alg_choice, const std::array<double, 4>& boundary)
{
    double An{};
    double Bn{};
    double NormaAn{};
    double NormaBn{};
    double temp_V{};

    std::thread::id this_id = std::this_thread::get_id();
    std::cout<<"Thread number: "<< this_id << " values of m: "<<temp_m<<", n: "<<temp_n<<'\n';
    //todo an bn true true neni potreba upravit ty boundary takhle jsou globalni a to je na H...O
    An = gaussINT(temp_m, temp_n, 1, M1, system, true, choice, alg_choice, boundary);
    NormaAn = gaussINT(temp_m, temp_n, 3, M1, system, true, choice, alg_choice, boundary);
    An = (1/NormaAn)*An;

    if(temp_n)
    {
        Bn = gaussINT(temp_m, temp_n, 2, M1, system, true, choice, alg_choice, boundary);
        NormaBn = gaussINT(temp_m, temp_n, 3, M1, system, false, choice, alg_choice, boundary);
        Bn = (1/NormaBn)*Bn;
    }

    if(choice)
    {
        temp_V = (An*cos(temp_n*p.y)+Bn*sin(temp_n*p.y))*algorithms::legenderNorm(p.x,temp_m,temp_n);
    }
    else
    {
        if(temp_m == 0 && temp_n == 0)
        {
            temp_V = An;
        }
        else
        {
            temp_V = An*system.getOrthoFunc(temp_m, temp_n, p) + Bn*system.getOrthoFunc(temp_m, temp_n, p, false);
        }
    }

    std::lock_guard<std::mutex> guard(safety::stop);
    M1.addPotential(temp_V, temp_m, temp_n);
    M1.setCoefficient(temp_m, temp_n, An, true);
    if(temp_n != 0)
        M1.setCoefficient(temp_m, temp_n, Bn, false);

}
/**
 * @brief Compute value of Fourier coefficients An, Bn at point p for local gravitational potential model
 * @param p - point of Gauss-Legendre numerical integration
 * @param moz - 1 - An, 2 - Bn, 3 - Fourier coefficient denominator
 * @param m - degree of Legender associated function
 * @param k - order of Legendre associated function
 * @param choice - true - An coefficient denominator, false - Bn coefficient denominator
 * @return growth of integration in Point p
 */
double gravAlgorithms::functionLocal(const Point &p, const int moz, const int m, const int k, bool choice, localOrthogonalSystem& system)
{
    double Int = 0;

    switch (moz)
    {
        case 1:
        {
            Int = p.potencial * cos(p.x) * system.getOrthoFunc(m, k, p, true);
            break;
        }
        case 2:
        {
            Int = p.potencial * cos(p.x) * system.getOrthoFunc(m, k, p, false);
            break;
        }
        case 3:
        {
            Int = system.getOrthoFunc(m, k, p, choice) * system.getOrthoFunc(m, k, p, choice) * cos(p.x);
            break;
        }
    }
    return Int;
}
