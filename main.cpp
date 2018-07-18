#include <iostream>
#include <vector>
#include <fstream>
#include <thread>
#include "gravAlgorithms.h"

// TODO - VSTUP Z PRIKAZOVY RADKY NENI KOMPLET ZATIM NEFUNGUJE !!!!
// TODO - GRAFIKA V QT

#if __cplusplus < 201103L
    #warning "USE C++11 STANDARD TO RUN THE PROGRAM"
    #error
#endif

/*
//std::string file_path = "/home/michal/Desktop/Mymatrix903602018.txt";
std::string file_path = "/home/michal/Desktop/LokalniPotencial.txt";
std::string filepath = "/home/michal/Desktop/Lokal10.txt";
//std::string file_path = "/home/michal/Desktop/MymatrixSmall.txt";
//std::string file_path = "/home/michal/Desktop/91359.txt";
//std::string file_path = "/home/michal/Desktop/Zdroje_diplomka/Hodnoty7.txt";
model Gravitational_model(file_path, 2405, 10);
std::array<double , 4> boundary ={PI/8, PI/6, PI/10, PI/2};
//std::array<double , 4> boundary ={-PI/2, PI/2, 0, 2*PI};
localOrthogonalSystem System(10, boundary, filepath);
//model Gravitational_model(file_path, 91, 359, 1);
*/

int main(int argv, char** argc)
{
    std::string file_path;
    std::string output_path;
    std::ofstream output;
    Point p;
    const std::array<double, 4>& arr={-PI/2, PI/2, 0, 2*PI};

    int alg_choice {}; // What of spatial algorithm interpolation to use ?
    bool global;

    if(argv < 2)
    {
        std::cerr<<"No input file provided !!!"<<"\n";
        exit(EXIT_FAILURE);
    }
    else
    {
        if (argv == 8 && argc[1] == "-g")
        {
            global = true;
            file_path = argc[2];
            output_path = argc[4];
            p.potencial = 0;
            p.x = std::stod((std::string)argc[6]);
            p.y = std::stod((std::string)argc[7]);
        }
    }

    model Gravitational_model(file_path, 2405, 10);
    //localOrthogonalSystem System(10, );

    output.open(output_path, std::ofstream::out | std::ofstream::app);

    if(!output.is_open())
    {
        std::cerr<<"Fail to open output data file !!!"<<'\n';
        exit(EXIT_FAILURE);
    }

    auto start = std::chrono::high_resolution_clock::now();

    // return number of cores available (since c++11)
    u_int core_num = std::thread::hardware_concurrency();

    // create vector of threads
    std::vector<std::thread> threads;

    // Based on the area of computaion need to specify gravitationalPotential bool - false - local, true - Gloabal
    if(core_num > 7)
    {
        for (size_t i = 0; i < 4; ++i)
            threads.emplace_back(std::thread(gravAlgorithms::gravitationalPotential, std::ref(Gravitational_model), std::ref(System), p, global, alg_choice, arr));
    }
    else if (core_num <= 4)
    {
        for (size_t i = 0; i < 2; ++i)
            threads.emplace_back(std::thread(gravAlgorithms::gravitationalPotential, std::ref(Gravitational_model), std::ref(System), p, global, alg_choice, arr));
    }

    for(auto &thread:threads)
    {
        thread.join();
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;
    std::cout<<"Ellapsed time: "<<time.count()<<'\n';
    Gravitational_model.writeToTxt(output);

    output.close();

    return 0;
}
