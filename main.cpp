#include <iostream>
#include "System.h"
#include "DiningPhilosophers.h"
#include "DiningPhilosophersNew.h"
#include "readVLTS.h"
int main() {
    //Abstraction test
    //test G16
    /*
    std::cout << "Test G16" << std::endl;
    System::size_t_in_vec state_g16 = {0, 1, 2, 3, 4};
    System::size_t_in_vec_in_vec trans_g16 = {{0, 2, 1},
                            {1, 1, 2},
                            {2, 1, 2},
                            {0, 1, 3},
                            {3, 1, 4},
                            {4, 1, 4}};
    std::vector<size_t > partitions_new_g16 = {2, 2, 2, 2, 2};
    System G16(state_g16, trans_g16, partitions_new_g16, {1});
    G16.transReduce();

    std::cout << "Test G17" << std::endl;
    System::size_t_in_vec state_g17 = {0, 1, 2};
    System::size_t_in_vec_in_vec trans_g17 = {{0, 1, 1},
                            {0, 1, 2},
                            {1, 2, 1},
                            {2, 3, 2}};
    std::vector<size_t> partitions_new_17 = {};
    System G17(state_g17, trans_g17, partitions_new_17, {1});
    G17.transReduce();

    //test G19
    std::cout << "Test G19" << std::endl;
    System::size_t_in_vec state_g19 = {0, 1, 2, 3};
    System::size_t_in_vec_in_vec trans_g19 = {{0, 1, 1},
                            {1, 2, 2},
                            {0, 2, 2},
                            {0, 1, 3},
                            {2, 1, 3},
                            {3, 2, 3}};
    std::vector<size_t> partitions_new_19 = {2, 2, 2, 1};
    System G19(state_g19, trans_g19, partitions_new_19, {1});
    G19.transReduce();
    std::cout << "Reduced Transitions #: " << G19.transitions.size()
    << "\nReduced States #: " << *std::max_element(G19.partitions_new.begin(), G19.partitions_new.end()) + 1<<std::endl;

    //Sync test
    System G1({0, 1}, {{0, 2, 1}, {1, 3, 1}}, {1, 2}, {0});
    System G1_({0, 1}, {{0, 3, 1}, {1, 4, 1}}, {1, 2}, {0});
    G1_.sync(G1);

    System G4({0, 1, 2}, {{0, 2, 1}, {0, 4, 1}, {0, 3, 2}}, {1, 1, 2}, {0});
    System G4_({0, 1}, {{0, 2, 1}, {0, 3, 1}}, {1, 2}, {0});
    G4.sync(G4_);

    System G5({0, 1}, {{0, 2, 1}, {1, 3, 0}}, {1, 2}, {0});
    System G5_({0, 1}, {{0, 3, 1}, {1, 4, 0}}, {1, 2}, {0});
    G5.sync(G5_);
*/
    /*
    clock_t t1 = clock();
    syncAndAbstract(10);
    clock_t t2 = clock();
    std::cout<<"Incremental Approach Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;

    t1 = clock();
    syncThenAbstract(10);
    t2 = clock();
    std::cout<<"Brute-Force Approach Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
    //testPerformance("vasy_18_73");
     */

    //syncAndAbsNew(14);
    syncThenAbsNew(8);
    return 0;
}
