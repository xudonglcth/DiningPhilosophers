#include <iostream>
#include <fstream>
#include "System.h"
#include "DiningPhilosophers.h"
#include "DiningPhilosophersNew.h"
#include "readVLTS.h"
#include "SyncInT.h"
#include "DiningPhilosophersEff.h"
#include "test.h"
#include "GroupDining.h"
int main() {
/*
    //Sync test
    {
        //System G1({0, 1}, {{0, 2, 1},
        //                   {1, 3, 1}}, {1, 2}, {0});
        //System G1_({0, 1}, {{0, 3, 1},
        //                    {1, 4, 1}}, {1, 2}, {0});
        //G1_.sync(G1);
        //G1.syncInT(G1_);

        System G4({0, 1, 2}, {{0, 2, 1},
                               {0, 4, 1},
                               {0, 3, 2}}, {1, 1, 2}, {0});
        System G4_({0, 1}, {{0, 2, 1},
                              {0, 3, 1}}, {1, 2}, {0});
        G4.init = {0};
        G4_.init = {0};
        //G4.sync(G4_);
        G4_.syncInT(G4);
        for (const auto & iter : G4_.transitions){
            std::cout << iter[0] << " " << iter[1] << " " << iter[2] << std::endl;
        }

        std::cout << std::endl;
        System G5({0, 10}, {{0, 2, 10},
                            {10, 3, 0}}, {1, 2}, {0});
        System G5_({0, 10}, {{0, 3, 10},
                             {10, 4, 0}}, {1, 2}, {0});
        G5.init = {0};
        G5_.init = {0};
        G5_.syncInT(G5);
        for (const auto & iter : G5_.transitions){
            std::cout << iter[0] << " " << iter[1] << " " << iter[2] << std::endl;
        }

        System G0 ({0, 1, 2, 3, 4}, {{0, 0, 1}, {1, 1, 2}, {2, 2, 3}, {3, 3, 4}}, {0});
        System G0_({0, 1, 2, 3, 4}, {{0, 4, 1}, {1, 1, 2}, {2, 4, 3}, {3, 2, 4}}, {0});
        G0.syncInT(G0_);

        System G1 ({0, 1}, {{0, 2, 1}, {1, 3, 1}}, {0});
        System G1_ ({0, 1}, {{0, 3, 1}, {1, 4, 1}}, {0});
        G1.syncInT(G1_);

        System G6({0, 1}, {{0, 2, 0}, {0, 3, 1}}, {0});
        System G6_({0, 1, 2}, {{0, 2, 1}, {1, 3, 2}}, {0});
        G6.syncInT(G6_);

        System G7({0, 1, 2}, {{0, 5, 1}, {1, 6, 2}, {2, 3, 2}}, {0});
        System G7_({0, 1, 2}, {{0, 6, 1}, {1, 8, 2}, {2, 3, 2}}, {0});
        G7.syncInT(G7_);

        System G8({0, 1}, {{0, 6, 1}, {1, 3, 1}}, {0});
        System G8_({0, 1, 2}, {{0, 6, 1}, {1, 8, 2}}, {0});
        G8.syncInT(G8_);

        System G10({0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                   {{0, 2, 1}, {1, 3, 2}, {2, 4, 3}, {3, 5, 4}, {4, 6, 5}, {5, 7, 6}, {6, 8, 7}, {7, 9, 8}, {8, 10, 9},
                    {0, 12, 1}, {1, 13, 3}}, {0});
        System G10_({0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                    {{0, 2, 1}, {1, 3, 2}, {2, 4, 3}, {3, 5, 4}, {4, 6, 5}, {5, 7, 6}, {6, 8, 7}, {7, 9, 8},
                     {8, 10, 9}}, {0});
        G10.syncInT(G10_);

        System G11({0, 1, 2}, {{0, 6, 1}, {1, 10 ,2}}, {0});
        System G11_({0, 1, 2}, {{0, 7, 2}, {2, 6, 1}}, {0});
        G11.syncInT(G11_);

    }
*/
    //syncThenAbsEff(8);
    //syncThenAbsNew(10);
    clock_t t1 = clock();
    groupIn3(21);
    clock_t t2 = clock();
    std::cout<<"Incremental Approach Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
    //syncAndAbsNew(14);
    /*
    clock_t t1 = clock();
    syncAndAbstract(10);
    clock_t t2 = clock();
    std::cout<<"Incremental Approach Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;

    t1 = clock();
    syncThenAbstract(10);
    t2 = clock();
    std::cout<<"Brute-Force Approach Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
     */
    ///testPerformance("vasy_40_60");
    //syncAndAbsNew(14);
    //syncThenAbsNew(7);
    return 0;
}
