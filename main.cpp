#include <iostream>
#include <fstream>
#include "System.h"
#include "ReachBoolTest.h"
#include "DiningPhilosophers.h"
#include "DiningPhilosophersNew.h"
#include "readVLTS.h"
#include "SyncInT.h"
#include "DiningPhilosophersEff.h"
#include "test.h"
#include "GroupDining.h"
#include "IncSyncBuffer.h"
#include "twoPhils.h"
#include "diningTest2.h"
#include "GVAlgorithmClass.h"
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
/*
    {

        System dsv({0, 1, 2, 3, 4, 5}, {{0, 0, 1}, {1, 0, 0}, {1, 0, 2}, {2, 0, 3}, {3, 0, 0}, {3, 0, 3}, {2, 0, 4},
                                        {4, 0, 5}, {5, 1, 5}}, {0});
        dsv.transReduce();
        for (const auto& i : dsv.transitions){
            std::cout << i[0] <<" "<<i[1] <<" "<<i[2] << std::endl;
        }
        std::cout << std::endl;

        System dsv2({0, 1, 2, 3, 4, 5, 6, 7}, {{0, 0, 1}, {1, 0, 0}, {0, 0, 2}, {2, 0, 3},
                                                    {3, 0, 1}, {2, 0, 2}, {1, 0, 4}, {4, 0, 5},
                                                    {5, 0, 5}, {5, 0, 6}, {6, 0, 6}, {6, 0, 7}, {7, 1, 7}}, {0});
        dsv2.lmd = {0, 0, 1, 1, 0, 0, 0, 0};
        dsv2.divergent = {false, false, false, false, false, false, false, false};
        dsv2.transReduceDSV();
        for (const auto& i : dsv2.divergent){
            std::cout << i << " ";
        }
        std::cout << std::endl;

        System dsv3 = {{0, 1, 2}, {{0, 1, 1}, {1, 0, 2}, {2, 0, 1}}, {0}};
        dsv3.lmd = {0, 0, 0};
        dsv3.divergent = {false, false, false};
        dsv3.transReduceDSV();
        for (const auto& i : dsv3.divergent){
            std::cout << i <<" ";
        }
        std::cout << std::endl;

        System dsv4({0, 1, 2, 3, 4, 5}, {{0, 0, 1}, {1, 0, 3}, {3, 0, 0},
                                              {0, 0, 2}, {2, 0, 3}, {2, 0, 4}, {4, 0, 5}, {3, 0, 5}}, {0});
        dsv4.lmd = {0, 0, 0, 0, 0, 0};
        dsv4.divergent = {false, false, false, false, false, false};
        dsv4.transReduceDSV();
        for (const auto &i : dsv4.divergent){
            std::cout << i << " ";
        }
        std::cout << std::endl;

        ///void tarjan(size_t s, size_t_in_vec dfn, size_t_in_vec low, std::stack<size_t> stack, size_t_in_vec inStack,
        ///            size_t_in_vec visited, size_t index = 0)
        //std::cout<< "Co-reach: " << std::endl;
        //System c1({0, 1, 2, 3, 4}, {{0, 1, 1}, {1, 2, 2}, {2, 4, 3}, {3, 6, 4}, {2, 3, 0}, {3, 5, 1}}, {3});
        //c1.coreach();
        //System c2({0, 1, 2, 3, 4, 5}, {{0, 1, 1}, {1, 4, 2}, {0, 3, 4}, {1, 2, 5}, {3, 3, 4}, {4, 4, 5}}, {0});
        //c2.coreach();
    }
    */
    //syncThenAbsEff(8);
    //syncThenAbsNew(10);
    //incSyncBufferNonDsv(3, 200);
    //incSyncBuffer(20, 3);
    ///clock_t t1 = clock();
    ///groupIn2(4);
    ///clock_t t2 = clock();
    ///std::cout<<"Incremental Approach Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
    //syncAndAbsEff(4);
    //clock_t t1 = clock();
    //testPerformance("dsv_test");
    //incSyncBuffer(50, 3);
    //incSyncBufferNonDsv(50, 3);
    //clock_t t2 = clock();
    //std::cout<<"inc sync buffer Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
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
//    testPerformance("G_dsv_test");
    //syncAndAbsNew(14);
    //syncThenAbsNew(7);

    //dpGroupTest2(8);
    //supDPtestNonInc(4);
    //syncAndAbsEff(4);
    //supDPtestNonInc(6);
    ///dpGroupTest2(8);
    //dpGroupTest2NonDSV(4);
    //dpTest2(22);
    ///dpTest2NonDSV(22);
    //syncThenAbsEff(4);
    //incSyncBuffer(20, 3);
    //incSyncBufferNonDsv(20, 3);
    ///reach test
    /*
    {
        std::set<size_t> result;
        System c1({0, 1, 2, 3, 4}, {{0, 1, 1}, {1, 2, 2}, {2, 4, 3}, {3, 6, 4}, {2, 3, 0}, {3, 5, 1}}, {0});

        result = c1.reachInT();
        for(const auto &iter : result){
            std::cout << iter << " ";
        }
        std::cout << std::endl;
        System c2({0, 1, 2, 3, 4, 5}, {{0, 1, 1}, {1, 4, 2}, {0, 3, 4}, {1, 2, 5}, {3, 3, 4}, {4, 4, 5}}, {4});
        result = c2.reachInT();
        for(const auto &iter : result){
            std::cout << iter << " ";
        }
        std::cout << std::endl;

        System c3({0, 1, 2, 3, 4, 5, 6, 7, 8}, {{0, 1, 1}, {0, 3, 5}, {1, 2, 6},
                                                     {1, 4, 2}, {4, 3, 5}, {5, 4, 6}, {7, 1, 8}}, {4});
        std::vector<size_t > r = c3.reachInTSetOp();
        for(const auto &iter : r){
            std::cout << iter << " ";
        }
        std::cout << std::endl;
    }
     */
    //reachPerformance("/Users/liangxudong/PycharmProjects/mat/reach_500**2.txt");
    /*
    System r0;
    r0.deltaR = {{2, 4, 7, 4, 5, 6, 7, 8, 9, 2}, {3, 5, 8, 0, 0, 0, 0, 0, 0, 0}, {0, 6, 9, 0, 0, 0, 0, 0, 0, 0}};
    r0.ndeltaR = {2, 3, 3, 1, 1, 1, 1, 1, 1, 1};
    r0.init = {10};
    r0.reachBoolean();

    System r2;
    r2.deltaR = {{2, 3, 4, 5, 6, 7, 7}, {3, 0, 0, 7, 0, 0, 0}};
    r2.ndeltaR = {2, 1, 1, 2, 1, 1, 1};
    r2.init = {1};
    r2.reachBoolean();

    System r3;
    r3.deltaR = {{2, 3, 4, 5, 6, 7, 7}, {3, 0, 0, 7, 0, 0, 0}};
    r3.ndeltaR = {2, 1, 1, 2, 1, 1, 1};
    r3.init = {2};
    r3.reachBoolean();

    System r4;
    r4.deltaR = {{2, 3, 4, 5, 6, 7, 2, 0, 0, 0}, {3, 0, 0, 7, 0, 0, 8, 0, 0, 0},
                 {0, 0, 0, 0, 0, 0, 9, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 10, 0, 0, 0}};
    r4.ndeltaR = {2, 1, 1, 2, 1, 1, 4, 0, 0, 0};
    r4.init = {6};
    r4.reachBoolean();
*/
    //System rt;
    //rt.ndeltaR = std::get<0> (reachTest2(500));
    //rt.deltaR = std::get<1> (reachTest2(500));
    //rt.init = {1};
    //rt.reachBoolean();
    GV gv1;
    gv1.states = {0, 1, 2, 3, 4};
    gv1.transitions = {{0, 2, 1}, {1, 0, 2}, {2, 0, 2}, {0, 0, 3}, {3, 0, 4}, {4, 0, 4}};
    gv1.reduceGV(1);
    return 0;
}
