//
// Created by 梁旭东 on 2020-03-05.
//

#ifndef DININGPHILOSOPHERS_GROUPDINING_H
#define DININGPHILOSOPHERS_GROUPDINING_H

#include "SyncInT.h"
void groupIn3(size_t n){
    std::vector<System> diningPhils, diningPhilPair, diningGroup;
    size_t cnt = 2;
    System::size_t_in_set taus = {};
    System G;
    //create n dining phils
    diningPhils = createDiningPhilosophersNew(n);
    //synch in pair
    for(size_t i = 0; i < 2*n; i += 2){
        diningPhilPair.push_back(diningPhils[i].syncInT(diningPhils[i + 1]));
    }

    //sync and abstract
    //first step: sych first three and reduce
    G = diningPhilPair[0].syncInT(diningPhilPair[1]).syncInT(diningPhilPair[2]);
    std::cout <<"Sync With G" << cnt ++ <<"\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    taus.insert({1, 2, 3, 4, 5, 2*n + 2, 2*n + 3, 2*n +4, 2*n +5, 2*n +6});
    for (auto &t : G.transitions){
        if (taus.find(t[1]) != taus.end()){
            t[1] = 0;
        }
    }
    G.transReduce();
    //std::cout << "After Reduction State#: " << G.states.size() << "\nAfter reduction Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    //std::cout << "After reduction: " << std::endl;
    //for (const auto &i : G.transitions){
    //    std::cout << i[0] << " " << i[1] << " "<< i[2] << std::endl;
    //}
    //std::cout << std::endl;
    //second step: sync and reduce until the last pair
    for (size_t i = 3; i < n - 1; i += 3){
        G.syncInT(diningPhilPair[i]).syncInT(diningPhilPair[i + 1]).syncInT(diningPhilPair[i + 2]);
        std::cout <<"Sync With G" << cnt ++ << "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
        taus.insert({2*i, 2*i+1, 2*i+2, 2*i+3, 2*i+4, 2*i+5, 2*n+1+2*i, 2*n+1+2*i+1, 2*n+1+2*i+2, 2*n+1+2*i+3, 2*n+1+2*i+4, 2*n+1+2*i+5});
        for (auto &t : G.transitions){
            if (taus.find(t[1]) != taus.end()){
                t[1] = 0;
            }
        }
        G.lmd[0] = 1;
        clock_t t1 = clock();
        G.transReduce();
        clock_t t2 = clock();
        std::cout<<"One reduction run time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;

        std::cout << "After Reduction State#: " << G.states.size() << "\nAfter Reduction Trans#: " << G.transitions.size() <<"\n"<< std::endl;
        //std::cout << "After reduction: " << std::endl;
        //for (const auto &i : G.transitions){
        //    std::cout << i[0] << " " << i[1] << " "<< i[2] << std::endl;
        //}
        std::cout << std::endl;
    }
    //third step: sync with the last pair and reduce
    G.syncInT(diningPhilPair[n - 3]).syncInT(diningPhilPair[n - 2]).syncInT(diningPhilPair[n - 1]);
    std::cout <<"Sync With G" << cnt ++ << "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    //at last all transitions are local and invisible.
    for (auto &t : G.transitions){
        t[1] = 0;
    }
    G.lmd[0] = 1;
    clock_t t1 = clock();
    G.transReduce();
    clock_t t2 = clock();
    std::cout<<"One reduction run time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
    for (const auto &iter : G.transitions){
        std::cout << iter[0] << " "<< iter[1] << " " << iter[2] << std::endl;
    }
}
#endif //DININGPHILOSOPHERS_GROUPDINING_H
