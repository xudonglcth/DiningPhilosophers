//
// Created by 梁旭东 on 2020-02-24.
//

#ifndef DININGPHILOSOPHERS_DININGPHILOSOPHERSEFF_H
#define DININGPHILOSOPHERS_DININGPHILOSOPHERSEFF_H
#include "SyncInT.h"
void syncThenAbsEff(size_t n){
    std::vector<System> diningPhils, diningPhilPair;
    System::size_t_in_set taus = {};
    System G;
    size_t cnt = 2;
    //create n dining phils
    diningPhils = createDiningPhilosophersNew(n);
    //synch in pair
    for(size_t i = 0; i < 2*n; i += 2){
        //diningPhils[i].transInfo(); diningPhils[i + 1].transInfo();
        diningPhilPair.push_back(diningPhils[i].syncInT(diningPhils[i + 1]));
    }

    //sync all the pairs
    G = diningPhilPair[0];
    for (size_t i = 1; i < n; ++i){
        G.syncInT(diningPhilPair[i]);
        std::cout <<"Sync With G" << cnt ++<< "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    }
    //set all event labels in transitions to zero
    for (auto &t : G.transitions){
        t[1] = 0;
    }
    G.lmd[0] = 1;
    G.transReduce();
}

void syncAndAbsEff(size_t n){
    std::vector<System> diningPhils, diningPhilPair;
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
    //first step: sych first two and reduce
    G = diningPhilPair[0].syncInT(diningPhilPair[1]);
    std::cout <<"Sync With G" << cnt ++ <<"\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    taus.insert({1, 2*n+1});
    for (auto &t : G.transitions){
        if (taus.find(t[1]) != taus.end()){
            t[1] = 0;
        }
    }
    G.transReduce();
    std::cout << "After Reduction State#: " << G.states.size() << "\nAfter reduction Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    //std::cout << "After reduction: " << std::endl;
    //for (const auto &i : G.transitions){
    //    std::cout << i[0] << " " << i[1] << " "<< i[2] << std::endl;
    //}
    std::cout << std::endl;
    //second step: sync and reduce until the last pair
    for (size_t i = 2; i < n - 1; ++i){
        G.syncInT(diningPhilPair[i]);
        std::cout <<"Sync With G" << cnt ++ << "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
        taus.insert({2*i-2, 2*i-1, 2*n+2*i-2, 2*n+2*i-1});
        for (auto &t : G.transitions){
            if (taus.find(t[1]) != taus.end()){
                t[1] = 0;
            }
        }
        G.lmd[0] = 1;
        if (cnt == 6){
            std::ofstream outFile("/Users/liangxudong/Desktop/VLTS/NewFile.txt");
            outFile <<"des (" << 0 <<", " <<G.transitions.size() <<", " <<G.states.size() <<")\n";
            for (const auto &e : G.transitions) {
                if((e[1] == 0) && (e[0] != 0) && (e[2] != 0)){
                    outFile <<"(" << e[0] << ", "  <<"\""<< "tau" <<"\""<<", " << e[2] <<")"<< "\n";
                }
                else{
                    outFile <<"(" << e[0] << ", " << "\"" << e[1] << "\""<<", " << e[2] <<")"<< "\n";
                }
            }
        }
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
    G.syncInT(diningPhilPair[n - 1]);
    std::cout <<"Sync With G" << cnt ++ << "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    //at last all transitions are local and invisible.
    for (auto &t : G.transitions){
        t[1] = 0;
    }
    G.lmd[0] = 1;
    std::ofstream outFile("/Users/liangxudong/Desktop/VLTS/LastFile.txt");
    outFile <<"des (" << 0 <<", " <<G.transitions.size() <<", " <<G.states.size() <<")\n";
    for (const auto &e : G.transitions) {
        if((e[1] == 0) && (e[0] != 0) && (e[2] != 0)){
            outFile <<"(" << e[0] << ", "  <<"\""<< "tau" <<"\""<<", " << e[2] <<")"<< "\n";
        }
        else{
            outFile <<"(" << e[0] << ", " << "\"" << e[1] << "\""<<", " << e[2] <<")"<< "\n";
        }
    }
    clock_t t1 = clock();
    G.transReduce();
    clock_t t2 = clock();
    std::cout<<"One reduction run time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
    for (const auto &iter : G.transitions){
        std::cout << iter[0] << " "<< iter[1] << " " << iter[2] << std::endl;
    }
    //std::cout << "After Reduction State#: " << G.states.size() << "\nAfter Reduction Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    //std::cout << "After reduction: " << std::endl;
    //for (const auto &i : G.transitions){
    //    std::cout << i[0] << " " << i[1] << " "<< i[2] << std::endl;
    //}
    //std::cout << std::endl;

}
#endif //DININGPHILOSOPHERS_DININGPHILOSOPHERSEFF_H
