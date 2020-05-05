//
// Created by 梁旭东 on 2020-03-05.
//

#ifndef DININGPHILOSOPHERS_GROUPDINING_H
#define DININGPHILOSOPHERS_GROUPDINING_H

#include "SyncInT.h"
void groupIn2(size_t n){
    std::vector<System> diningPhils, diningPhilPair, diningGroup;
    size_t cnt = 1;
    System::size_t_in_set taus = {};
    System G;
    //create n dining phils
    diningPhils = createDiningPhilosophersNew(n);
    //synch in pair
    for(size_t i = 0; i < 2*n; i += 2){
        diningPhilPair.push_back(diningPhils[i].syncInT(diningPhils[i + 1]));
    }

    //sync and abstract
    //first step: sych first four and reduce
    G = diningPhilPair[0].syncInT(diningPhilPair[1]);
    std::cout <<"Sync of "
              << cnt <<", "<< cnt + 1
              <<"\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    cnt += 2;
    taus.insert({1, 2, 3,
                 2*n + 2, 2*n + 3, 2*n +4});
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
    //std::cout << std::endl;
    //second step: sync and reduce until the last pair
    //third step: sync with the last pair and reduce
    G.syncInT(diningPhilPair[n - 2]).syncInT(diningPhilPair[n - 1]);
    std::cout <<"Sync with "
              << cnt <<", "<< cnt + 1 <<", "<< cnt + 2 << ", and"<< cnt + 3
              << "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
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

void groupIn3(size_t n){
    std::vector<System> diningPhils, diningPhilPair, diningGroup;
    size_t cnt = 1;
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
    std::cout <<"Sync of "
            << cnt <<", "<< cnt + 1 <<", and "<< cnt + 2
            <<"\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    cnt += 3;
    taus.insert({1, 2, 3, 4, 5, 2*n + 2, 2*n + 3, 2*n +4, 2*n +5, 2*n +6});
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
    //std::cout << std::endl;
    //second step: sync and reduce until the last pair
    for (size_t i = 3; i < n - 3; i += 3){
        G.syncInT(diningPhilPair[i]).syncInT(diningPhilPair[i + 1]).syncInT(diningPhilPair[i + 2]);
        std::cout <<"Sync with "
            << cnt <<", "<< cnt + 1 <<", and "<< cnt + 2
            << "\nSync State#: "
            << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
        cnt += 3;
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
    std::cout <<"Sync with "
            << cnt <<", "<< cnt + 1 <<", and "<< cnt + 2
            << "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
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

void groupIn4(size_t n){
    std::vector<System> diningPhils, diningPhilPair, diningGroup;
    size_t cnt = 1;
    System::size_t_in_set taus = {};
    System G;
    //create n dining phils
    diningPhils = createDiningPhilosophersNew(n);
    //synch in pair
    for(size_t i = 0; i < 2*n; i += 2){
        diningPhilPair.push_back(diningPhils[i].syncInT(diningPhils[i + 1]));
    }

    //sync and abstract
    //first step: sych first four and reduce
    G = diningPhilPair[0].syncInT(diningPhilPair[1]).syncInT(diningPhilPair[2]).syncInT(diningPhilPair[3]);
    std::cout <<"Sync of "
              << cnt <<", "<< cnt + 1 <<", and "<< cnt + 2 << cnt + 3
              <<"\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    cnt += 4;
    taus.insert({1, 2, 3, 4, 5, 6, 7, 2*n + 2, 2*n + 3, 2*n +4, 2*n +5, 2*n +6, 2*n + 7, 2*n + 8});
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
    //std::cout << std::endl;
    //second step: sync and reduce until the last pair
    //third step: sync with the last pair and reduce
    G.syncInT(diningPhilPair[n - 4]).syncInT(diningPhilPair[n - 3]).syncInT(diningPhilPair[n - 2]).syncInT(diningPhilPair[n - 1]);
    std::cout <<"Sync with "
            << cnt <<", "<< cnt + 1 <<", "<< cnt + 2 << ", and"<< cnt + 3
              << "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
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

void groupIn5(size_t n){
    std::vector<System> diningPhils, diningPhilPair, diningGroup;
    size_t cnt = 1;
    System::size_t_in_set taus = {};
    System G;
    //create n dining phils
    diningPhils = createDiningPhilosophersNew(n);
    //synch in pair
    for(size_t i = 0; i < 2*n; i += 2){
        diningPhilPair.push_back(diningPhils[i].syncInT(diningPhils[i + 1]));
    }

    //sync and abstract
    //first step: sych first four and reduce
    G = diningPhilPair[0].syncInT(diningPhilPair[1]).syncInT(diningPhilPair[2]).syncInT(diningPhilPair[3]).syncInT(diningPhilPair[4]);
    std::cout <<"Sync of "
              << cnt <<", "<< cnt + 1 <<", and "<< cnt + 2 << cnt + 3 << cnt + 5
              <<"\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    cnt += 5;
    taus.insert({1, 2, 3, 4, 5, 6, 7, 8, 9,
                 2*n + 2, 2*n + 3, 2*n +4, 2*n +5, 2*n +6, 2*n + 7, 2*n + 8, 2*n + 9, 2*n + 10});
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
    //std::cout << std::endl;
    //second step: sync and reduce until the last pair
    //third step: sync with the last pair and reduce
    G.syncInT(diningPhilPair[n - 5]).syncInT(diningPhilPair[n - 4]).syncInT(diningPhilPair[n - 3]).syncInT(diningPhilPair[n - 2]).syncInT(diningPhilPair[n - 1]);
    std::cout <<"Sync with "
              << cnt <<", "<< cnt + 1 <<", "<< cnt + 2 << ", and"<< cnt + 3
              << "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
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
void syncIn5(){
    std::vector<System> diningPhils, diningPhilPair, diningGroup;
    diningPhils = createDiningPhilosophersNew(10);
    for(size_t i = 0; i < 20; i += 2){
        diningPhilPair.push_back(diningPhils[i].syncInT(diningPhils[i + 1]));
    }
    std::cout << "Starting sync the first 5" << std::endl;
    diningPhilPair[0].syncInT(diningPhilPair[1]).syncInT(diningPhilPair[2]).syncInT(diningPhilPair[3]).syncInT(diningPhilPair[4]);
    std::cout << "Starting sync the last 5" << std::endl;
    diningPhilPair[5].syncInT(diningPhilPair[6]).syncInT(diningPhilPair[7]).syncInT(diningPhilPair[8]).syncInT(diningPhilPair[9]);
    std::cout << "Sync two big models"<<std::endl;
    clock_t t1 = clock();
    diningPhilPair[0].syncInT(diningPhilPair[5]);
    clock_t t2 = clock();
    std::cout<<"One Big Sync Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
}
#endif //DININGPHILOSOPHERS_GROUPDINING_H
