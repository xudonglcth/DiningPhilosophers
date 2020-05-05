//
// Created by 梁旭东 on 2020-02-20.
//

#ifndef DININGPHILOSOPHERS_DININGPHILOSOPHERSNEW_H
#define DININGPHILOSOPHERS_DININGPHILOSOPHERSNEW_H

#include <iostream>
#include "System.h"
std::vector<System> createDiningPhilosophersNew(size_t n){
    System::size_t_in_vec r_up, r_down;
    System G;
    std::vector<System> Gv, Rv, diningPhilosophersV;
    size_t iter_Gn = 0, iter_Rn = 0, cnt = 0;
    for (size_t i = 0; i < n; ++i){
        Gv.push_back(G); Rv.push_back(G);
    }

    for (; cnt <= 2*n; ++cnt){
        r_up.push_back(cnt);
    }

    for (; cnt <= 4*n + 1; ++cnt){
        r_down.push_back(cnt);
    }

    for (auto &iter : Gv){
        iter.states = {0, 1, 2, 3, 4, 5};
        iter.transitions = {
                {0, r_up[(iter_Gn == 0) ? 2*n : 2 * iter_Gn],     1},
                {0, r_up[2 * iter_Gn + 1],                        2},
                {1, r_up[2 * iter_Gn + 1],                        3},
                {2, r_up[(iter_Gn == 0) ? 2*n : 2 * iter_Gn],     3},
                {3, r_down[(iter_Gn == 0) ? 2*n : 2 * iter_Gn],   4},
                {3, r_down[2 * iter_Gn + 1],                      5},
                {4, r_down[2 * iter_Gn + 1],                      0},
                {5, r_down[(iter_Gn == 0) ? 2*n : 2 * iter_Gn],   0}
        };
        iter.lmd = {1, 0, 0, 0, 0, 0};
        iter.init = {0};
        ++iter_Gn;
    }

    for (auto &iter : Rv){
        iter.states = {0, 1, 2};
        iter.transitions = {
                {0, r_up[2 * iter_Rn + 1],   1},
                {1, r_down[2 * iter_Rn + 1], 0},
                {0, r_up[2 * iter_Rn + 2],   2},
                {2, r_down[2 * iter_Rn + 2], 0}
        };
        iter.lmd = {1, 0, 0};
        iter.init = {0};
        ++iter_Rn;
    }
    /*
    std::cout << "Up: "<<std::endl;
    for (auto const iter : r_up){
        std::cout << iter <<" ";
    }
    std::cout << "\n";

    std::cout << "Down: "<<std::endl;
    for (auto const iter : r_down){
        std::cout << iter <<" ";
    }
    std::cout << "\n";

    std::cout << "G" <<std::endl;
    for (const auto &iter : Gv){
        for (const auto &t : iter.transitions){
            std::cout << t[0]<<" "<<t[1]<<" "<<t[2]<<std::endl;
        }
        std::cout << '\n';
    }

    std::cout << "R" << std::endl;
    for (const auto &iter : Rv){
        for (const auto &t : iter.transitions){
            std::cout << t[0]<<" "<<t[1]<<" "<<t[2]<<std::endl;
        }
        std::cout << '\n';
    }
     */
    for (size_t i = 0; i < n; ++i){
        diningPhilosophersV.push_back(Gv[i]);
        diningPhilosophersV.push_back(Rv[i]);
    }

    return diningPhilosophersV;
}

std::vector<System> createDiningPhilosophers2(size_t n){
    System::size_t_in_vec take, put, eat;
    System G;
    std::vector<System> Gv, Rv, diningPhilosophersV;
    size_t iter_Gn = 0, iter_Rn = 0, cnt = 0;
    for (size_t i = 0; i < n; ++i){
        Gv.push_back(G); Rv.push_back(G);
    }


    while(cnt < 2*n){
        take.push_back(++cnt);
        put.push_back(2*n + cnt);
        eat.push_back(4*n + cnt);
    }

    for (auto &iter : Gv){
        iter.states = {0, 1, 2, 3, 4, 5, 6};
        iter.transitions = {
                {0, take[2 * iter_Gn + 1], 2},
                {0, take[2 * iter_Gn],     1},
                {1, take[2 * iter_Gn + 1], 3},
                {2, take[2 * iter_Gn],     3},
                {3, eat[iter_Gn],          4},
                {4, put[2 * iter_Gn + 1],  6},
                {4, put[2 * iter_Gn],      5},
                {5, put[2 * iter_Gn + 1],  0},
                {6, put[2 * iter_Gn],      0}
        };
        iter.lmd = {1, 0, 0, 0, 0, 0, 0};
        iter.init = {0};
        ++iter_Gn;
    }

    for (auto &iter : Rv){
        iter.states = {0, 1};
        iter.transitions = {
                {0, take[2 * iter_Rn],                                1},
                {1, put[2 * iter_Rn],                                 0},
                {0, take[(iter_Rn == 0) ? 2*n - 1 : 2 * iter_Rn - 1],   1},
                {1, put[(iter_Rn == 0) ? 2*n - 1 : 2 * iter_Rn - 1],    0}
        };
        iter.lmd = {1, 0};
        iter.init = {0};
        ++iter_Rn;
    }
    for (size_t i = 0; i < n; ++i){
        diningPhilosophersV.push_back(Gv[i]);
        diningPhilosophersV.push_back(Rv[i]);
    }

    return diningPhilosophersV;
}

void syncAndAbsNew(size_t n){
    std::vector<System> diningPhils, diningPhilPair;
    size_t cnt = 2;
    System::size_t_in_set taus = {};
    System G;
    //create n dining phils
    diningPhils = createDiningPhilosophersNew(n);
    //synch in pair
    for(size_t i = 0; i < 2*n; i += 2){
        diningPhilPair.push_back(diningPhils[i].sync(diningPhils[i + 1]));
    }

    //sync and abstract
    //first step: sych first two and reduce
    G = diningPhilPair[0].sync(diningPhilPair[1]);
    std::cout <<"Sync With G" << cnt ++ <<"\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    taus.insert({1, 2*n+1});
    for (auto &t : G.transitions){
        if (taus.find(t[1]) != taus.end()){
            t[1] = 0;
        }
    }
    G.transReduceDSV();
    std::cout << "After Reduction State#: " << G.states.size() << "\nAfter reduction Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    //std::cout << "After reduction: " << std::endl;
    //for (const auto &i : G.transitions){
    //    std::cout << i[0] << " " << i[1] << " "<< i[2] << std::endl;
    //}
    std::cout << std::endl;
    //second step: sync and reduce until the last pair
    for (size_t i = 2; i < n - 1; ++i){
        G.sync(diningPhilPair[i]);
        std::cout <<"Sync With G" << cnt ++ << "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
        taus.insert({2*i-2, 2*i-1, 2*n+2*i-2, 2*n+2*i-1});
        for (auto &t : G.transitions){
            if (taus.find(t[1]) != taus.end()){
                t[1] = 0;
            }
        }
        G.lmd[0] = 1;
        G.transReduceDSV();
        std::cout << "After Reduction State#: " << G.states.size() << "\nAfter Reduction Trans#: " << G.transitions.size() <<"\n"<< std::endl;
        //std::cout << "After reduction: " << std::endl;
        //for (const auto &i : G.transitions){
        //    std::cout << i[0] << " " << i[1] << " "<< i[2] << std::endl;
        //}
        std::cout << std::endl;
    }
    //third step: sync with the last pair and reduce
    G.sync(diningPhilPair[n - 1]);
    std::cout <<"Sync With G" << cnt ++ << "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    //at last all transitions are local and invisible.
    for (auto &t : G.transitions){
        t[1] = 0;
    }
    G.lmd[0] = 1;
    G.transReduceDSV();
    for (const auto &iter : G.transitions){
        std::cout << iter[0] << " "<< iter[1] << " " << iter[2] << std::endl;
    }
    std::cout << "After Reduction State#: " << G.states.size() << "\nAfter Reduction Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    //std::cout << "After reduction: " << std::endl;
    //for (const auto &i : G.transitions){
    //    std::cout << i[0] << " " << i[1] << " "<< i[2] << std::endl;
    //}
    //std::cout << std::endl;

}

void syncThenAbsNew(size_t n){
    std::vector<System> diningPhils, diningPhilPair;
    System::size_t_in_set taus = {};
    System G;
    size_t cnt = 2;
    //create n dining phils
    diningPhils = createDiningPhilosophersNew(n);
    //synch in pair
    for(size_t i = 0; i < 2*n; i += 2){
        //diningPhils[i].transInfo(); diningPhils[i + 1].transInfo();
        diningPhilPair.push_back(diningPhils[i].sync(diningPhils[i + 1]));
    }

    //sync all the pairs
    G = diningPhilPair[0];
    for (size_t i = 1; i < n; ++i){
        G.sync(diningPhilPair[i]);
        std::cout <<"Sync With G" << cnt ++<< "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    }
    //set all event labels in transitions to zero
    for (auto &t : G.transitions){
        t[1] = 0;
    }
    G.lmd[0] = 1;
    /*
    G.transitions.clear();
    for (const auto & t : G.delta){
        for (const auto & i : t.second){
            G.transitions.push_back({t.first[0], 0, i});
        }
    }
     */
    //clock_t t_1 = clock();
    //G.transReduceDSV();
    //clock_t t_2 = clock();
    //std::cout<<"Brute-Force Approach Abstraction Run Time: "<<(double)(t_2 - t_1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
    //for (const auto &iter : G.transitions){
    //    std::cout << iter[0] << " "<< iter[1] << " " << iter[2] << std::endl;
    //}
    //std::cout << "After Reduction State#: " << G.states.size() << "\nAfter Reduction Trans#: " << G.transitions.size() <<"\n"<< std::endl;

}
#endif //DININGPHILOSOPHERS_DININGPHILOSOPHERSNEW_H
