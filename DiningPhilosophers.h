//
//  DiningPhilosophers.h
//  system
//
//  Created by 梁旭东 on 2020-02-12.
//  Copyright © 2020 xudongl. All rights reserved.
//

#ifndef DiningPhilosophers_h
#define DiningPhilosophers_h
#include "System.h"
std::vector<System> createDiningPhilosophers(size_t n){
    System::size_t_in_vec r_up, r_down, eat, think;
    std::vector<System> Gv, Rv, diningPhilosophersV;
    size_t i = 0, j = 0;
    System G, R;
    for (int i = 1; i <= n; ++i){
        r_up.push_back(i);
        r_down.push_back(n + i);
        eat.push_back(2*n + 1);
        think.push_back(3*n + i);
    }
    for (int i = 0; i != n; ++i){
        Gv.push_back(G);
        Rv.push_back(R);
    }
    for (auto &iter : Gv){
        iter.states = {0, 1, 2, 3, 4, 5, 6};
        iter.transitions = {
                {0, think[i], 0},
                {0, r_up[i], 2},
                {0, r_up[(i + 1 == n) ? 1 : i + 1], 3},
                {1, r_down[(i + 1 == n) ? 1 : i + 1], 0},
                {2, r_up[(i + 1 == n) ? 1 : i + 1], 5},
                {3, r_up[i], 5},
                {4, r_down[i], 0},
                {5, eat[i], 6},
                {6, r_down[i], 1},
                {6, r_down[(i + 1 == n) ? 1 : i + 1], 4}
        };
        ++i;
        iter.lmd = {1, 0, 0, 0, 0, 0, 0};
        iter.init = {0};
    }
    for (auto & iter : Rv){
        iter.states = {0, 1};
        iter.transitions = {
                {0, r_up[j], 1},
                {1, r_down[j], 0}
        };
        ++j;
        iter.init = {0};
        iter.lmd = {1, 0};
    }

    for (size_t k = 0; k != n; ++k){
        diningPhilosophersV.push_back(Gv[k]);
        diningPhilosophersV.push_back(Rv[k]);
    }
    return diningPhilosophersV;
}

void syncAndAbstract(size_t n){
    std::vector<System> diningPhils, diningPhilPair;
    System G;
    //create n dining phils
    diningPhils = createDiningPhilosophers(n);
    //synch in pair
    for(size_t i = 0; i < 2*n; i += 2){
        diningPhilPair.push_back(diningPhils[i].sync(diningPhils[i + 1]));
    }
    //sync and abstract
    //first step: sych first two and reduce
    G = diningPhilPair[0].sync(diningPhilPair[1]);
    std::cout << "Sync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<'\n'<< std::endl;
    G.tau = {2, n+2, 2*n+1, 2*n+2, 3*n+1, 3*n+2};
    G.transReduce();
    std::cout << "After Reduction State#: " << G.states.size() << "\nAfter reduction Trans#: " << G.transitions.size() <<'\n'<< std::endl;
    //std::cout << "After reduction: " << std::endl;
    //for (const auto &i : G.transitions){
    //    std::cout << i[0] << " " << i[1] << " "<< i[2] << std::endl;
    //}
    //std::cout << std::endl;
    //second step: sync and reduce until the last pair
    for (size_t i = 2; i < n - 1; ++i){
        G.sync(diningPhilPair[i]);
        std::cout << "Sync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() << '\n'<<std::endl;
        G.tau.insert({i+1, n+i+1, 2*n+i+1, 3*n+i+1});
        G.lmd[0] = 1;
        G.transReduce();
        std::cout << "After Reduction State#: " << G.states.size() << "\nAfter Reduction Trans#: " << G.transitions.size() <<'\n'<< std::endl;
        //std::cout << "After reduction: " << std::endl;
        //for (const auto &i : G.transitions){
        //    std::cout << i[0] << " " << i[1] << " "<< i[2] << std::endl;
        //}
        //std::cout << std::endl;
    }
    //third step: sync with the last pair and reduce
    G.sync(diningPhilPair[n - 1]);
    std::cout << "Sync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() << '\n'<<std::endl;
    G.tau.insert({1, n, 2*n, n+1, 3*n, 4*n});
    G.lmd[0] = 1;
    G.transReduce();
    std::cout << "After Reduction State#: " << G.states.size() << "\nAfter Reduction Trans#: " << G.transitions.size() <<'\n'<< std::endl;
    //std::cout << "After reduction: " << std::endl;
    //for (const auto &i : G.transitions){
    //    std::cout << i[0] << " " << i[1] << " "<< i[2] << std::endl;
    //}
    //std::cout << std::endl;
}

#endif /* DiningPhilosophers_h */
