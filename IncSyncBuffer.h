//
// Created by 梁旭东 on 2020-03-10.
//

#ifndef DININGPHILOSOPHERS_INCSYNCBUFFER_H
#define DININGPHILOSOPHERS_INCSYNCBUFFER_H

#include "System.h"
void incSyncBuffer(size_t n_cap, size_t n_buffer){
    System::size_t_in_set sigma_local;
    System G1, G2;
    size_t a = 5, mark = 1, model = 1;
    sigma_local = G1.TSmodel(model, a, n_cap, mark);
    for (auto & iter : G1.transitions){
        if (sigma_local.find(iter[1]) != sigma_local.end()){
            iter[1] = 0;
        }
    }
    std::cout << "#state: " << G1.states.size() <<" #trans: "<< G1.transitions.size() << std::endl;
    clock_t t1 = clock();
    G1.transReduceDSV2();
    clock_t t2 = clock();
    std::cout << "#reduced state: " << G1.states.size() <<" #reduced trans: "<< G1.transitions.size() <<std::endl;
    std::cout<<"One Sigref Reduction Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
    ++a;
    for(size_t i = 1; i <= n_buffer; ++i){
        model = 2;
        sigma_local = G2.TSmodel(model, a, n_cap, mark);
        ++a;
        G1.syncInT(G2);
        for (auto & iter : G1.transitions){
            if (sigma_local.find(iter[1]) != sigma_local.end()){
                iter[1] = 0;
            }
        }
        std::cout << "#state: " << G1.states.size() <<" #trans: "<< G1.transitions.size() << std::endl;



        t1 = clock();
        G1.transReduceDSV2();
        t2 = clock();
        std::cout << "#reduced state: " << G1.states.size() <<" #reduced trans: "<< G1.transitions.size() <<std::endl;
        std::cout<<"One Sigref Reduction Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
        std::cout << "#reduced state: " << G1.states.size() << std::endl;
        model = 3;
        if(i == n_buffer) model = 4;
        sigma_local = G2.TSmodel(model, a, n_cap, mark);
        if(i < n_buffer) ++a;
        G1.syncInT(G2);
        if(i < n_buffer){
            for (auto & iter : G1.transitions){
                if (sigma_local.find(iter[1]) != sigma_local.end()){
                    iter[1] = 0;
                }
            }
            std::cout << "#state: " << G1.states.size() <<" #trans: "<< G1.transitions.size() << std::endl;
            t1 = clock();
            G1.transReduceDSV2();
            t2 = clock();
            std::cout << "#reduced state: " << G1.states.size() <<" #reduced trans: "<< G1.transitions.size() <<std::endl;
            std::cout<<"One Sigref Reduction Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
        }

    }
    model = 5;
    sigma_local = G2.TSmodel(model, a, n_cap, mark);
    G1.syncInT(G2);
    for (auto & iter : G1.transitions){
        if (sigma_local.find(iter[1]) != sigma_local.end()){
            iter[1] = 0;
        }
    }
    std::cout << "#state: " << G1.states.size() <<" #trans: "<< G1.transitions.size() << std::endl;
    t1 = clock();
    G1.transReduceDSV2();
    t2 = clock();
    std::cout << "#reduced state: " << G1.states.size() <<" #reduced trans: "<< G1.transitions.size() <<std::endl;
    std::cout<<"One Sigref Reduction Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
    std::cout <<G1.states.size()<<", "<< G1.transitions.size() << std::endl;
    for (const auto &iter : G1.transitions){
        std::cout << iter[0] << iter[1] << iter[2] << std::endl;
    }
}



void incSyncBufferNonDsv(size_t n_cap, size_t n_buffer){
    std::string fileName;
    System::size_t_in_set sigma_local;
    System G1, G2;
    size_t a = 5, mark = 1, model = 1;
    sigma_local = G1.TSmodel(model, a, n_cap, mark);
    for (auto & iter : G1.transitions){
        if (sigma_local.find(iter[1]) != sigma_local.end()){
            iter[1] = 0;
        }
    }
    std::cout << "#state: " << G1.states.size() <<" #trans: "<< G1.transitions.size() << std::endl;
    fileName = "/Users/liangxudong/Desktop/VLTS/inc_test/inc_" + std::to_string(G1.states.size()) + ".txt";
    std::ofstream outFile(fileName);
    outFile <<"des (" << 0 <<", " <<G1.transitions.size() <<", " <<G1.states.size() <<")\n";
    for (const auto &e : G1.transitions) {
        if((e[1] == 0) && (e[0] != 0) && (e[2] != 0)){
            outFile <<"(" << e[0] << ", "  <<"\""<< "tau" <<"\""<<", " << e[2] <<")"<< "\n";
        }
        else{
            outFile <<"(" << e[0] << ", " << "\"" << e[1] << "\""<<", " << e[2] <<")"<< "\n";
        }
    }
    clock_t t1 = clock();
    G1.transReduce();
    clock_t t2 = clock();
    std::cout << "#reduced state: " << G1.states.size() <<" #reduced trans: "<< G1.transitions.size() <<std::endl;
    std::cout<<"One Sigref Reduction Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
    ++a;
    for(size_t i = 1; i <= n_buffer; ++i){
        model = 2;
        sigma_local = G2.TSmodel(model, a, n_cap, mark);
        ++a;
        G1.syncInT(G2);
        for (auto & iter : G1.transitions){
            if (sigma_local.find(iter[1]) != sigma_local.end()){
                iter[1] = 0;
            }
        }
        std::cout << "#state: " << G1.states.size() <<" #trans: "<< G1.transitions.size() << std::endl;

        fileName = "/Users/liangxudong/Desktop/VLTS/inc_test/inc_" + std::to_string(G1.states.size()) + ".txt";
        std::ofstream outFile(fileName);
        outFile <<"des (" << 0 <<", " <<G1.transitions.size() <<", " <<G1.states.size() <<")\n";
        for (const auto &e : G1.transitions) {
            if((e[1] == 0) && (e[0] != 0) && (e[2] != 0)){
                outFile <<"(" << e[0] << ", "  <<"\""<< "tau" <<"\""<<", " << e[2] <<")"<< "\n";
            }
            else{
                outFile <<"(" << e[0] << ", " << "\"" << e[1] << "\""<<", " << e[2] <<")"<< "\n";
            }
        }
        t1 = clock();
        G1.transReduce();
        t2 = clock();
        std::cout << "#reduced state: " << G1.states.size() <<" #reduced trans: "<< G1.transitions.size() <<std::endl;
        std::cout<<"One Sigref Reduction Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;

        model = 3;
        if(i == n_buffer) model = 4;
        sigma_local = G2.TSmodel(model, a, n_cap, mark);
        if(i < n_buffer) ++a;
        G1.syncInT(G2);
        if(i < n_buffer){
            for (auto & iter : G1.transitions){
                if (sigma_local.find(iter[1]) != sigma_local.end()){
                    iter[1] = 0;
                }
            }
            std::cout << "#state: " << G1.states.size() <<" #trans: "<< G1.transitions.size() << std::endl;
            fileName = "/Users/liangxudong/Desktop/VLTS/inc_test/inc_" + std::to_string(G1.states.size()) + ".txt";
            std::ofstream outFile(fileName);
            outFile <<"des (" << 0 <<", " <<G1.transitions.size() <<", " <<G1.states.size() <<")\n";
            for (const auto &e : G1.transitions) {
                if((e[1] == 0) && (e[0] != 0) && (e[2] != 0)){
                    outFile <<"(" << e[0] << ", "  <<"\""<< "tau" <<"\""<<", " << e[2] <<")"<< "\n";
                }
                else{
                    outFile <<"(" << e[0] << ", " << "\"" << e[1] << "\""<<", " << e[2] <<")"<< "\n";
                }
            }
            t1 = clock();
            G1.transReduce();
            t2 = clock();
            std::cout << "#reduced state: " << G1.states.size() <<" #reduced trans: "<< G1.transitions.size() <<std::endl;
            std::cout<<"One Sigref Reduction Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;

        }

    }
    model = 5;
    sigma_local = G2.TSmodel(model, a, n_cap, mark);
    G1.syncInT(G2);
    for (auto & iter : G1.transitions){
        if (sigma_local.find(iter[1]) != sigma_local.end()){
            iter[1] = 0;
        }
    }
    std::cout << "#state: " << G1.states.size() <<" #trans: "<< G1.transitions.size() << std::endl;
    fileName = "/Users/liangxudong/Desktop/VLTS/inc_test/inc_" + std::to_string(G1.states.size()) + ".txt";
    std::ofstream outFile1(fileName);
    outFile1 <<"des (" << 0 <<", " <<G1.transitions.size() <<", " <<G1.states.size() <<")\n";
    for (const auto &e : G1.transitions) {
        if((e[1] == 0) && (e[0] != 0) && (e[2] != 0)){
            outFile <<"(" << e[0] << ", "  <<"\""<< "tau" <<"\""<<", " << e[2] <<")"<< "\n";
        }
        else{
            outFile <<"(" << e[0] << ", " << "\"" << e[1] << "\""<<", " << e[2] <<")"<< "\n";
        }
    }
    t1 = clock();
    G1.transReduce();
    t2 = clock();
    std::cout << "#reduced state: " << G1.states.size() <<" #reduced trans: "<< G1.transitions.size() <<std::endl;
    std::cout<<"One Sigref Reduction Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;

}


System::size_t_in_set System::TSmodel(size_t model, size_t a, size_t ncap, size_t mark){
    size_t_in_set sigma_local;
    transitions.clear(); states.clear();lmd.clear();
    if(model == 1){
        transitions = {{0, a, 1}, {1, a + 1, 0}};
        sigma_local = {a};
        states = {0, 1};
        lmd = {0, 0};
    }
    else if(model == 2){
        for (size_t i = 0; i < ncap; ++i){
            transitions.push_back({i, a, i + 1});
            transitions.push_back({i + 1, a + 1, i});
            states.push_back(i);
            lmd.push_back(0);
        }
        sigma_local = {a};
        states.push_back(ncap + 1);
    }
    else if(model == 3){
        transitions = {{0, a, 1}, {1, a + 1, 0}};
        sigma_local = {a};
        states = {0, 1};
        lmd = {0, 0};
    }
    else if(model == 4){
        transitions = {{0, a, 1}, {1, a + 1, 0}};
        sigma_local = {};
        states = {0, 1};
        lmd = {0, 0};
    }
    else if(model == 5){
        transitions = {{0, a, 1}, {1, a+1, 0}};
        //transitions = {{0, a, 1}, {1, a+1, 0}, {1, a+2, 0}};
        sigma_local = {a, a+1};
        //sigma_local = {a, a+1, a+2};
        states = {0, 1};
        lmd = {0, 0};
    }
    if (mark == 1){
        lmd[0] = 1;
    }
    init = {0};
    return sigma_local;
}
#endif //DININGPHILOSOPHERS_INCSYNCBUFFER_H
