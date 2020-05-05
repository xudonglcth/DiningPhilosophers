//
// Created by 梁旭东 on 2020-03-23.
//

#ifndef DININGPHILOSOPHERS_DININGTEST2_H
#define DININGPHILOSOPHERS_DININGTEST2_H

#include "System.h"
#include "DiningPhilosophersNew.h"
#include "SyncInT.h"
void dpTest2(size_t n){
    std::vector<System> diningPhils, diningPhilPair;
    size_t cnt = 2;
    System::size_t_in_set taus = {};
    System G;
    //create n dining phils
    diningPhils = createDiningPhilosophers2(n);
    //synch in pair
    for(size_t i = 0; i < 2*n; i += 2){
        diningPhilPair.push_back(diningPhils[i].syncInT(diningPhils[i + 1]));
    }
    //sync and abstract
    //first step: sych first two and reduce
    G = diningPhilPair[0].syncInT(diningPhilPair[1]);
    std::cout <<"Sync With G" << cnt ++ <<"\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() << std::endl;
    taus.insert({1, 2, 3, 2*n + 1, 2*n + 2, 2*n + 3, 4*n + 1, 4*n + 2});
    for (auto &t : G.transitions){
        if (taus.find(t[1]) != taus.end()){
            t[1] = 0;
        }
    }
    /*
    std::ofstream outFile("/Users/liangxudong/Desktop/VLTS/G_dsv_test.aut");
    outFile <<"des (" << 0 <<", " <<G.transitions.size() <<", " <<G.states.size() <<")\n";
    for (const auto &e : G.transitions) {
        if((e[1] == 0) && (e[0] != 0) && (e[2] != 0)){
            outFile <<"(" << e[0] << ", "  <<"\""<< "tau" <<"\""<<", " << e[2] <<")"<< "\n";
        }
        else{
            outFile <<"(" << e[0] << ", " << "\"" << e[1] << "\""<<", " << e[2] <<")"<< "\n";
        }
    }
     */
    G.transReduceDSV2();
    std::cout << "After Reduction State#: " << G.states.size() << "\nAfter reduction Trans#: " << G.transitions.size() << std::endl;
    std::cout << std::endl;
    //second step: sync and reduce until the last pair
    for (size_t i = 2; i < n - 1; ++i){
        G.syncInT(diningPhilPair[i]);
        std::cout <<"Sync With G" << cnt ++ << "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() << std::endl;
        taus.insert({2*i, 2*i + 1, 2*i + 2*n, 2*i + 1 + 2*n, 4*n + i + 1});
        for (auto &t : G.transitions){
            if (taus.find(t[1]) != taus.end()){
                t[1] = 0;
            }
        }
        if (i == 35){
        std::ofstream outFile("/Users/liangxudong/Desktop/VLTS/G_dsv_test.aut");
        outFile <<"des (" << 0 <<", " <<G.transitions.size() <<", " <<G.states.size() <<")\n";
        for (const auto &e : G.transitions) {
            if((e[1] == 0) && (e[0] != 0) && (e[2] != 0)){
                outFile <<"(" << e[0] << ", "  <<"\""<< "tau" <<"\""<<", " << e[2] <<")"<< "\n";
            }
            else if((e[1] == 0) && (e[0] != 0) && (e[2] == 0)){
                outFile <<"(" << e[0] << ", "  <<"\""<< "50002" <<"\""<<", " << e[2] <<")"<< "\n";
            }
            else if((e[1] == 0) && (e[0] == 0) && (e[2] != 0)){
                outFile <<"(" << e[0] << ", "  <<"\""<< "50001" <<"\""<<", " << e[2] <<")"<< "\n";
            }
            else{
                outFile <<"(" << e[0] << ", " << "\"" << e[1] << "\""<<", " << e[2] <<")"<< "\n";
            }
        }
        }
        System Gt;
        Gt.states = G.states; Gt.transitions = G.transitions;
        G = Gt;
        clock_t t1 = clock();
        G.transReduceDSV2();
        clock_t t2 = clock();
        std::cout << "After Reduction State#: " << G.states.size() << "\nAfter Reduction Trans#: " << G.transitions.size() << std::endl;
        std::cout<<"One reduction run time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s"<<std::endl;
        std::cout << std::endl;
    }
    //third step: sync with the last pair and reduce
    G.syncInT(diningPhilPair[n - 1]);
    std::cout <<"Sync With G" << cnt ++ << "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() << std::endl;
    //at last all transitions are local and invisible.
    for (auto &t : G.transitions){
        t[1] = 0;
    }
    //G.lmd[0] = 1;
    G.transReduceDSV2();
    for (const auto &iter : G.transitions){
        std::cout << iter[0] << " "<< iter[1] << " " << iter[2] << std::endl;
    }
}

void supDPtestNonInc(size_t n){
    std::vector<System> diningPhils, diningPhilPair;
    System::size_t_in_set taus = {};
    System G;
    size_t cnt = 2;
    //create n dining phils
    diningPhils = createDiningPhilosophers2(n);
    //synch in pair
    for(size_t i = 0; i < 2*n; i += 2){
        //diningPhils[i].transInfo(); diningPhils[i + 1].transInfo();
        /*
        std::string s1 = "/Users/liangxudong/Desktop/VLTS/dp", s2 = ".txt";
        std::string fileName =  s1 + std::to_string(i) +s2;
        std::ofstream outFile(fileName);
        outFile <<"des (" << 0 <<", " <<diningPhils[i].transitions.size() <<", " <<diningPhils[i].states.size() <<")\n";
        for (const auto &e : diningPhils[i].transitions) {
            if((e[1] == 0) && (e[0] != 0) && (e[2] != 0)){
                outFile <<"(" << e[0] << ", "  <<"\""<< "tau" <<"\""<<", " << e[2] <<")"<< "\n";
            }
            else{
                outFile <<"(" << e[0] << ", " << "\"" << e[1] << "\""<<", " << e[2] <<")"<< "\n";
            }
        }

        std::string s3 = "/Users/liangxudong/Desktop/VLTS/dp", s4 = ".txt";
        std::string fileName2 =  s3 + std::to_string(i + 1) +s4;
        std::ofstream outFile2(fileName2);
        outFile2 <<"des (" << 0 <<", " <<diningPhils[i+1].transitions.size() <<", " <<diningPhils[i+1].states.size() <<")\n";
        for (const auto &e : diningPhils[i+1].transitions) {
            if((e[1] == 0) && (e[0] != 0) && (e[2] != 0)){
                outFile2 <<"(" << e[0] << ", "  <<"\""<< "tau" <<"\""<<", " << e[2] <<")"<< "\n";
            }
            else{
                outFile2 <<"(" << e[0] << ", " << "\"" << e[1] << "\""<<", " << e[2] <<")"<< "\n";
            }
        }
    */
        diningPhilPair.push_back(diningPhils[i].syncInT(diningPhils[i + 1]));
    }

    //sync all the pairs except the last one
    G = diningPhilPair[0];
    for (size_t i = 1; i < n - 1; ++i){
        G.syncInT(diningPhilPair[i]);
        std::cout <<"Sync With G" << cnt ++<< "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() <<"\n"<< std::endl;
    }
    std::cout << "Last sync: " << std::endl;
    clock_t t1 = clock();
    G.syncInT(diningPhilPair[n - 1]);
    clock_t t2 = clock();
    std::cout<<"Last Sync Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
    //set all event labels in transitions to zero
    for (auto &t : G.transitions){
        t[1] = 0;
    }
    G.lmd[0] = 1;
    //G.transReduceDSV();
}

void dpGroupTest2(size_t n){
    std::vector<System> diningPhils, diningPhilPair;
    size_t cnt = 2;
    System::size_t_in_set tausf = {}, tauss = {};
    System Gf, Gs;
    //create n dining phils
    diningPhils = createDiningPhilosophers2(n);
    //synch in pair
    for(size_t i = 0; i < 2*n; i += 2){
        diningPhilPair.push_back(diningPhils[i].syncInT(diningPhils[i + 1]));
    }

    //first half
    Gf = diningPhilPair[0];
    //tausf.insert({1, 2*n + 1, 4*n + 1});
    clock_t t1 = clock();
    for (size_t i = 1; i < n; ++i){
        Gf.syncInT(diningPhilPair[i]);
    //    tausf.insert({2*i, 2*i+1, 2*n+2*i, 2*n+2*i+1,4*n + i + 1});
    }
    clock_t t2 = clock();
    std::cout<<"Sub sync runtime: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
    //for (auto &t : Gf.transitions){
    //    if (tausf.find(t[1]) != tausf.end()){
    //        t[1] = 0;
    //    }
    //}
    std::cout <<"\nSync State#: " << Gf.states.size() << "\nSync Trans#: " << Gf.transitions.size() << std::endl;
    //Gf.transReduceDSV2();
    //std::cout <<"\nreduced State#: " << Gf.states.size() << "\nreduced Trans#: " << Gf.transitions.size() << std::endl;

    //second half
    //Gs = diningPhilPair[n-1];
    //tauss.insert({2*n - 1, 4*n - 1, 5*n});
    //t1 = clock();
    //for (size_t i = n - 2; i >= n/2; --i){

    //    Gs.syncInT(diningPhilPair[i]);
        //tauss.insert({2*i + 1, 2*i+2, 2*n+2*i+1, 2*n+2*i+2,4*n + i + 1});
    //}
   // t2 = clock();
    //std::cout<<"Sub sync runtime: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;

    //for (auto &t : Gs.transitions){
    //    if (tauss.find(t[1]) != tauss.end()){
    //        t[1] = 0;
    //    }
    //}
    //Gs.transReduceDSV2();
    //t1 = clock();
    //Gf.syncInT(Gs);
    //t2 = clock();
    //std::cout<<"Big sync runtime: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;

    //std::cout <<"\nSync State#: " << Gf.states.size() << "\nSync Trans#: " << Gf.transitions.size() << std::endl;
    for(auto &iter : Gf.transitions){
        iter[1] = 0;
    }
    //Gf.transReduceDSV2();
    //for (const auto &iter : Gf.transitions){
    //    std::cout << iter[0] << " "<<iter[1] <<" " << iter[2] << std::endl;
    //}
    //tauss.insert({})
}

void dpGroupTest2NonDSV(size_t n){
    std::vector<System> diningPhils, diningPhilPair;
    size_t cnt = 2;
    System::size_t_in_set tausf = {}, tauss = {};
    System Gf, Gs;
    //create n dining phils
    diningPhils = createDiningPhilosophers2(n);
    //synch in pair
    for(size_t i = 0; i < 2*n; i += 2){
        diningPhilPair.push_back(diningPhils[i].syncInT(diningPhils[i + 1]));
    }
    //first half
    Gf = diningPhilPair[0];
    tausf.insert({1, 2*n + 1, 4*n + 1});
    for (size_t i = 1; i < n/2; ++i){
        Gf.syncInT(diningPhilPair[i]);
        tausf.insert({2*i, 2*i+1, 2*n+2*i, 2*n+2*i+1,4*n + i + 1});
    }
    for (auto &t : Gf.transitions){
        if (tausf.find(t[1]) != tausf.end()){
            t[1] = 0;
        }
    }
    std::cout <<"\nSync State#: " << Gf.states.size() << "\nSync Trans#: " << Gf.transitions.size() << std::endl;
    Gf.transReduce();
    std::cout <<"\nreduced State#: " << Gf.states.size() << "\nreduced Trans#: " << Gf.transitions.size() << std::endl;

    //second half
    Gs = diningPhilPair[n-1];
    tauss.insert({2*n - 1, 4*n - 1, 5*n});
    for (size_t i = n - 2; i >= n/2; --i){
        Gs.syncInT(diningPhilPair[i]);
        tauss.insert({2*i + 1, 2*i+2, 2*n+2*i+1, 2*n+2*i+2,4*n + i + 1});
    }
    for (auto &t : Gs.transitions){
        if (tauss.find(t[1]) != tauss.end()){
            t[1] = 0;
        }
    }
    Gs.transReduce();

    Gf.syncInT(Gs);
    std::cout <<"\nSync State#: " << Gf.states.size() << "\nSync Trans#: " << Gf.transitions.size() << std::endl;

    for(auto &iter : Gf.transitions){
        iter[1] = 0;
    }
    Gf.transReduce();
    for (const auto &iter : Gf.transitions){
        std::cout << iter[0] << " "<<iter[1] <<" " << iter[2] << std::endl;
    }
    //tauss.insert({})
}

void dpTest2NonDSV(size_t n){
    std::vector<System> diningPhils, diningPhilPair;
    size_t cnt = 2;
    System::size_t_in_set taus = {};
    System G;
    //create n dining phils
    diningPhils = createDiningPhilosophers2(n);
    //synch in pair
    for(size_t i = 0; i < 2*n; i += 2){
        diningPhilPair.push_back(diningPhils[i].syncInT(diningPhils[i + 1]));
    }
    //sync and abstract
    //first step: sych first two and reduce
    G = diningPhilPair[0].syncInT(diningPhilPair[1]);
    std::cout <<"Sync With G" << cnt ++ <<"\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() << std::endl;
    taus.insert({1, 2, 3, 2*n + 1, 2*n + 2, 2*n + 3, 4*n + 1, 4*n + 2});
    for (auto &t : G.transitions){
        if (taus.find(t[1]) != taus.end()){
            t[1] = 0;
        }
    }
    /*
    std::ofstream outFile("/Users/liangxudong/Desktop/VLTS/G_dsv_test");
    outFile <<"des (" << 0 <<", " <<G.transitions.size() <<", " <<G.states.size() <<")\n";
    for (const auto &e : G.transitions) {
        if((e[1] == 0) && (e[0] != 0) && (e[2] != 0)){
            outFile <<"(" << e[0] << ", "  <<"\""<< "tau" <<"\""<<", " << e[2] <<")"<< "\n";
        }
        else{
            outFile <<"(" << e[0] << ", " << "\"" << e[1] << "\""<<", " << e[2] <<")"<< "\n";
        }
    }
     */
    G.transReduce();
    std::cout << "After Reduction State#: " << G.states.size() << "\nAfter reduction Trans#: " << G.transitions.size() << std::endl;
    std::cout << std::endl;
    //second step: sync and reduce until the last pair
    for (size_t i = 2; i < n - 1; ++i){
        G.syncInT(diningPhilPair[i]);
        std::cout <<"Sync With G" << cnt ++ << "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() << std::endl;
        taus.insert({2*i, 2*i + 1, 2*i + 2*n, 2*i + 1 + 2*n, 4*n + i + 1});
        for (auto &t : G.transitions){
            if (taus.find(t[1]) != taus.end()){
                t[1] = 0;
            }
        }
        //G.lmd[0] = 1;
        clock_t t1 = clock();
        G.transReduce();
        clock_t t2 = clock();
        std::cout << "After Reduction State#: " << G.states.size() << "\nAfter Reduction Trans#: " << G.transitions.size() << std::endl;
        std::cout<<"One reduction run time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s"<<std::endl;
        std::cout << std::endl;
    }
    //third step: sync with the last pair and reduce
    G.syncInT(diningPhilPair[n - 1]);
    std::cout <<"Sync With G" << cnt ++ << "\nSync State#: " << G.states.size() << "\nSync Trans#: " << G.transitions.size() << std::endl;
    //at last all transitions are local and invisible.
    for (auto &t : G.transitions){
        t[1] = 0;
    }
    //G.lmd[0] = 1;
    G.transReduce();
    for (const auto &iter : G.transitions){
        std::cout << iter[0] << " "<< iter[1] << " " << iter[2] << std::endl;
    }
}
#endif //DININGPHILOSOPHERS_DININGTEST2_H
