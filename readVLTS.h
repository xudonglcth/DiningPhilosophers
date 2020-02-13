//
// Created by 梁旭东 on 2019-12-10.
//

#ifndef LABELREF_READVLTS_H
#define LABELREF_READVLTS_H
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include "System.h"
typedef std::vector<size_t > state;
typedef std::vector<std::vector<size_t >> transition;
typedef std::vector<size_t > partition;

std::string read_vlts(const std::string& file_path){
    std::ifstream inFile;
    inFile.open(file_path); //open the input file
    std::stringstream strStream;
    strStream << inFile.rdbuf(); //read the file
    std::string str = strStream.str(); //str holds the content of the file
    return str;
}

std::vector<std::string> StringSplit(const std::string& aut, const std::string& delim){
    std::vector <std::string>  VecString;
    auto start = 0U;
    auto end = aut.find(delim);
    while(end != std::string::npos){
        VecString.push_back(aut.substr(start, end - start));
        start = end + 1;
        end = aut.find(delim, start);
    }
    return VecString;
}

std::tuple<state, transition> BigTest(const std::string& TestCase){
    std::vector <std::string>  VecString, VecStringTemp;
    std::string VLTS_Case = read_vlts("/Users/liangxudong/Desktop/VLTS/" + TestCase + ".aut");
    std::string delim = "\n";
    std::vector<size_t > SourceVec, EventVec, TargetVec;
    std::vector<std::string> EventVecRaw;
    VecString = StringSplit(VLTS_Case, delim);
    for (size_t i = 1; i < VecString.size(); ++i){
        VecString[i].erase(0, 1);
        VecString[i].erase(VecString[i].length() - 1, 1);
        VecString[i] += ", ";
        VecStringTemp = StringSplit(VecString[i], ", ");
        SourceVec.push_back(std::stoi(VecStringTemp[0]));
        TargetVec.push_back(std::stoi(VecStringTemp.back().erase(0, 1)));
        EventVecRaw.push_back(VecStringTemp[1]);
    }

    std::string description = VecString[0] + ", ";
    std::vector<std::string> VecDes = StringSplit(description, ", ");
    std::string statesNumber = VecDes.back();
    statesNumber.erase(0, 1);
    statesNumber.erase(statesNumber.length() - 1, 1);
    size_t n = std::stoi(statesNumber);
    std::map<std::string, int > EventMap;
    EventMap[" \"tau\""] = 0;
    int count = 0;
    for(const auto& j : EventVecRaw){
        if (EventMap.find(j) == EventMap.end()){
            count++;
            EventMap[j] = count;
        }
        EventVec.push_back(EventMap[j]);
    }
    transition VLTS_Trans;
    std::vector<size_t > subTrans;
    for(size_t i = 0; i < EventVec.size(); ++i){
        subTrans.push_back(SourceVec[i]);
        subTrans.push_back(EventVec[i]);
        subTrans.push_back(TargetVec[i]);
        VLTS_Trans.push_back(subTrans);
        subTrans.clear();
    }
    std::vector<size_t > VLTS_State;
    VLTS_State.reserve(n);
    for(size_t i = 0; i < n; i++){
        VLTS_State.push_back(i);
    }
    std::tuple<state, transition> res = {VLTS_State, VLTS_Trans};
    return res;
    clock_t t1 = clock();
    clock_t t2 = clock();
    std::cout<<"Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
}
void testPerformance(const std::string& testCase){
    System t;
    t.states = std::get<0> (BigTest(testCase));
    t.transitions = std::get<1> (BigTest(testCase));
    t.tau = {0};
    clock_t t1 = clock();
    t.transReduce();
    clock_t t2 = clock();
    std::cout<<"Sigref: "<<std::endl;
    std::cout<<"After Reduction State#: "<<*std::max_element(t.partitions_new.begin(), t.partitions_new.end()) + 1<< "\nAfter Reduction Trans#: "<< t.transitions.size()<<std::endl;
    std::cout<<"Run Time: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<"s\n"<<std::endl;
}
#endif //LABELREF_READVLTS_H
