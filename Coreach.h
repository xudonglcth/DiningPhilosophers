//
// Created by 梁旭东 on 2020-03-18.
//

#ifndef DININGPHILOSOPHERS_COREACH_H
#define DININGPHILOSOPHERS_COREACH_H

#include "System.h"
#include "stack"
std::vector<bool> System::coreach(){
    std::vector<bool> visited(states.size(), false);
    std::stack<size_t > toConsider;
    boundaries();
    for (const auto &iter : init){
        toConsider.push(iter);
    }
    while(!toConsider.empty()){
        std::size_t currentState = toConsider.top();
        visited[currentState] = true;
        toConsider.pop();
        for (size_t i = lower_bound[currentState]; i < upper_bound[currentState]; ++i){
            if (!visited[transInStateOrder[i][0]]){
                visited[transInStateOrder[i][0]] = true;
                toConsider.push(transInStateOrder[i][0]);
            }
        }

    }
    for (const auto &i : visited){
        std::cout << i << " ";
    }
    std::cout << std::endl;
    return visited;
}
#endif //DININGPHILOSOPHERS_COREACH_H
