//
// Created by 梁旭东 on 2020-03-17.
//

#ifndef DININGPHILOSOPHERS_TWOPHILS_H
#define DININGPHILOSOPHERS_TWOPHILS_H

#include "DiningPhilosophersNew.h"
void twoPhils(){
    std::vector<System> twoP;
    twoP = createDiningPhilosophersNew(2);
    for (size_t i = 1; i < 4; ++i){
        twoP[0].syncInT(twoP[i]);
    }
    for (auto & iter : twoP[0].transitions){
        iter[1] = 0;
    }
    twoP[0].lmd[0] = 1;
    twoP[0].transReduceDSV();
    for (const auto & iter : twoP[0].transitions){
        std::cout << iter[0] <<" "<< iter[1] <<" "<< iter[2] << std::endl;
    }
}
#endif //DININGPHILOSOPHERS_TWOPHILS_H
