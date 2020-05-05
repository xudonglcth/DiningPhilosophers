//
// Created by 梁旭东 on 2020-04-07.
//

#ifndef DININGPHILOSOPHERS_REACHBOOLTEST_H
#define DININGPHILOSOPHERS_REACHBOOLTEST_H

#include "System.h"
std::tuple<std::vector<size_t>, std::vector<std::vector<size_t> >, std::set<size_t> > reachTest1(size_t n){
    std::vector<size_t > ndelta = {2, 1, 1, 2, 1, 1, 4, 0, 0, 0};
    std::vector<std::vector<size_t> > delta = {{2, 3, 4, 5, 6, 7, 2, 0, 0, 0}, {3, 0, 0, 7, 0, 0, 8, 0, 0, 0},
                                               {0, 0, 0, 0, 0, 0, 9, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 10, 0, 0, 0}};
    size_t q = 0;
    std::vector<size_t > init(n, 0);
    std::vector<size_t> ndeltaQ(10 * n, 0);
    std::vector<size_t> zeroV(10*n, 0);
    std::vector<std::vector<size_t > > deltaQ({zeroV, zeroV, zeroV, zeroV});
    for (size_t i = 1; i < n+1; ++i){
        for (size_t j = 0; j < 7; ++j){
            ndeltaQ[10*(i - 1) + j] = ndelta[j];
        }
        for (size_t j = 0; j < 7; ++j){
            deltaQ[0][10*(i - 1) + j] = delta[0][j] + q;
        }
        deltaQ[1][10*(i - 1) + 0] = delta[1][0] + q;
        deltaQ[1][10*(i - 1) + 3] = delta[1][3] + q;
        deltaQ[1][10*(i - 1) + 6] = delta[1][6] + q;
        deltaQ[2][10*(i - 1) + 6] = delta[2][6] + q;
        deltaQ[3][10*(i - 1) + 6] = delta[3][6] + q;
        init[i - 1] = 1 + 10 * (i - 1);
        q += 10;
    }
    std::set<size_t > initQ(init.begin(), init.end());
    return {ndeltaQ, deltaQ, initQ};
}

std::tuple<std::vector<size_t>, std::vector<std::vector<size_t> >> reachTest2(size_t n){
    size_t q = 0;
    std::vector<size_t> ndeltaQ(n * n, 0);
    std::vector<size_t> zeroV(n*n, 0);
    std::vector<std::vector<size_t > > deltaQ({zeroV, zeroV});
    for (size_t ir = 0; ir < n - 1; ++ir){
        for (size_t ic = 0; ic < n - 1; ++ic){
            ndeltaQ[n * ir + ic] = 2;
            deltaQ[0][n * ir + ic] = n*ir +ic + 2;
            deltaQ[1][n * ir + ic] = n*(ir + 1) + ic + 1;
        }
        for (size_t i = 0; i < n - 1; ++i){
            ndeltaQ[n * (n - 1) + i] = 1;
            deltaQ[0][n * (n - 1) + i] = n * (n - 1) + i + 2;
            ndeltaQ[n * i + n - 1] = 1;
            deltaQ[0][n * i + n - 1] = n * (i + 1) + n;
        }
    }
    return {ndeltaQ, deltaQ};
}


#endif //DININGPHILOSOPHERS_REACHBOOLTEST_H
