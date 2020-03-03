//
// Created by 梁旭东 on 2020-02-24.
//

#ifndef DININGPHILOSOPHERS_SYNCINT_H
#define DININGPHILOSOPHERS_SYNCINT_H

#include "DiningPhilosophersNew.h"
System& System::syncInT(System& g_to_sync){
    /*
     * Note:Extra sb(shared events) are ignored.
     */
    size_t non_deterministic = 0;
    std::map<std::vector<size_t>, size_t > synced_state_map;
    std::vector<std::set<std::vector<size_t > > > delta_x;
    size_t_in_vec v1, v2, v_temp, v_sb, v_union, v_diff;
    size_t k = 0, x1, x2, cnt = 0;
    init = {0};
    transSigma();
    g_to_sync.transSigma();
    outgoingBoundaries();
    g_to_sync.outgoingBoundaries();
    size_t_in_vec vt1(sigma.begin(), sigma.end()), vt2(g_to_sync.sigma.begin(), g_to_sync.sigma.end());
    std::set_intersection(vt1.begin(),vt1.end(), vt2.begin(), vt2.end(), std::back_inserter(v_sb));
    std::set_union(vt1.begin(), vt1.end(), vt2.begin(), vt2.end(), std::back_inserter(v_union));
    std::set_difference(v_union.begin(), v_union.end(), v_sb.begin(), v_sb.end(), std::back_inserter(v_diff));
    transitions.clear();
    states.clear();
    X.clear();
    setCrossProduct(init, g_to_sync.init, X);
    delta_x.push_back(X);
    do{
        ++k;
        delta_x.emplace_back();
        for (const auto &state_pair : delta_x[k - 1]){
            v1.clear(); v2.clear(); v_temp.clear();
            x1 = state_pair[0]; x2 = state_pair[1];
            for (size_t i = lower_bound[x1]; i < upper_bound[x1]; ++i){
                v1.push_back(transInStateOrder[i][1]);
            }
            for (size_t i = g_to_sync.lower_bound[x2]; i < g_to_sync.upper_bound[x2]; ++i){
                v2.push_back(g_to_sync.transInStateOrder[i][1]);
            }
            std::sort(v1.begin(), v1.end()); std::sort(v2.begin(), v2.end());
            //first for loop:
            std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(v_temp));
            for (const auto &a : v_temp){
                size_t_in_vec_in_vec target;
                /*
                    for (size_t i = lower_bound[x1]; i < upper_bound[x1]; ++i){
                        size_t_in_vec target_instance(2, 0);
                        if (transInStateOrder[i][1] == a){
                            target_instance[0] = transInStateOrder[i][2];
                            for (size_t i = g_to_sync.lower_bound[x2]; i < g_to_sync.upper_bound[x2]; ++i){
                                if (g_to_sync.transInStateOrder[i][1] == a){
                                    target_instance[1] = g_to_sync.transInStateOrder[i][2];
                                    target.push_back(target_instance);
                                }
                            }
                        }
                    }
                */
                findTarget(g_to_sync, a, x1, x2, v_temp, target);
                if (synced_state_map.find({x1, x2}) == synced_state_map.end()){
                    synced_state_map[{x1, x2}] = cnt;
                    states.push_back(cnt++);
                }
                for (const auto &target_iter : target){
                    if (synced_state_map.find(target_iter) == synced_state_map.end()){
                        synced_state_map[target_iter] = cnt;
                        states.push_back(cnt++);
                    }
                    transitions.push_back({synced_state_map[{x1, x2}], a, synced_state_map[target_iter]});
                    delta_x[k].insert(target_iter);
                }
            }
            //second for loop:
            v_temp.clear();
            std::set_intersection(v1.begin(), v1.end(), v_diff.begin(), v_diff.end(), std::back_inserter(v_temp));
            for (const auto &a : v_temp){
                size_t_in_vec target_instance = {0, x2};
                size_t_in_vec_in_vec target;
                for (size_t i = lower_bound[x1]; i < upper_bound[x1]; ++i){
                    if (transInStateOrder[i][1] == a){
                        target_instance[0] = transInStateOrder[i][2];
                        target.push_back(target_instance);
                    }
                }
                if (synced_state_map.find({x1, x2}) == synced_state_map.end()){
                    synced_state_map[{x1, x2}] = cnt;
                    states.push_back(cnt++);
                }
                for (const auto &target_iter : target){
                    if (synced_state_map.find(target_iter) == synced_state_map.end()){
                        synced_state_map[target_iter] = cnt;
                        states.push_back(cnt++);
                    }
                    transitions.push_back({synced_state_map[{x1, x2}], a, synced_state_map[target_iter]});
                    delta_x[k].insert(target_iter);
                }
            }
            //third for loop
            v_temp.clear();
            std::set_intersection(v2.begin(), v2.end(), v_diff.begin(), v_diff.end(), std::back_inserter(v_temp));
            for (const auto &a : v_temp){
                size_t_in_vec target_instance = {x1, 0};
                size_t_in_vec_in_vec target;
                for (size_t i = g_to_sync.lower_bound[x2]; i < g_to_sync.upper_bound[x2]; ++i){
                    if (g_to_sync.transInStateOrder[i][1] == a){
                        target_instance[1] = g_to_sync.transInStateOrder[i][2];
                        target.push_back(target_instance);
                    }
                }
                if (synced_state_map.find({x1, x2}) == synced_state_map.end()){
                    synced_state_map[{x1, x2}] = cnt;
                    states.push_back(cnt++);
                }
                for (const auto &target_iter : target){
                    if (synced_state_map.find(target_iter) == synced_state_map.end()){
                        synced_state_map[target_iter] = cnt;
                        states.push_back(cnt++);
                    }
                    transitions.push_back({synced_state_map[{x1, x2}], a, synced_state_map[target_iter]});
                    delta_x[k].insert(target_iter);
                }
            }
        }
        assert((synced_state_map[{0, 0}]) == 0);
        for (const auto& ele: X){
            delta_x[k].erase(ele);
        }
        for (const auto& ele : delta_x[k]){
            X.insert(ele);
        }
    }while(!delta_x[k].empty());
    lmd.clear();
    for (size_t i = 0; i < states.size(); ++i){
        lmd.push_back(0);
    }
    lmd[0] = 1;
    /*
    for(const auto &iter : transitions){
        std::cout << iter[0] << " " << iter[1] << " " << iter[2]<<std::endl;
    }

    for(const auto &iter : states){
        std::cout << iter << " ";
    }*/
    //std::cout << "\n"<<std::endl;
    return *this;
}

void System::outgoingBoundaries(){
    transInStateOrder.clear();
    lower_bound.clear();
    upper_bound.clear();
    size_t max_state = *std::max_element(states.begin(), states.end());
    for(size_t i = 0; i < max_state + 1; i++){
        lower_bound.push_back(0);
        upper_bound.push_back(0);
    }
    for(size_t j = 0; j < transitions.size(); j++){
        transInStateOrder.push_back({0, 0, 0});
    }
    std::vector<size_t > indices(max_state + 2, 0);
    for(const auto &i : transitions){
        ++indices[i[0]];
    }
    size_t sum = 0;
    for (auto &j : indices){
        sum += j;
        j = sum;
    }
    for(const auto &k : transitions){
        transInStateOrder[--indices[k[0]]] = {k[0], k[1], k[2]};
    }
    for(size_t i = 0; i < lower_bound.size(); ++i){
        lower_bound[i] = indices[i];
        upper_bound[i] = indices[i + 1];
    }
}


void System::findTarget(const System &g_to_sync, const size_t a, const size_t x1, const size_t x2, const size_t_in_vec &v_temp, size_t_in_vec_in_vec& target){
    for (size_t i = lower_bound[x1]; i < upper_bound[x1]; ++i){
        size_t_in_vec target_instance(2, 0);
        if (transInStateOrder[i][1] == a){
            target_instance[0] = transInStateOrder[i][2];
            for (size_t i = g_to_sync.lower_bound[x2]; i < g_to_sync.upper_bound[x2]; ++i){
                if (g_to_sync.transInStateOrder[i][1] == a){
                    target_instance[1] = g_to_sync.transInStateOrder[i][2];
                    target.push_back(target_instance);
                }
            }
        }
    }
}


void System::transSigma(){
    sigma.clear();
    for (const auto &t : transitions){
        sigma.insert(t[1]);
    }
}

#endif //DININGPHILOSOPHERS_SYNCINT_H
