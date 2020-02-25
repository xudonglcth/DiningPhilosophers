//
// Created by 梁旭东 on 2020-02-11.
//

#ifndef DININGPHILOSOPHERS_SYSTEM_H
#define DININGPHILOSOPHERS_SYSTEM_H
#include <set>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
class System{
public:
    using size_t_in_vec = std::vector<size_t >;
    using size_t_in_vec_in_vec =  std::vector<std::vector<size_t >>;
    using size_t_in_set = std::set<size_t >;
    using hashtable_t = std::map<std::pair<size_t, std::set<std::pair<size_t , size_t >>>, size_t >;
    using delta_t = std::map<std::vector<size_t >, std::set<size_t >>;
    using synced_delta_t = std::map<std::vector<size_t>, std::set<std::vector<size_t > > >;
    using sigmaFeasible_t = std::map <size_t, std::set<size_t>>;
    using synced_states_map_t = std::map<std::vector<size_t >, size_t >;
    using label_t = std::vector<std::set<std::pair<size_t, size_t>>>;
    using size_t_in_vec_in_set = std::set<std::vector<size_t> >;
    using synced_sigmaFeasible_t = std::map <std::vector<size_t>, std::set<size_t>>;
    size_t_in_vec states, partitions, partitions_new, lmd, lower_bound, upper_bound;
    size_t_in_vec_in_vec transitions, transPerState;
    size_t_in_set tau = {0}, sigma, init, sigma_sb;
    size_t_in_vec_in_set X;
    delta_t delta;
    synced_delta_t synced_delta;
    sigmaFeasible_t sigmaFeasible;
    synced_sigmaFeasible_t synced_sigmaFeasible;
    label_t label;
    std::map<size_t_in_vec, size_t> newStates;
    //Constructors;
    System() = default;

    System(size_t_in_vec s, size_t_in_vec_in_vec t, size_t_in_set i)
            :
    states(s), transitions(t), init(i)
    {}

    System(size_t_in_vec s, size_t_in_vec_in_vec t, size_t_in_vec p, size_t_in_set tau_)
            :
            states(std::move(s)), transitions(std::move(t)), lmd(std::move(p)), tau(std::move(tau_))
    {}
    //Some internal functions
    void transInfo(){
        delta.clear();
        synced_delta.clear();
        sigmaFeasible.clear();
        sigma.clear();
        X.clear();sigma_sb.clear();synced_sigmaFeasible.clear();
        for (const auto &i : transitions){
            //trans to delta
            if (delta.find({i[0], i[1]}) == delta.end()){
                delta[{i[0], i[1]}] = {i[2]};
            }
            else{
                delta[{i[0], i[1]}].insert(i[2]);
            }
            //build sigma_f
            if (sigmaFeasible.find(i[0]) == sigmaFeasible.end()){
                sigmaFeasible[i[0]] = {i[1]};
            }
            else{
                sigmaFeasible[i[0]].insert(i[1]);
            }
            //build sigma
            sigma.insert(i[1]);
            //default init
            init = {0};
        }
    }
    //Sync
    System& sync(System& G_toSyncWith){
        transInfo(); G_toSyncWith.transInfo();
        size_t_in_set sigma_il = sigma, set1, set2, set3, res, sb;
        size_t_in_vec_in_set delta_val;
        size_t_in_vec delta_key;
        size_t k = 0;
        delta_t delta_temp = G_toSyncWith.delta;
        sigmaFeasible_t sf_temp = G_toSyncWith.sigmaFeasible;
        std::set_intersection(sigma.begin(), sigma.end(),
                              G_toSyncWith.sigma.begin(), G_toSyncWith.sigma.end(), std::inserter(sb, sb.begin()));
        sigma_sb = sb;
        std::vector<std::set<std::vector<size_t > > > delta_X;
        for (auto ele : G_toSyncWith.sigma){
            sigma_il.insert(ele);
        }
        for (auto ele : sigma_sb){
            sigma_il.erase(ele);
        }
        setCrossProduct(init, G_toSyncWith.init, X);
        delta_X.push_back(X);

        do{
            ++k;
            delta_X.emplace_back();
            for (const auto& i : delta_X[k - 1]){
                set1.clear(); set2.clear(); set3.clear();
                synced_sigmaFeasible[i] = {};

                //first for loop
                res.clear();

                std::set_intersection(sigmaFeasible[i[0]].begin(), sigmaFeasible[i[0]].end(),
                                      G_toSyncWith.sigmaFeasible[i[1]].begin(), G_toSyncWith.sigmaFeasible[i[1]].end(),
                                      std::inserter(set1, set1.begin()));
                std::set_intersection(set1.begin(), set1.end(),
                                      sigma_sb.begin(), sigma_sb.end(),
                                      std::inserter(res, res.begin()));
                set1 = res;
                for (auto a : set1){
                    synced_sigmaFeasible[i].insert(a);
                    delta_key = {i[0], i[1], a};
                    delta_val.clear();
                    setCrossProduct(delta[{i[0], a}], G_toSyncWith.delta[{i[1], a}], delta_val);
                    if (synced_delta[delta_key].empty()){
                        synced_delta[delta_key] = delta_val;
                    }
                    else{
                        for (const auto& ele : delta_val){
                            synced_delta[delta_key].insert(ele);
                        }
                    }
                    for (const auto& val:delta_val){
                        delta_X[k].insert(val);
                    }
                }

                //second for loop
                res.clear();
                std::set_intersection(sigmaFeasible[i[0]].begin(), sigmaFeasible[i[0]].end(),
                                      sigma_il.begin(), sigma_il.end(),
                                      std::inserter(res, res.end()));
                set2 = res;
                for (auto a : set2){
                    synced_sigmaFeasible[i].insert(a);
                    delta_key = {i[0], i[1], a};
                    delta_val.clear();
                    setCrossProduct(delta[{i[0], a}], {i[1]}, delta_val);
                    if (synced_delta[delta_key].empty()){
                        synced_delta[delta_key] = delta_val;
                    }
                    else{
                        for (const auto& ele : delta_val){
                            synced_delta[delta_key].insert(ele);
                        }
                    }
                    for (const auto& val:delta_val){
                        delta_X[k].insert(val);
                    }
                }

                //third for loop
                res.clear();
                std::set_intersection(G_toSyncWith.sigmaFeasible[i[1]].begin(), G_toSyncWith.sigmaFeasible[i[1]].end(),
                                      sigma_il.begin(), sigma_il.end(),
                                      std::inserter(res, res.end()));
                set3 = res;
                for (auto a : set3){
                    synced_sigmaFeasible[i].insert(a);
                    delta_key = {i[0], i[1], a};
                    delta_val.clear();
                    setCrossProduct({i[0]}, G_toSyncWith.delta[{i[1], a}], delta_val);
                    if (synced_delta[delta_key].empty()){
                        synced_delta[delta_key] = delta_val;
                    }
                    else{
                        for (const auto& ele : delta_val){
                            synced_delta[delta_key].insert(ele);
                        }
                    }
                    for (const auto& val:delta_val){
                        delta_X[k].insert(val);
                    }
                }
            }
            for (const auto& ele: X){
                delta_X[k].erase(ele);
            }
            for (const auto& ele : delta_X[k]){
                X.insert(ele);
            }
        }while(! delta_X[k].empty());

        syncResult(G_toSyncWith);
        return *this;
    }
    static void setCrossProduct(const size_t_in_set& s1, const size_t_in_set& s2, std::set<std::vector<size_t> > &s){
        for (auto i : s1){
            for (auto j : s2){
                s.insert({i, j});
            }
        }
    }
    void syncResult(System G_toSyncWith){
        std::map<size_t_in_vec, size_t> state_map;
        size_t_in_vec sync_lmd;
        states.clear();
        transitions.clear();
        size_t cnt = 0;
        for (const auto &synced_state : X){
            state_map[synced_state] = cnt;
            states.push_back(cnt++);
            sync_lmd.push_back(std::min(lmd[synced_state[0]], G_toSyncWith.lmd[synced_state[1]]));
        }
        lmd = sync_lmd;
        for (const auto &synced_state : X){
            for (const auto &feasible_event : synced_sigmaFeasible[synced_state]){
                for (const auto &synced_target : synced_delta[{synced_state[0], synced_state[1], feasible_event}]){
                    transitions.push_back({state_map[synced_state], feasible_event, state_map[synced_target]});
                }
            }
        }
    }

    System& syncDelta(System& G_toSyncWith){
        size_t_in_set sigma_il = sigma, set1, set2, set3, res, sb;
        size_t_in_vec_in_set delta_val;
        size_t_in_vec delta_key;
        size_t k = 0;
        delta_t delta_temp = G_toSyncWith.delta;
        sigmaFeasible_t sf_temp = G_toSyncWith.sigmaFeasible;
        std::set_intersection(sigma.begin(), sigma.end(),
                              G_toSyncWith.sigma.begin(), G_toSyncWith.sigma.end(), std::inserter(sb, sb.begin()));
        sigma_sb = sb;
        std::vector<std::set<std::vector<size_t > > > delta_X;
        for (auto ele : G_toSyncWith.sigma){
            sigma_il.insert(ele);
        }
        for (auto ele : sigma_sb){
            sigma_il.erase(ele);
        }
        setCrossProduct(init, G_toSyncWith.init, X);
        delta_X.push_back(X);

        do{
            ++k;
            delta_X.emplace_back();
            for (const auto& i : delta_X[k - 1]){
                set1.clear(); set2.clear(); set3.clear();
                synced_sigmaFeasible[i] = {};

                //first for loop
                res.clear();

                std::set_intersection(sigmaFeasible[i[0]].begin(), sigmaFeasible[i[0]].end(),
                                      G_toSyncWith.sigmaFeasible[i[1]].begin(), G_toSyncWith.sigmaFeasible[i[1]].end(),
                                      std::inserter(set1, set1.begin()));
                std::set_intersection(set1.begin(), set1.end(),
                                      sigma_sb.begin(), sigma_sb.end(),
                                      std::inserter(res, res.begin()));
                set1 = res;
                for (auto a : set1){
                    synced_sigmaFeasible[i].insert(a);
                    delta_key = {i[0], i[1], a};
                    delta_val.clear();
                    setCrossProduct(delta[{i[0], a}], G_toSyncWith.delta[{i[1], a}], delta_val);
                    if (synced_delta[delta_key].empty()){
                        synced_delta[delta_key] = delta_val;
                    }
                    else{
                        for (const auto& ele : delta_val){
                            synced_delta[delta_key].insert(ele);
                        }
                    }
                    for (const auto& val:delta_val){
                        delta_X[k].insert(val);
                    }
                }

                //second for loop
                res.clear();
                std::set_intersection(sigmaFeasible[i[0]].begin(), sigmaFeasible[i[0]].end(),
                                      sigma_il.begin(), sigma_il.end(),
                                      std::inserter(res, res.end()));
                set2 = res;
                for (auto a : set2){
                    synced_sigmaFeasible[i].insert(a);
                    delta_key = {i[0], i[1], a};
                    delta_val.clear();
                    setCrossProduct(delta[{i[0], a}], {i[1]}, delta_val);
                    if (synced_delta[delta_key].empty()){
                        synced_delta[delta_key] = delta_val;
                    }
                    else{
                        for (const auto& ele : delta_val){
                            synced_delta[delta_key].insert(ele);
                        }
                    }
                    for (const auto& val:delta_val){
                        delta_X[k].insert(val);
                    }
                }

                //third for loop
                res.clear();
                std::set_intersection(G_toSyncWith.sigmaFeasible[i[1]].begin(), G_toSyncWith.sigmaFeasible[i[1]].end(),
                                      sigma_il.begin(), sigma_il.end(),
                                      std::inserter(res, res.end()));
                set3 = res;
                for (auto a : set3){
                    synced_sigmaFeasible[i].insert(a);
                    delta_key = {i[0], i[1], a};
                    delta_val.clear();
                    setCrossProduct({i[0]}, G_toSyncWith.delta[{i[1], a}], delta_val);
                    if (synced_delta[delta_key].empty()){
                        synced_delta[delta_key] = delta_val;
                    }
                    else{
                        for (const auto& ele : delta_val){
                            synced_delta[delta_key].insert(ele);
                        }
                    }
                    for (const auto& val:delta_val){
                        delta_X[k].insert(val);
                    }
                }
            }
            for (const auto& ele: X){
                delta_X[k].erase(ele);
            }
            for (const auto& ele : delta_X[k]){
                X.insert(ele);
            }
        }while(! delta_X[k].empty());
        //result in delta
        std::map<size_t_in_vec, size_t> states_map;
        size_t cnt = 0; size_t_in_vec s_lmd;
        states.clear(); delta.clear(); sigmaFeasible.clear();sigma.clear();
        sigma_sb.clear();synced_sigmaFeasible.clear();
        for (const auto &x : X){
            states_map[x] = cnt;
            states.push_back(cnt++);
            s_lmd.push_back(std::min(lmd[x[0]], G_toSyncWith.lmd[x[1]]));
        }
        lmd = s_lmd;
        for (const auto &ele : synced_delta){
            size_t s = states_map[{ele.first[0], ele.first[1]}];
            size_t t = ele.first[2];
            sigma.insert(t);
            sigmaFeasible[s].insert(ele.first[2]);
            for (const auto &target : ele.second){
                delta[{s, t}].insert(states_map[target]);
            }
        }
        X.clear(); synced_delta.clear();
        return *this;
    }
    //Abstraction
    void boundaries(){
        lower_bound.clear();
        upper_bound.clear();
        for(size_t i = 0; i < states.size(); i++){
            lower_bound.push_back(0);
            upper_bound.push_back(0);
        }
        for(size_t j = 0; j < transitions.size(); j++){
            transPerState.push_back({0, 0, 0});
        }
        std::vector<size_t > indices(states.size() + 1, 0);

        for(const auto &i : transitions){
            ++indices[i[2]];
        }
        size_t sum = 0;
        for (auto &j : indices){
            sum += j;
            j = sum;
        }
        for(const auto &k : transitions){
            indices[k[2]]--;
            transPerState[indices[k[2]]] = {k[0], k[1], k[2]};
        }
        for(size_t i = 0; i < lower_bound.size(); i++){
            lower_bound[i] = indices[i];
            upper_bound[i] = indices[i + 1];
        }
    }
    void labelInsert(size_t target, size_t event, size_t block){
        //size_t size_before = label[target].size(), size_after;
        //std::vector<std::pair<size_t, size_t> > to_insert, vec_res, vec_temp(label[target].begin(), label[target].end());
        //to_insert.push_back({std::make_pair(event, block)});
        //std::set_union(vec_temp.begin(), vec_temp.end(), to_insert.begin(), to_insert.end(),
        //       std::back_inserter(vec_res));
        //std::set<std::pair<size_t, size_t> > set_temp(vec_res.begin(), vec_res.end());
        //label[target] = set_temp;
        //size_after = label[target].size();
        //if(size_before != size_after)
        if (label[target].insert(std::make_pair(event, block)).second)
        {
            for(size_t i = lower_bound[target]; i < upper_bound[target]; i++){
                //in_going_trans ingoingTrans = inTransPerState[i];
                if (tau.find(transPerState[i][1]) != tau.end()
                    && partitions[target] == partitions[transPerState[i][0]]){
                    labelInsert(transPerState[i][0], event, block);
                }
            }
        }
    }
    void transReduce(){
        std::set<size_t > sts;
        //std::cout << "Max capacity: "<<sts.max_size()<<std::endl;
        size_t n = states.size(), count;
        hashtable_t hashtable;
        boundaries();
        partitions_new = lmd;
        if (partitions_new.empty()){
            for(size_t i = 0; i < n; i++){
                partitions_new.push_back(0);
            }
        }

        do{
            count = 0;
            partitions = partitions_new;
            label.clear();
            hashtable.clear();
            for(size_t i = 0; i < n; i++){
                label.emplace_back();
            }

            for(auto current_trans : transitions){
                size_t source = current_trans[0], event = current_trans[1], target = current_trans[2];
                if (!(tau.find(event) != tau.end() && partitions[source] == partitions[target])){
                    labelInsert(source, event, partitions[target]);
                }
            }

            for (auto current_state : states){
                if(hashtable.find(std::make_pair(partitions[current_state], label[current_state])) == hashtable.end()){
                    hashtable[std::make_pair(partitions[current_state], label[current_state])] = count++;
                }
            }

            for (auto current_state : states){
                partitions_new[current_state] = hashtable[std::make_pair(partitions[current_state], label[current_state])];
            }
        }while(partitions_new != partitions);
        size_t_in_vec_in_set sPi;
        for (const auto &iter : transitions) {
            if (!(partitions[iter[0]] == partitions[iter[2]] && (tau.find(iter[1]) != tau.end()))) {
                sPi.insert({partitions[iter[0]], iter[1], partitions[iter[2]]});
            }
        }
        transitions.clear();
        states.clear();
        lmd.clear();
        for (const auto &iter : sPi){
            transitions.push_back(iter);
        }
        size_t nPi = *std::max_element(partitions_new.begin(), partitions_new.end()) + 1;
        for (size_t i =  0; i < nPi; ++i){
            states.push_back(i);
            lmd.push_back(0);
        }
        //for (const auto& i : partitions_new){
        //    std::cout << i << " ";
        //}
        //std::cout << std::endl;
    }

    System & syncInT(System &g_to_sync);

    void outgoingBoundaries();

    void transSigma();
};
#endif //DININGPHILOSOPHERS_SYSTEM_H
//
//  System.h
//  system
//
//  Created by 梁旭东 on 2020-02-12.
//  Copyright © 2020 xudongl. All rights reserved.
//

#ifndef System_h
#define System_h


#endif /* System_h */
