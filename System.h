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
#include <unordered_map>
#include <stack>
#include <algorithm>
#include <iterator>
#include <boost/functional/hash.hpp>
class System{
public:
    using size_t_in_vec = std::vector<size_t >;
    using size_t_in_vec_in_vec =  std::vector<std::vector<size_t >>;
    using size_t_in_set = std::set<size_t >;
    using hashtable_t = std::unordered_map<std::string, size_t >;
    using delta_t = std::map<std::vector<size_t >, std::set<size_t >>;
    using synced_delta_t = std::map<std::vector<size_t>, std::set<std::vector<size_t > > >;
    using sigmaFeasible_t = std::map <size_t, std::set<size_t>>;
    using synced_states_map_t = std::map<std::vector<size_t >, size_t >;
    using label_t = std::vector<std::set<std::pair<size_t, size_t>>>;
    using size_t_in_vec_in_set = std::set<std::vector<size_t> >;
    using synced_sigmaFeasible_t = std::map <std::vector<size_t>, std::set<size_t>>;
    using label_each_t = std::set<std::pair<size_t, size_t>>;
    std::vector<std::vector<size_t> > deltaR;
    std::vector<size_t> ndeltaR;

    size_t_in_vec states, partitions, partitions_new, lmd, lower_bound, upper_bound;
    size_t_in_vec_in_vec transitions, transInStateOrder;
    size_t_in_set tau = {0}, sigma, init, sigma_sb;
    size_t_in_vec_in_set X;
    delta_t delta, deltaC;
    synced_delta_t synced_delta;
    sigmaFeasible_t sigmaFeasible, sigmaFeasibleC;
    synced_sigmaFeasible_t synced_sigmaFeasible;
    label_t label;
    std::map<size_t_in_vec, size_t> newStates;
    std::vector<bool> divergent;
    template <typename Container> // we can make this generic for any container [1]
    struct container_hash {
        std::size_t operator()(Container const& c) const {
            return boost::hash_range(c.begin(), c.end());
        }
    };
    std::unordered_map<label_each_t, size_t, container_hash<label_each_t>> label_map;
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

    ///a constructor for coreachable states computation
    System(delta_t dc, sigmaFeasible_t sc)
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
    void tauSccs(){
        for (size_t i = 0; i < states.size(); ++i){
            divergent.push_back(false);
        }
        size_t unused = 1, lastScc = 1;
        size_t_in_vec scc(states.size(), 0);
        size_t_in_vec lowValue(states.size(), 0);
        std::stack<size_t > stack, sccStack;
        outgoingBoundaries();
        for (size_t i = 0; i < states.size(); ++i){
            if (scc[i] == 0){
                stack.push(i);
            }
            while(!stack.empty()){
                const size_t current_state = stack.top();

                if(lowValue[current_state] == 0 && scc[current_state] == 0){
                    scc[current_state] = unused;
                    lowValue[current_state] = unused++;
                    sccStack.push(current_state);

                    for (size_t j = lower_bound[current_state]; j < upper_bound[current_state]; ++j){
                        if(lowValue[transInStateOrder[j][2]] == 0
                            && scc[transInStateOrder[j][2]] == 0
                            && transInStateOrder[j][1] == 0){
                            std::cout << transInStateOrder[j][2] << " In" << std::endl;
                            stack.push(transInStateOrder[j][2]);
                        }
                    }
                }
                else{
                    for (size_t j = lower_bound[current_state]; j < upper_bound[current_state]; ++j){
                        if (lowValue[transInStateOrder[j][2]] != 0 && transInStateOrder[j][1] == 0){
                            lowValue[current_state] =
                                    lowValue[current_state] < lowValue[transInStateOrder[j][2]] ?
                                    lowValue[current_state] : lowValue[transInStateOrder[j][2]];
                        }
                    }
                    if (lowValue[current_state] == scc[current_state]){
                        std::size_t tos, sccId = lastScc++;
                        size_t_in_vec thisScc;
                        do{
                            tos = sccStack.top();
                            lowValue[tos] = 0;
                            scc[tos] = sccId;
                            sccStack.pop();
                            thisScc.push_back(tos);
                        }while(tos != current_state);

                        //non-trivial
                        if (thisScc.size() == 1){
                            for (size_t j = lower_bound[current_state]; j < upper_bound[current_state]; ++j){
                                if (current_state == transInStateOrder[j][2] && transInStateOrder[j][1] == 0){
                                    divergent[tos] = true;
                                    break;
                                }
                            }
                        }
                        else{
                            divergent[tos] = true;
                        }

                    }
                    stack.pop();
                }
            }
        }
    }
    void boundaries(){
        lower_bound.clear();
        upper_bound.clear();
        for(size_t i = 0; i < states.size(); i++){
            lower_bound.push_back(0);
            upper_bound.push_back(0);
        }
        for(size_t j = 0; j < transitions.size(); j++){
            transInStateOrder.push_back({0, 0, 0});
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
            transInStateOrder[indices[k[2]]] = {k[0], k[1], k[2]};
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
                if (tau.find(transInStateOrder[i][1]) != tau.end()
                    && partitions[target] == partitions[transInStateOrder[i][0]]){
                    labelInsert(transInStateOrder[i][0], event, block);
                }
            }
        }
    }
    void transReduceDSV2(){
        //First move all states into one block, and replace tau event between different blocks with non-tau event
        //lmd[0] = 0;

        for (auto & trans : transitions){
            if (trans[2] == 0 && trans[0] != 0 && trans[1] == 0){
                trans[1] = 5*10000 + 2;
            }
            else if(trans[0] == 0 && trans[2] != 0 && trans[1] == 0){
                trans[1] = 5*10000 + 1;
            }
        }
        outgoingBoundaries();
        if(divergent.size() != states.size()){
            std::vector<bool> divergentInit(states.size(), false);
            divergent = divergentInit;
        }
        if(lmd.empty()){
            std::vector<size_t > lmdInit(states.size(), 0);
            lmd = lmdInit;
        }
        //size_t_in_vec dfn(states.size(), 0), low(states.size(), 0), inStack(states.size(), 0),
        //              visited(states.size(), 0);
        //std::stack<size_t > stack;
        //size_t index = 0;
        //for (const auto s : states){
        //    tarjan(0, dfn, low, stack, inStack, visited, index);
        //}
        //for (const auto &i : states){
        //    if (!divergent[i]){
        //        for (size_t j = lower_bound[i]; j < upper_bound[i]; ++j){
        //            if (divergent[transInStateOrder[j][2]] && transInStateOrder[j][1] == 0 && lmd[i] == lmd[transInStateOrder[j][2]]){
        //                divergent[i] = true;
        //            }
        //        }
        //    }
        //}
        compute_tau_sccs(10000);
        std::set<size_t > sts;
        size_t n = states.size(), count;
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
            label_map.clear();
            for(size_t i = 0; i < n; i++){
                label.emplace_back();
            }

            for(auto current_trans : transitions){
                size_t source = current_trans[0], event = current_trans[1], target = current_trans[2];
                if (!(tau.find(event) != tau.end() && partitions[source] == partitions[target])
                    || divergent[current_trans[2]]){
                    labelInsert(source, event, partitions[target]);
                }
            }
            size_t m = states.size() + 1;
            for (auto current_state : states){
                label[current_state].insert({partitions[current_state], m});
                if (label_map.find(label[current_state]) == label_map.end()){
                    label_map[label[current_state]] = count++;
                }
            }

            for (auto current_state : states){
                partitions_new[current_state] = label_map[label[current_state]];
            }
        }while(partitions_new != partitions);
        size_t_in_vec_in_set sPi;
        for (const auto &iter : transitions) {
            if (!(partitions[iter[0]] == partitions[iter[2]] && (tau.find(iter[1]) != tau.end()))
                ||divergent[iter[2]]) {
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
        for (auto &iter : transitions){
            if(iter[1] == 5*10000+1 || iter[1] == 5*10000 +2){
                iter[1] = 0;
            }
        }
    }
    void transReduceDSV(){
        outgoingBoundaries();
        if(divergent.size() != states.size()){
            std::vector<bool> divergentInit(states.size(), false);
            divergent = divergentInit;
        }
        if(lmd.empty()){
            std::vector<size_t > lmdInit(states.size(), 0);
            lmd = lmdInit;
        }
        size_t_in_vec dfn(states.size(), 0), low(states.size(), 0), inStack(states.size(), 0),
                      visited(states.size(), 0);
        std::stack<size_t > stack;
        size_t index = 0;
        for (const auto s : states){
            tarjan(s, dfn, low, stack, inStack, visited, index);
        }
        for (const auto &i : states){
            if (!divergent[i]){
                for (size_t j = lower_bound[i]; j < upper_bound[i]; ++j){
                    if (divergent[transInStateOrder[j][2]] && transInStateOrder[j][1] == 0 && lmd[i] == lmd[transInStateOrder[j][2]]){
                        divergent[i] = true;
                    }
                }
            }
        }
        std::set<size_t > sts;
        size_t n = states.size(), count;
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
            label_map.clear();
            for(size_t i = 0; i < n; i++){
                label.emplace_back();
            }

            for(auto current_trans : transitions){
                size_t source = current_trans[0], event = current_trans[1], target = current_trans[2];
                if (!(tau.find(event) != tau.end() && partitions[source] == partitions[target])
                    || divergent[current_trans[2]]){
                    labelInsert(source, event, partitions[target]);
                }
            }
            size_t m = states.size() + 1;
            for (auto current_state : states){
                label[current_state].insert({partitions[current_state], m});
                if (label_map.find(label[current_state]) == label_map.end()){
                    label_map[label[current_state]] = count++;
                }
            }

            for (auto current_state : states){
                partitions_new[current_state] = label_map[label[current_state]];
            }
        }while(partitions_new != partitions);
        size_t_in_vec_in_set sPi;
        for (const auto &iter : transitions) {
            if (!(partitions[iter[0]] == partitions[iter[2]] && (tau.find(iter[1]) != tau.end()))
                ||divergent[iter[2]]) {
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
        for (auto &iter : transitions){
            if(iter[1] == 5*10000+1 || iter[1] == 5*10000 +2){
                iter[1] = 0;
            }
        }
    }

    void transReduce(){
        transitions.size();
        std::set<size_t > sts;
        size_t n = states.size(), count;
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
            label_map.clear();
            for(size_t i = 0; i < n; i++){
                label.emplace_back();
            }

            for(auto current_trans : transitions){
                size_t source = current_trans[0], event = current_trans[1], target = current_trans[2];
                if (!(tau.find(event) != tau.end() && partitions[source] == partitions[target])){
                    labelInsert(source, event, partitions[target]);
                }
            }
            size_t m = states.size() + 1;
            for (auto current_state : states){
                label[current_state].insert({partitions[current_state], m});
                if (label_map.find(label[current_state]) == label_map.end()){
                    label_map[label[current_state]] = count++;
                }
            }

            for (auto current_state : states){
                partitions_new[current_state] = label_map[label[current_state]];
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
    }

    std::vector<bool> coreach(){
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
    void reachBoolean(){
        size_t_in_vec x;
        for (size_t i = 1; i < ndeltaR.size() + 1; ++i){
            x.push_back(i);
        }
        std::vector<bool> b(ndeltaR.size(), false);
        std::vector<bool> bdFalse = b;
        for (const auto &i : init){
            b[i - 1] = true;
        }
        std::vector<bool> bd = b;
        size_t dMax = ndeltaR.size();

        while (std::find(bd.begin(), bd.end(), true) != bd.end()){
            std::vector<bool> bd1 = bdFalse;
            for(size_t i = 0; i < ndeltaR.size(); ++i){
                if (bd[i]){
                    size_t q = x[i];
                    for(size_t j = 0; j < ndeltaR[q - 1]; ++j){
                        bd1[deltaR[j][q - 1] - 1] = true;
                    }
                }
            }
            for(size_t i = 0; i < ndeltaR.size(); ++i){
                bd[i] = bd1[i] && !b[i];
            }
            for(size_t i = 0; i < ndeltaR.size(); ++i){
                b[i] = b[i] || bd[i];
            }
        }
        for (const auto &i : b){
            std::cout << i << "";
        }
        std::cout << std::endl;
    }
    size_t_in_set reachInT(){
        outgoingBoundaries();
        size_t_in_set reachableStateSet(init);
        std::vector<std::set<size_t> > deltaX;
        deltaX.push_back(init);
        size_t k = 0;
        do{
            ++k;
            deltaX.emplace_back();
            for (const auto currentState : deltaX[k - 1]){
                for (size_t i = lower_bound[currentState]; i < upper_bound[currentState]; ++i){
                    deltaX[k].insert(transInStateOrder[i][2]);
                }
            }
            ///size_t_in_vec reachableStateVec(reachableStateSet.begin(), reachableStateSet.end()),
            ///                deltaV(deltaX[k].begin(), deltaX[k].end()), result;
            ///std::set_difference(deltaV.begin(), deltaV.end(), reachableStateVec.begin(), reachableStateVec.end(),
            ///        std::back_inserter(result));
            size_t_in_set resultSet;
            std::set_difference(deltaX[k].begin(), deltaX[k].end(), reachableStateSet.begin(), reachableStateSet.end(),
                    std::inserter(resultSet, resultSet.end()));
            //for (const auto &state : reachableStateSet){
            //    deltaX[k].erase(state);
            //}
            deltaX[k] = resultSet;
            for (const auto &state : deltaX[k]){
                reachableStateSet.insert(state);
            }
        }while(!deltaX[k].empty());
        return reachableStateSet;
    }

    size_t_in_set reachInTSetOperation2(){
        outgoingBoundaries();
        size_t_in_set reachableStateSet(init);
        std::vector<std::set<size_t> > deltaX;
        deltaX.push_back(init);
        size_t k = 0;
        do{
            ++k;
            deltaX.emplace_back();
            for (const auto currentState : deltaX[k - 1]){
                for (size_t i = lower_bound[currentState]; i < upper_bound[currentState]; ++i){
                    deltaX[k].insert(transInStateOrder[i][2]);
                }
            }
            ///size_t_in_vec reachableStateVec(reachableStateSet.begin(), reachableStateSet.end()),
            ///                deltaV(deltaX[k].begin(), deltaX[k].end()), result;
            ///std::set_difference(deltaV.begin(), deltaV.end(), reachableStateVec.begin(), reachableStateVec.end(),
            ///        std::back_inserter(result));
            size_t_in_set resultSet, resultSet1;
            std::set_difference(deltaX[k].begin(), deltaX[k].end(), reachableStateSet.begin(), reachableStateSet.end(),
                                std::inserter(resultSet, resultSet.end()));
            //for (const auto &state : reachableStateSet){
            //    deltaX[k].erase(state);
            //}
            deltaX[k] = resultSet;
            std::set_union(reachableStateSet.begin(), reachableStateSet.end(), deltaX[k].begin(), deltaX[k].end(),
                    std::inserter(resultSet1, resultSet1.end()));
            reachableStateSet = resultSet1;
        }while(!deltaX[k].empty());
        return reachableStateSet;
    }

    std::vector<bool> reachInStack(){
        std::vector<bool> visited(states.size(), false);
        std::stack<size_t > toConsider;
        outgoingBoundaries();
        for (const auto &iter : init){
            toConsider.push(iter);
        }
        while(!toConsider.empty()){
            std::size_t currentState = toConsider.top();
            visited[currentState] = true;
            toConsider.pop();
            for (size_t i = lower_bound[currentState]; i < upper_bound[currentState]; ++i){
                if (!visited[transInStateOrder[i][2]]){
                    visited[transInStateOrder[i][2]] = true;
                    toConsider.push(transInStateOrder[i][2]);
                }
            }
        }
        for (const auto &i : visited){
            std::cout << i << " ";
        }
        std::cout << std::endl;
        return visited;
    }

    size_t_in_vec reachInTSetOp(){
        outgoingBoundaries();
        size_t_in_vec reachableStateSet(init.begin(), init.end());
        std::vector<std::vector<size_t> > deltaX;
        size_t_in_vec t, target;
        deltaX.push_back(reachableStateSet);
        size_t k = 0;
        do{
            ++k;
            deltaX.emplace_back();
            for (const auto currentState : deltaX[k - 1]){
                size_t_in_set targetS;
                for (size_t i = lower_bound[currentState]; i < upper_bound[currentState]; ++i){
                    targetS.insert(transInStateOrder[i][2]);
                }
                size_t_in_vec t0;
                size_t_in_vec targetV(targetS.begin(), targetS.end());
                std::set_union(deltaX[k].begin(), deltaX[k].end(), targetV.begin(), targetV.end(),
                               std::back_inserter(t0));
                deltaX[k] = t0;
            }
            size_t_in_vec t1, t2;
            std::set_difference(deltaX[k].begin(), deltaX[k].end(), reachableStateSet.begin(), reachableStateSet.end(),
                    std::back_inserter(t1));
            deltaX[k] = t1;
            std::set_union(reachableStateSet.begin(), reachableStateSet.end(), deltaX[k].begin(), deltaX[k].end(),
                    std::back_inserter(t2));
            reachableStateSet = t2;
        }while(!deltaX[k].empty());
        return reachableStateSet;
    }

    size_t_in_set reachInTNew(){
        outgoingBoundaries();
        size_t_in_set reachableStateSet(init);
        size_t_in_set todo = reachableStateSet;
        bool inserted;
        do{
            size_t_in_set newElements;
            inserted = false;
            for (const auto currentState : todo){
                for (size_t i = lower_bound[currentState]; i < upper_bound[currentState]; ++i){
                    if(reachableStateSet.insert(transInStateOrder[i][2]).second){
                        inserted = true;
                        newElements.insert(transInStateOrder[i][2]);
                    }
                }
            }
            todo = newElements;
        }while(inserted);
        return reachableStateSet;
    }

    size_t_in_set reachBreakPoint(){
        outgoingBoundaries();
        size_t_in_set reachableStateSet(init);
        size_t_in_set todo = reachableStateSet;
        size_t_in_set xr0 = {};
        while(xr0 != reachableStateSet){
            xr0 = reachableStateSet;
            size_t_in_set newElements;
            for (const auto currentState : todo){
                for (size_t i = lower_bound[currentState]; i < upper_bound[currentState]; ++i){
                    if(reachableStateSet.insert(transInStateOrder[i][2]).second){
                        newElements.insert(transInStateOrder[i][2]);
                    }
                }
            }
            todo = newElements;
        }
        return reachableStateSet;
    }
    void tarjanScc(size_t s){
        //DFS,always starts from state 0
        outgoingBoundaries();
        std::stack<size_t> todo;
        bool pushed;
        std::vector<bool> visited(states.size(), false);
        todo.push(s);
        visited[s] = true;
        size_t currentState;
        while(!todo.empty()){
            pushed = false;
            currentState = todo.top();
            size_t indexLow = lower_bound[currentState], indexHigh = upper_bound[currentState];
            while(indexLow < indexHigh){
                if(!visited[transInStateOrder[indexLow][2]]){
                    todo.push(transInStateOrder[indexLow][2]);
                    visited[transInStateOrder[indexLow][2]] = true;
                    pushed = true;
                    break;
                }
                else{
                    ++indexLow;
                }
            }
            if (!pushed){
                std::cout << "Pop: " << todo.top()<< std::endl;
                todo.pop();
            }
        }
    }

    void tarjan(size_t s, size_t_in_vec &dfn, size_t_in_vec &low, std::stack<size_t> &stack, size_t_in_vec &inStack,
                size_t_in_vec &visited, size_t &index){
        dfn[s] = ++index;
        low[s] = index;
        stack.push(s);
        inStack[s] = true;
        visited[s] = true;
        for (size_t i = lower_bound[s]; i < upper_bound[s]; ++i){
            if(!visited[transInStateOrder[i][2]] && transInStateOrder[i][1] == 0 && lmd[transInStateOrder[i][0]] == lmd[transInStateOrder[i][2]]){
                tarjan(transInStateOrder[i][2], dfn, low, stack, inStack, visited, index);
                low[s] = std::min(low[s], low[transInStateOrder[i][2]]);
            }
            else if(inStack[transInStateOrder[i][2]] && dfn[transInStateOrder[i][2]] < low[s]){
                low[s] = dfn[transInStateOrder[i][2]];
            }
        }
        if (dfn[s] == low[s]){
            size_t t;
            size_t_in_vec scc;
            do{
                t = stack.top();
                stack.pop();
                inStack[t] = false;
                scc.push_back(t);
                //std::cout << t << std::endl;
            }while(s != t);
            if (scc.size() == 1) {
                for (size_t i = lower_bound[s]; i < upper_bound[s]; ++i){
                    if(transInStateOrder[i][0] == transInStateOrder[i][2] && transInStateOrder[i][1] == 0){
                        divergent[s] = true;
                        break;
                    }
                }
            }
            else{
                for (const auto &i : scc){
                    divergent[i] = true;
                }
            }
        }
    }

    std::set<size_t > coreachDelta(size_t markedState){
        size_t k = 0;
        std::set<size_t > coreachableStates;
        std::vector<std::set<size_t > > deltaX;
        coreachableStates.insert(markedState);
        deltaX.push_back(coreachableStates);
        do{
            ++k;
            deltaX.emplace_back();
            for (const auto &s : deltaX[k - 1]){
                for (const auto &e : sigmaFeasibleC[s]){
                    for (const auto &c : deltaC[{s, e}]){
                        deltaX[k].insert(c);
                    }
                }
            }
            for (const auto &c : coreachableStates){
                deltaX[k].erase(c);
            }
            for (const auto &s : deltaX[k]){
                coreachableStates.insert(s);
            }
        }while(!deltaX[k].empty());
        return coreachableStates;
    }

    void compute_tau_sccs(size_t n)
    {
        std::size_t unused = 1, lastscc = 1;
        std::vector<std::size_t> scc(states.size(), 0);
        std::vector<std::size_t> low(states.size(), 0);
        std::stack<std::size_t> stack;
        std::stack<std::size_t> sccstack;

        // Record forward transition relation sorted by state.
        //outgoing_transitions_per_state_t m_lts_succ_transitions(m_lts.get_transitions(),m_lts.num_states(),true);

        for (std::size_t i = 0; i < states.size(); ++i)
        {
            if (scc[i] == 0)
                stack.push(i);
            while (!stack.empty())
            {
                const std::size_t vi = stack.top();

                // Outgoing transitions of vi.
                // std::pair<outgoing_transitions_per_state_t::const_iterator, outgoing_transitions_per_state_t::const_iterator> succ_range
                //  = m_lts_succ_transitions.equal_range(vi);

                if (low[vi] == 0 && scc[vi] == 0)
                {
                    scc[vi] = unused;
                    low[vi] = unused++;
                    sccstack.push(vi);

                    // for (outgoing_transitions_per_state_t::const_iterator t = succ_range.first; t != succ_range.second; ++t)
                    // for (const outgoing_pair_t& t: m_lts_succ_transitions[vi])
                    // for (const outgoing_pair_t& t: m_lts_succ_transitions[vi])
                    for (std::size_t i=lower_bound[vi]; i<upper_bound[vi]; ++i)
                    {
                        const size_t_in_vec & t=transInStateOrder[i];
                        if ((low[t[2]] == 0) && (scc[t[2]] == 0) && t[1] == 0)
                        {
                            stack.push(t[2]);
                        }
                    }
                }
                else
                {
                    // for (outgoing_transitions_per_state_t::const_iterator t = succ_range.first; t != succ_range.second; ++t)
                    for (std::size_t i=lower_bound[vi]; i<upper_bound[vi]; ++i)
                    {
                        const size_t_in_vec & t=transInStateOrder[i];
                        if ((low[t[2]] != 0) && t[1] == 0)
                            low[vi] = low[vi] < low[t[2]] ? low[vi] : low[t[2]];
                    }
                    if (low[vi] == scc[vi])
                    {
                        std::size_t tos, scc_id = lastscc++;
                        std::vector<std::size_t> this_scc;
                        do
                        {
                            tos = sccstack.top();
                            low[tos] = 0;
                            scc[tos] = scc_id;
                            sccstack.pop();
                            this_scc.push_back(tos);
                        }
                        while (tos != vi);

                        // We only consider non-trivial sccs, hence
                        // if the scc consists of a single schate, check whether it has a tau-loop

                        if(this_scc.size() == 1)

                        {
                            // for(outgoing_transitions_per_state_t::const_iterator i = succ_range.first; i != succ_range.second; ++i)
                            // for (const outgoing_pair_t& i: m_lts_succ_transitions[vi])
                            for (std::size_t i_=lower_bound[vi]; i_<upper_bound[vi]; ++i_)
                            {
                                const size_t_in_vec & i=transInStateOrder[i_];
                                if(vi == i[2] && i[1] == 0)
                                {
                                    divergent[tos] = true;
                                    break;
                                }
                            }
                        }

                        else
                        {
                            //std::cout << tos << std::endl;
                            divergent[tos] = true;
                        }
                    }
                    stack.pop();
                }
            }
        }
    }

    System& syncInT(System &g_to_sync);
    void outgoingBoundaries();
    void findTarget(const System&, const size_t, const size_t, const size_t, const size_t_in_vec&, size_t_in_vec_in_vec&);
    void transSigma();
    size_t_in_set TSmodel(size_t model, size_t a, size_t ncap, size_t mark);


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
