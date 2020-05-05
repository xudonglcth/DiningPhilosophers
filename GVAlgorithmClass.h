//
// Created by 梁旭东 on 2020-04-14.
//

#ifndef DININGPHILOSOPHERS_GVALGORITHMCLASS_H
#define DININGPHILOSOPHERS_GVALGORITHMCLASS_H

#include "System.h"
#include <unordered_set>
class GV: public System{
public:
    size_t eqIndex, maxStateIndex;
    size_t_in_vec dfs, blockIndex;
    size_t_in_vec_in_vec tgtSortedTrans, srcSortedTrans;
    void sortTauTrans(const bool outgoing){
            lower_bound.clear();
            upper_bound.clear();
            srcSortedTrans.clear();
            tgtSortedTrans.clear();
            size_t max_state = *std::max_element(states.begin(), states.end());
            lower_bound.resize(max_state + 1);
            upper_bound.resize(max_state + 1);
            std::vector<size_t > indices(max_state + 2, 0);
            for(const auto &i : transitions){
                if (i[1] == 0){
                    ++indices[outgoing ? i[0] : i[2]];
                }
            }
            size_t sum = 0;
            for (auto &j : indices){
                sum += j;
                j = sum;
            }
        outgoing ? srcSortedTrans.resize(sum) : tgtSortedTrans.resize(sum);

        for(const auto &k : transitions){
            if(k[1] == 0){
                outgoing ?
                        srcSortedTrans[--indices[k[0]]] = {k[0], k[1], k[2]}
                         :
                        tgtSortedTrans[--indices[k[2]]] = {k[0], k[1], k[2]};
            }
        }
        for(size_t i = 0; i < lower_bound.size(); ++i){
            lower_bound[i] = indices[i];
            upper_bound[i] = indices[i + 1];
        }
    }
    void reduceGV(const bool divergenceSensitive){
        //First, remove tau loops
        sccReduce(divergenceSensitive);
        replaceTransSystem(divergenceSensitive);
    }

    void sccReduce(const bool divergenceSensitive){
        sortTauTrans(true);
        dfs.reserve(states.size());
        std::vector<bool> visited(states.size(), false);
        for (size_t i = 0; i < states.size(); ++i){
            dfsEnumerate(i, visited);
        }
        sortTauTrans(false);
        for (size_t i = dfs.size() - 1; i > 0; --i){
            if (visited[i]){
                groupComponents(i, visited);
                ++eqIndex;
            }
        }
    }

    void dfsEnumerate(size_t s, std::vector<bool> & visited){
        if (visited[s]){
            return;
        }
        visited[s] = true;
        for (size_t i = lower_bound[s]; i < upper_bound[s]; ++i){
            dfsEnumerate(srcSortedTrans[i][2], visited);
        }
        dfs.push_back(s);
    }


    void groupComponents(size_t s, std::vector<bool> &visited){
        blockIndex.resize(s);
        if (!visited[s]){
            return;
        }
        visited[s] = false;
        for (size_t i = lower_bound[s]; i < upper_bound[s]; ++i){
            groupComponents(tgtSortedTrans[i][0], visited);
        }
        blockIndex[s] = eqIndex;
    }

    void replaceTransSystem(const bool divergenceSensitive){
        std::set<size_t_in_vec > resultTrans;
        for (const auto & t: transitions){
            if (t[1] != 0 || divergenceSensitive || blockIndex[t[0]] != blockIndex[t[2]]){
                resultTrans.insert({blockIndex[t[0]], blockIndex[t[1]], blockIndex[t[2]]});
            }
        }
        size_t_in_vec_in_vec swapV(resultTrans.begin(), resultTrans.end());
        transitions.swap(swapV);
        std::vector<size_t > numState(eqIndex);
        for (size_t i = 0; i < eqIndex; ++i){
            numState[i] = i;
        }
    }
    void bisimPartitioner(const bool divergenceSensitive){
        createInitialParition(divergenceSensitive);
        refinePartition(divergenceSensitive);
    }
    std::vector<size_t> blockIndexState;
    std::vector<bool> blockFlag;
    std::vector<bool> blockInToBeProcessed;
    std::vector<bool> stateFlag;

    std::vector<size_t> toBeProcessed;
    std::vector<size_t> BL;
    struct nonBottomState
    {
        size_t state;
        std::vector < size_t > inertTransitions; // Only the target state is interesting.

        nonBottomState(const size_t s)
                : state(s)
        {}
        nonBottomState(const size_t s, const std::vector < size_t > &it)
                : state(s), inertTransitions(it)
        {}
    };
    struct block{
        size_t stateIndex;
        size_t blockIndex;
        size_t parentBlockIndex;

        std::pair <size_t, size_t > splitter;

        std::vector<size_t> bottomStates;
        std::vector<nonBottomState> nonBottomStates;
        std::vector<std::vector<size_t> > nonInertTrans;

        void swap(block& b){
            size_t stateIndex1 = b.stateIndex;
            b.stateIndex = stateIndex;
            stateIndex = stateIndex1;

            size_t blockIndex1 = b.blockIndex;
            b.blockIndex = blockIndex;
            blockIndex = blockIndex1;

            size_t parentBlockIndex1 = b.parentBlockIndex;
            b.parentBlockIndex = parentBlockIndex;
            parentBlockIndex = parentBlockIndex1;

            std::pair < size_t , size_t > splitter1 = b.splitter;
            b.splitter=splitter;
            splitter = splitter1;

            bottomStates.swap(b.bottomStates);
            nonBottomStates.swap(b.nonBottomStates);
            nonInertTrans.swap(b.nonInertTrans);
        }
    };
    std::vector<block> blocks;

    void createInitialParition(const bool divergenceSensitive){
        toBeProcessed.clear();
        block initialPartition;
        sortTrans(1);
        size_t lastNonStoredStateNumber = 0;
        bool bottomState = true;
        size_t_in_vec currentInertTrans;

        {
            //reserve space, avoiding reallocation for vectors
            size_t initialPartitionNonInertCounter = 0;
            size_t currentInertTransCounter = 0;
            for (const auto& t : transitions){
                if (t[1] == 0){
                    if (divergenceSensitive && t[0] == t[2]){
                        ++initialPartitionNonInertCounter;
                    }
                    else{
                        ++currentInertTransCounter;
                    }
                }
            }
            currentInertTrans.reserve(initialPartitionNonInertCounter);
            initialPartition.nonInertTrans.reserve(currentInertTransCounter);
        }
        size_t currentTransNumber = 0;
        for (const auto& t : transitions){
            ++currentTransNumber;
            if (t[1] == 0){
                if (divergenceSensitive && t[0] == t[2]){
                    initialPartition.nonInertTrans.push_back({t[0], 0, t[2]});
                }
                else{
                    currentInertTrans.push_back(t[2]);
                    bottomState = false;
                }
            }
            if (currentTransNumber == transitions.size() ||
            transitions[currentTransNumber][0] == transitions[currentTransNumber][2]){
                for ( ; lastNonStoredStateNumber < t[0]; ++lastNonStoredStateNumber){
                    initialPartition.bottomStates.push_back(lastNonStoredStateNumber);
                }
                if (bottomState){
                    initialPartition.bottomStates.push_back(t[0]);
                }
                else{
                    initialPartition.nonBottomStates.push_back(nonBottomState(t[0]));
                    currentInertTrans.swap(initialPartition.nonBottomStates.back().inertTransitions);
                    bottomState = true;
                }
                assert(lastNonStoredStateNumber == t[0]);
                ++lastNonStoredStateNumber;
            }
        }

        for (; lastNonStoredStateNumber < states.size(); ++lastNonStoredStateNumber){
            initialPartition.bottomStates.push_back(lastNonStoredStateNumber);
        }

        orderOnTauReachability(initialPartition.nonBottomStates);
        sortTrans(0);
        const size_t_in_vec_in_vec &trans1 = transitions;
        for (const auto & t : trans1){
            if (t[1] != 0){
                initialPartition.nonInertTrans.push_back(t);
            }
        }

        initialPartition.blockIndex = 0;
        initialPartition.stateIndex = 0;
        maxStateIndex = 1;
        initialPartition.parentBlockIndex = 0;
        initialPartition.splitter = std::pair<size_t, size_t > (0, 0);
        blocks.push_back(block());
        blocks.back().swap(initialPartition);
        blockIndexState = std::vector<size_t > (states.size(), 0);
        blockFlag.push_back(false);
        stateFlag = std::vector<bool> (states.size(), false);
        blockInToBeProcessed.push_back(false);
        toBeProcessed.clear();
        BL.clear();
    }//End Create Initial Partition

    void refinePartition(const bool divergenceSensitive){
        size_t consistencyCheckCounter = 1;
        size_t consistencyCheckBarrier = 1;
        bool partitionIsUnstable = true;
        while (!toBeProcessed.empty() || partitionIsUnstable){
            ++consistencyCheckCounter;
            //Check temporarily skipped
            if (consistencyCheckCounter >= consistencyCheckBarrier){
                consistencyCheckCounter = 0;
                consistencyCheckBarrier *= 2;
                checkInternalConsistencyOfThePartitioningDataStruct();
            }
            if (toBeProcessed.empty() && partitionIsUnstable){
                //All blocks put in toBeProcessed
                for (size_t i = 0 ; i < blocks.size(); ++i){
                    toBeProcessed.push_back(i);
                }
                blockInToBeProcessed = std::vector<bool>(blocks.size(), true);
                partitionIsUnstable = false;
            }
            const size_t splitterIndex = toBeProcessed.back();
            toBeProcessed.pop_back();
            blockInToBeProcessed[splitterIndex] = false;

            for (size_t i = 0; i < blocks[splitterIndex].nonInertTrans.size(); ++i){
                stateFlag[blocks[splitterIndex].nonInertTrans[i][0]] = true;
                const size_t markedBlockIndex =
                        blockIndexState[blocks[splitterIndex].nonInertTrans[i][0]];
                if (!blockFlag[markedBlockIndex]){
                    blockFlag[markedBlockIndex] = true;
                    BL.push_back(markedBlockIndex);
                }

                if((i + 1) == blocks[splitterIndex].nonInertTrans.size()||
                blocks[splitterIndex].nonInertTrans[i][1] !=
                blocks[splitterIndex].nonInertTrans[i+1][1]){
                    splitBlocksInBL(partitionIsUnstable, blocks[splitterIndex].nonInertTrans[i][1],
                            splitterIndex);
                }
            }
        }
        blockFlag.clear();
        blockInToBeProcessed.clear();
        stateFlag.clear();
        toBeProcessed.clear();
        BL.clear();
    }
    void splitBlocksInBL(
            bool& partitionIsUnstable,
            const size_t splitterLabel,
            const size_t splitterBlock) {
        for (const auto &BLElement : BL) {
            blockFlag[BLElement] = false;
            std::vector<size_t> flaggedStates;
            std::vector<size_t> nonFlaggedStates;
            std::vector<size_t> i1BottomStates;
            i1BottomStates.swap(blocks[BLElement].bottomStates);

            for (const auto &i1Element : i1BottomStates) {
                if (stateFlag[i1Element]) {
                    flaggedStates.push_back(i1Element);
                } else {
                    nonFlaggedStates.push_back(i1Element);
                    blockIndexState[i1Element] = blocks.size();
                }
            }
            size_t resetStateFlagsBlock = BLElement;

            if (!nonFlaggedStates.empty()) {
                blocks[BLElement].splitter =
                        std::pair<size_t, size_t>(splitterLabel, splitterBlock);

                //create first new block
                blocks.push_back(block());
                size_t newBlock1 = blocks.size() - 1;
                blocks.back().blockIndex = newBlock1;
                blocks.back().parentBlockIndex = BLElement;

                //put the indices of first split block to toBeProcessed
                nonFlaggedStates.swap(blocks.back().bottomStates);
                toBeProcessed.push_back(blocks.back().blockIndex);
                blockInToBeProcessed.push_back(true);

                //create second new block
                blocks.push_back(block());
                size_t newBlock2 = blocks.size() - 1;
                blocks.back().stateIndex = blocks[BLElement].stateIndex;
                blocks.back().blockIndex = newBlock2;
                resetStateFlagsBlock = newBlock2;
                blocks.back().parentBlockIndex = BLElement;

                //Move flagged states to the second block,
                //and let block indices of these states refer to this block
                std::vector<size_t> &reference2FlaggedStatesOfBlock2 =
                        blocks.back().bottomStates;
                for (auto j = reference2FlaggedStatesOfBlock2.begin();
                     j != reference2FlaggedStatesOfBlock2.end(); ++j) {
                    blockIndexState[*j] = newBlock2;
                }

                toBeProcessed.push_back(blocks.back().blockIndex);
                blockInToBeProcessed.push_back(true);

                blockInToBeProcessed[BLElement] = false;

                blockFlag.push_back(false);
                blockFlag.push_back(false);

                std::vector<size_t_in_vec> flaggedNonInertTransitions;
                std::vector<size_t_in_vec> nonFlaggedNonInertTransitions;

                flaggedNonInertTransitions.reserve(blocks[BLElement].nonInertTrans.size());
                nonFlaggedNonInertTransitions.reserve(blocks[BLElement].nonInertTrans.size());

                std::vector<nonBottomState> flaggedNonBottomStates;
                std::vector<nonBottomState> nonFlaggedNonBottomStates;
                std::vector<nonBottomState> i1NonBottomStates;
                i1NonBottomStates.swap(blocks[BLElement].nonBottomStates);

                for (auto k = i1NonBottomStates.begin();
                     k != i1NonBottomStates.end(); ++k) {
                    const std::vector<size_t> &inertTransitions = k->inertTransitions;
                    if (!stateFlag[k->state]) {
                        bool allTransitionsEndInUnflaggedBlock = true;
                        for (auto l = inertTransitions.begin();
                             allTransitionsEndInUnflaggedBlock && l != inertTransitions.end();
                             ++l) {
                            if (blockIndexState[*l] != newBlock1) {
                                blockIndexState[*l] = newBlock2;
                                allTransitionsEndInUnflaggedBlock = false;
                            }
                        }
                        if (allTransitionsEndInUnflaggedBlock) {
                            nonBottomState s(k->state);
                            s.inertTransitions.swap(k->inertTransitions);
                            nonFlaggedNonBottomStates.push_back(s);
                            blockIndexState[k->state] = newBlock1;
                            continue;
                        }
                    }

                    size_t_in_vec remainingInertTransitions;
                    for (auto l = inertTransitions.begin();
                         l != inertTransitions.end(); ++l) {
                        if (blockIndexState[*l] == newBlock1) {
                            nonFlaggedNonInertTransitions.push_back({k->state, 0, *l});
                        } else {
                            blockIndexState[*l] = newBlock2;
                            remainingInertTransitions.push_back(*l);
                        }
                    }
                    if (remainingInertTransitions.empty()) {
                        blocks[newBlock2].bottomStates.push_back(k->state);
                        blockIndexState[k->state] = newBlock2;
                        partitionIsUnstable = true;
                    } else {
                        flaggedNonBottomStates.emplace_back(k->state, remainingInertTransitions);
                        blockIndexState[k->state] = newBlock2;
                    }
                }
                nonFlaggedNonBottomStates.swap(blocks[newBlock1].nonBottomStates);
                flaggedNonBottomStates.swap(blocks[newBlock2].nonBottomStates);

                std::vector<size_t_in_vec> i1NonInerttransitions;
                i1NonInerttransitions.swap(blocks[BLElement].nonInertTrans);
                for (const auto &k : i1NonInerttransitions) {
                    if (blockIndexState[k[2]] == newBlock1) {
                        nonFlaggedNonInertTransitions.push_back(k);
                    } else {
                        blockIndexState[k[2]] = newBlock2;
                        flaggedNonInertTransitions.push_back(k);
                    }
                }

                blocks[newBlock1].nonInertTrans.reserve(nonFlaggedNonInertTransitions.size());
                for (auto i = nonFlaggedNonInertTransitions.begin();
                     i != nonFlaggedNonInertTransitions.end(); ++i) {
                    blocks[newBlock1].nonInertTrans.push_back(*i);
                }

                blocks[newBlock2].nonInertTrans.reserve(flaggedNonInertTransitions.size());
                for (auto i = flaggedNonInertTransitions.begin();
                     i != flaggedNonInertTransitions.end(); ++i) {
                    blocks[newBlock2].nonInertTrans.push_back(*i);
                }
            } else {
                i1BottomStates.swap(blocks[BLElement].bottomStates);
            }
            size_t_in_vec &flaggedStatesToBeUnflagged = blocks[resetStateFlagsBlock].bottomStates;
            for (auto j = flaggedStatesToBeUnflagged.begin();
                 j != flaggedStatesToBeUnflagged.end(); ++j) {
                stateFlag[*j] = false;
            }

            std::vector<nonBottomState> &flaggedStatesToBeUnflagged1 =
                    blocks[resetStateFlagsBlock].nonBottomStates;
            for (auto j = flaggedStatesToBeUnflagged1.begin(); j != flaggedStatesToBeUnflagged1.end(); ++j) {
                stateFlag[j->state] = false;
            }
        }
        BL.clear();
    }
    void checkInternalConsistencyOfThePartitioningDataStruct(){

    }

    inline void sortTrans(const int sortType){
        class compareTransSrcTgt{
        public:
            bool operator()(const size_t_in_vec& t1, size_t_in_vec& t2){
                if (t1[0] != t2[0]){
                    return t1[0] < t2[0];
                }
                else{
                    if (t1[1] != t2[1]){
                        return t1[1] < t2[1];
                    }
                    else{
                        return t1[2] < t2[2];
                    }
                }
            }
        };

        class compareTransLabel{
        public:
            bool operator()(const size_t_in_vec& t1, size_t_in_vec& t2) {
                if (t1[1] != t2[1]){
                    return t1[1] < t2[1];
                }
                else if (t1[2] != t2[2]){
                    return t1[2] < t2[2];
                }
                else{
                    return t1[0] < t2[0];
                }
            }
        };

        switch (sortType){
            case 0:{
                compareTransLabel compare;
                std::sort(transitions.begin(), transitions.end(), compare);
                break;
            }
            case 1:
            default:{
                compareTransSrcTgt compare;
                std::sort(transitions.begin(), transitions.end(), compare);
                break;
            }
        }
    }

    void orderOnTauReachability(std::vector<nonBottomState> &nonBottomStates){
        std::set<size_t > visited;
        std::map<size_t, std::vector<size_t > > inertTransitionMap;
        for (auto &i : nonBottomStates){
            i.inertTransitions.swap(inertTransitionMap[i.state]);
        }
        std::vector<nonBottomState> newNonBottomStates;
        for (const auto& nbs : nonBottomStates){
            orderRecursivelyOnTauReachability(nbs.state, inertTransitionMap, newNonBottomStates, visited);
        }
        newNonBottomStates.swap(nonBottomStates);
    }

    void orderRecursivelyOnTauReachability(
            const size_t s, std::map<size_t, std::vector<size_t >> & inertTransitionMap,
            std::vector<nonBottomState>& newNonBottomStates,
            std::set<size_t > visited
            ){
        if (inertTransitionMap.count(s) > 0){
            if (visited.count(s) == 0){
                visited.insert(s);
                std::vector<size_t >& inertTransitions = inertTransitionMap[s];
                for (const auto& j : inertTransitions){
                    orderRecursivelyOnTauReachability(j, inertTransitionMap, newNonBottomStates, visited);
                }
                newNonBottomStates.push_back(nonBottomState(s));
                inertTransitions.swap(newNonBottomStates.back().inertTransitions);
            }
        }
    }

};
#endif //DININGPHILOSOPHERS_GVALGORITHMCLASS_H
