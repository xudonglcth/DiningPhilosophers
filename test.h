//
// Created by 梁旭东 on 2020-02-24.
//

#ifndef TEST_H
#define TEST_H

#include <utility>
#include <vector>
#include <set>
#include <unordered_map>
#include <boost/functional/hash.hpp>
using Label = std::set<std::pair<size_t, size_t >>;

    template <typename Container> // we can make this generic for any container [1]
    struct container_hash {
        std::size_t operator()(Container const& c) const {
            return boost::hash_range(c.begin(), c.end());
        }
    };

void testUnorderedMap(){
    std::unordered_map<Label, size_t, container_hash<Label>> map;
    std::set<std::pair<size_t, size_t>> s = {{1, 2}, {2, 3}, {3, 4}};
    Label ll(s);
    map[ll] = 10;
    std::cout << map[ll] << std::endl;
}
#endif //TEST_H
