#pragma once

#include <vector>

/// Union-find with non-comutative union operation
template<typename size_n_max>
class UnionFind {
    std::vector<size_n_max> roots;
    size_n_max component_count;
public:
    UnionFind(size_n_max size) : roots(size), component_count(size) {
        for (size_n_max i = 0; i < size; ++i) roots[i] = i;
    };
    UnionFind() = default;
    
    inline size_n_max find(size_n_max x) {
        size_n_max root = roots[x];
        if (roots[root] == root) return root;

        while (roots[root] != root) root = roots[root];
        while (x != root){
            size_n_max new_x = roots[x];
            roots[x] = root;
            x = new_x;
        }
        return root;
    }

    inline bool are_connected(size_n_max x, size_n_max y){
        return find(x) == find(y);
    }

    inline void connect(size_n_max to, size_n_max from){ // Second one points to the first one - points to the begining of a chain
        if (are_connected(from, to)) return;
        roots[from] = to;
        --component_count;
    }

    inline size_n_max count() const { return component_count; };
};