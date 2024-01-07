#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <tuple>
#include <memory>
#include <random>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <chrono>

std::mt19937 generator(std::chrono::steady_clock::now().time_since_epoch().count());
std::uniform_int_distribution<int64_t> prior(0, 1e15);

// dynamic euler tour tree using treaps with implicit keys

struct hash {
    int64_t operator()(const std::pair<int, int>& p) const {
        return static_cast<int64_t>((static_cast<int64_t>(p.first) << 20) ^ p.second);
    }
};


struct Node {

    /*
        key - edge u-v
        size - size of subtree
        size_of_min_level - number of edges u-v (u < v) in subtree with min level
        (actually it's maximal level, not minimal, but nvm)
        size_of_adjacent - number of vertices u in subtree such that 
        u has at least one edge to v and u-v has min level
        level - level of edge
        is_min_level - is this node has min level (only for u-v)
        is_has_adjacent - is this node u has adjacent vertex v so that
        level of u-v is minimal
        priority - node's priority in treap
        left, right, parent - pointers to left subtree / right subtree / parent
    */

    std::pair<int, int> key;
    int size;
    bool size_of_min_level;
    bool size_of_adjacent;
    bool is_min_level;
    bool is_has_adjacent;
    int level;
    int64_t priority;
    Node* left;
    Node* right;
    Node* parent;

    Node() : left(nullptr), right(nullptr), parent(nullptr) {
        size = 0;
        size_of_adjacent = false;
        size_of_min_level = false;
        is_min_level = false;
        is_has_adjacent = false;
        level = 0;
    }

    Node(std::pair<int, int> key, int64_t priority, int lvl)
        : key(key),
          priority(priority),
          left(nullptr),
          right(nullptr),
          parent(nullptr) {
        size = 1;
        size_of_adjacent = false;
        size_of_min_level = false;
        is_min_level = false;
        is_has_adjacent = false;
        level = lvl;
    }
};

inline int get_size(Node* root) {
    if (root == nullptr) {
        return 0;
    }
    return root->size;
}

inline bool get_size_min_level(Node* root) {
    if (root == nullptr) {
        return 0;
    }
    return root->size_of_min_level;
}

inline bool get_size_adjacent(Node* root) {
    if (root == nullptr) {
        return 0;
    }
    return root->size_of_adjacent;
}

inline void update_size(Node* root) {
    if (root) {
        root->size = get_size(root->left) + get_size(root->right) + 1;
    }
}

inline void update_size_flag(Node* root) {
    if (root) {
        root->size_of_min_level = (get_size_min_level(root->left) |
                                  get_size_min_level(root->right) |
                                  root->is_min_level);
        root->size_of_adjacent = (get_size_adjacent(root->left) |
                                 get_size_adjacent(root->right) |
                                 root->is_has_adjacent);
    }
}

inline void update_up(Node* root) {
    update_size_flag(root);
    if (root->parent) {
        update_up(root->parent);
    }
}

inline void split(Node* root, int key,
           Node*& left, Node*& right) {
    if (root == nullptr) {
        left = nullptr;
        right = nullptr;
        return;
    }
    if (get_size(root->left) >= key) {
        split(root->left, key, left, root->left);
        if (left) {
            left->parent = nullptr;
        }
        if (root->left) {
            root->left->parent = root;
        }
        right = root;
    } else {
        split(root->right, key - get_size(root->left) - 1, root->right, right);
        if (right) {
            right->parent = nullptr;
        }
        if (root->right) {
            root->right->parent = root;
        }
        left = root;
    }   
    update_size(left);
    update_size_flag(left);
    update_size(right);
    update_size_flag(right);
}

inline void merge(Node*& root, Node* left, Node* right) {
    if (left == nullptr) {
        root = right;
        if (root) {
            root->parent = nullptr;
        }
        update_size(root);
        update_size_flag(root);
        return;
    }
    if (right == nullptr) {
        root = left;
        if (root) {
            root->parent = nullptr;
        }
        update_size(root);
        update_size_flag(root);
        return;
    }
    if (left->priority < right->priority) {
        merge(left->right, left->right, right);
        root = left;
    } else {
        merge(right->left, left, right->left);
        root = right;
    }
    if (root->left) {
        root->left->parent = root;
    }
    if (root->right) {
        root->right->parent = root;
    }
    update_size(root);
    update_size_flag(root);
}

inline Node* lift(Node* root) {
    while (root && root->parent) {
        root = root->parent;
    }
    return root;
}

// finds implicit key

inline void get_normal_key(Node* parent,
                           Node* start_loop, int& result) {
    if (!parent) {
        return;
    }
    if (parent->right == start_loop) {
        result += get_size(parent->left) + 1;
    }
    get_normal_key(parent->parent, parent, result);
}

// reroot euler tour tree (helpful function for merging / spliting two ETT)

inline void reroot(Node*& root, int start, int end,
            std::unordered_map<std::pair<int, int>, Node*, hash>& map_edges) {
    auto start_loop = map_edges[{start, end}];
    int key = get_size(start_loop->left);
    get_normal_key(start_loop->parent, start_loop, key);
    Node* first = nullptr;
    split(root, key, root, first);
    merge(root, first, root);
}

/* 
   in accordance with the article, DynamicForest_{i} (= F_{i} in article) 
   is spanning tree consisting of edges u-v such that level(u-v) <= i
*/

struct DynamicForest {

    /*
        map_edges - map that stores pointer to place 
        of edge u-v inside spanning tree
        adjacent_edges - map that stores adjacent
        edges with needed level
        level - level of DynamicForest
    */

    std::unordered_map<std::pair<int, int>, Node*, hash> map_edges;
    std::unordered_map<int, std::unordered_set<int>> adjacent_edges;
    int level;
    DynamicForest(int nn, int level) : level(level) {
        for (int i = 0; i < nn; ++i) {
            int64_t value = prior(generator);
            Node* it = new Node({i, i}, value, -1);
            map_edges[{i, i}] = it;
        }
    }

    ~DynamicForest() {
        for (const auto& edge : map_edges) {
            delete edge.second;
        }
    }

    bool is_connected(int uu, int vv) {
        return lift(map_edges[{uu, uu}]) == lift(map_edges[{vv, vv}]);
    }

    void add_edge(int uu, int vv, int lvl) {
        auto left = lift(map_edges[{uu, uu}]);
        reroot(left, uu, uu, map_edges);
        auto right = lift(map_edges[{vv, vv}]);
        reroot(right, vv, vv, map_edges);
        int64_t value = prior(generator);
        Node* to = new Node({uu, vv}, value, lvl);
        to->is_min_level = (level == to->level && uu < vv); // u < v is important
        value = prior(generator);
        Node* from = new Node({vv, uu}, value, lvl);
        from->is_min_level = (level == from->level && vv < uu); // u < v is important
        map_edges[{uu, vv}] = to;
        map_edges[{vv, uu}] = from;
        merge(left, left, to);
        merge(left, left, right);
        merge(left, left, from);
    }

    void delete_edge(int uu, int vv) {
        auto treap = lift(map_edges[{uu, vv}]);
        reroot(treap, uu, vv, map_edges);
        Node* temporary = nullptr;
        split(treap, 1, temporary, treap);
        auto vu = map_edges[{vv, uu}];
        int key = get_size(vu->left);
        get_normal_key(vu->parent, vu, key);
        Node* left_first = nullptr;
        split(treap, key, left_first, treap);
        Node* left_second = nullptr;
        split(treap, 1, left_second, treap);
        delete vu;
        delete temporary;
        map_edges.erase({uu, vv});
        map_edges.erase({vv, uu});
    }
};

// graph G_{i} (as in article)

class DynamicGraph {
public:
    /*
        mx_level - maximal level across all edges
        n_ - number of vertices
        spanning_trees - vector of pointers to different DynamicForests
        spanning_edges_levels - map to store levels of spanning tree edges
        not_spanning_edges_levels - map to store levels of non-spanning tree edges
    */

    int mx_level = 0;
    int n_;
    int components;
    std::vector<std::unique_ptr<DynamicForest>> spanning_trees;
    std::unordered_map<std::pair<int, int>, int, hash> spanning_edges_levels;
    std::unordered_map<std::pair<int, int>, int, hash> not_spanning_edges_levels;

    explicit DynamicGraph(int nn) : n_(nn) {
        components = nn;
        build();
    }

    void build(int level = 0) {
        if (level == static_cast<int>(spanning_trees.size())) {
            spanning_trees.emplace_back(new DynamicForest(n_, level));
        }
    }

    // add edge (as in article)

    void AddEdge(int u_, int v_) {
        bool connected = spanning_trees[0]->is_connected(u_, v_);
        if (connected) {
            not_spanning_edges_levels[{u_, v_}] = 0;
            not_spanning_edges_levels[{v_, u_}] = 0;
            spanning_trees[0]->adjacent_edges[u_].insert(v_);
            spanning_trees[0]->adjacent_edges[v_].insert(u_);
            auto uu = spanning_trees[0]->map_edges[{u_, u_}];
            auto vv = spanning_trees[0]->map_edges[{v_, v_}];
            if (uu->is_has_adjacent == false) {
                uu->is_has_adjacent = true;
                update_up(uu);
            }
            if (vv->is_has_adjacent == false) {
                vv->is_has_adjacent = true;
                update_up(vv);
            }
        } else {
            --components;
            spanning_edges_levels[{u_, v_}] = 0;
            spanning_edges_levels[{v_, u_}] = 0;
            spanning_trees[0]->add_edge(u_, v_, 0);
        }
        return;
    }

    // in article, we should increase level of every edge that 
    // is not suitable but that we have visited

    void IncreaseLevel(Node* root) {
        if (!root) {
            return;
        }
        if (get_size_min_level(root)) {
            if (root->is_min_level) {
                root->is_min_level = false;
                auto key = root->key;
                int new_level = root->level + 1;
                mx_level = std::max(mx_level, new_level);
                if (new_level == static_cast<int>(spanning_trees.size())) {
                    build(new_level);
                }
                int u_ = key.first, v_ = key.second;
                spanning_trees[new_level]->add_edge(u_, v_, new_level);
                ++spanning_edges_levels[{u_, v_}];
                ++spanning_edges_levels[{v_, u_}];
            }
            IncreaseLevel(root->left);
            IncreaseLevel(root->right);
            update_up(root);
        } else {
            return;
        }
    }

    // walk around treap and visit only good nodes (nodes where we can find important edges)

    void BruteforceAdjacentEdges(Node* root,
                                 std::pair<int, int>& result, int level) {
        if (!root) {
            return;
        }
        if (result != std::make_pair(-1, -1)) {
            return;
        }
        if (get_size_adjacent(root)) {
            if (root->is_has_adjacent) {
                int u_ = root->key.first;
                std::vector<int> to_delete;
                for (auto to : spanning_trees[level]->adjacent_edges[u_]) {
                    if (spanning_trees[0]->is_connected(to, u_)) {
                        to_delete.emplace_back(to);
                    } else {
                        result = std::make_pair(u_, to);
                        spanning_trees[level]->adjacent_edges[u_].erase(to);
                        spanning_trees[level]->adjacent_edges[to].erase(u_);
                        if (spanning_trees[level]->adjacent_edges[u_].empty()) {
                            auto it = spanning_trees[level]->map_edges[{u_, u_}];
                            it->is_has_adjacent = false;
                        }
                        if (spanning_trees[level]->adjacent_edges[to].empty()) {
                            auto to_it = spanning_trees[level]->map_edges[{to, to}];
                            to_it->is_has_adjacent = false;
                            update_up(to_it);
                        }
                        break;
                    }
                }
                for (auto to : to_delete) {
                    int new_level = level + 1;
                    mx_level = std::max(mx_level, new_level);
                    if (new_level == static_cast<int>(spanning_trees.size())) {
                        build(new_level);
                    }
                    spanning_trees[level]->adjacent_edges[u_].erase(to);
                    spanning_trees[level]->adjacent_edges[to].erase(u_);
                    spanning_trees[new_level]->adjacent_edges[u_].insert(to);
                    spanning_trees[new_level]->adjacent_edges[to].insert(u_);
                    auto ffu = spanning_trees[new_level]->map_edges[{u_, u_}];
                    if (ffu->is_has_adjacent == false) {
                        ffu->is_has_adjacent = true;
                        update_up(ffu);
                    }
                    auto ffto = spanning_trees[new_level]->map_edges[{to, to}];
                    if (ffto->is_has_adjacent == false) {
                        ffto->is_has_adjacent = true;
                        update_up(ffto);
                    }
                    ++not_spanning_edges_levels[{u_, to}];
                    ++not_spanning_edges_levels[{to, u_}];
                    if (spanning_trees[level]->adjacent_edges[to].empty()) {
                        auto to_it = spanning_trees[level]->map_edges[{to, to}];
                        to_it->is_has_adjacent = false;
                        update_up(to_it);
                    }
                }
                if (spanning_trees[level]->adjacent_edges[u_].empty()) {
                    auto it = spanning_trees[level]->map_edges[{u_, u_}];
                    it->is_has_adjacent = false;
                }
            }
            BruteforceAdjacentEdges(root->left, result, level);
            BruteforceAdjacentEdges(root->right, result, level);
            update_up(root);
        } else {
            return;
        }
    }

    // finds new edge from not spanning edges that can replace deleted one

    void FindNewEdge(int u_, int v_, int level, bool& okay) {
        if (okay) {
            return;
        }
        Node* u_pointer = spanning_trees[level]->map_edges[{u_, u_}];
        Node* v_pointer = spanning_trees[level]->map_edges[{v_, v_}];
        u_pointer = lift(u_pointer);
        v_pointer = lift(v_pointer);
        if (get_size(u_pointer) > get_size(v_pointer)) {
            std::swap(u_pointer, v_pointer);
        }
        IncreaseLevel(u_pointer);
        std::pair<int, int> result = {-1, -1};
        BruteforceAdjacentEdges(u_pointer, result, level);
        if (result == std::make_pair(-1, -1)) {
            if (level == 0) {
                return;
            } else {
                FindNewEdge(u_, v_, level - 1, okay);
            }
        } else {
            okay = true;
            spanning_edges_levels[result] = level;
            spanning_edges_levels[std::make_pair(result.second, result.first)] = level;
            not_spanning_edges_levels.erase(result);
            not_spanning_edges_levels.erase({result.second, result.first});
            for (int lvl = level; lvl >= 0; --lvl) {
                spanning_trees[lvl]->add_edge(result.first, result.second, level);
            }
        }
    }

    // delete edge (as in article)

    void RemoveEdge(int u_, int v_) {
        if (not_spanning_edges_levels.count(std::make_pair(u_, v_))) {
            int current_level = not_spanning_edges_levels[std::make_pair(u_, v_)];
            not_spanning_edges_levels.erase({u_, v_});
            not_spanning_edges_levels.erase({v_, u_});
            spanning_trees[current_level]->adjacent_edges[u_].erase(v_);
            spanning_trees[current_level]->adjacent_edges[v_].erase(u_);
            if (spanning_trees[current_level]->adjacent_edges[u_].empty()) {
                auto uu = spanning_trees[current_level]->map_edges[{u_, u_}];
                uu->is_has_adjacent = false;
                update_up(uu);
            }
            if (spanning_trees[current_level]->adjacent_edges[v_].empty()) {
                auto vv = spanning_trees[current_level]->map_edges[{v_, v_}];
                vv->is_has_adjacent = false;
                update_up(vv);
            }
        } else if (spanning_edges_levels.count({u_, v_})) {
            int current_level = spanning_edges_levels[{u_, v_}];
            for (int lvl = current_level; lvl >= 0; --lvl) {
                spanning_trees[lvl]->delete_edge(u_, v_);
            }
            bool okay = false;
            FindNewEdge(u_, v_, current_level, okay);
            spanning_edges_levels.erase(std::make_pair(u_, v_));
            spanning_edges_levels.erase(std::make_pair(v_, u_));
            if (okay == false) {
                ++components;
            }
        } else {
            std::cout << "impossible" << '\n';
        }
        return;
    }

    bool IsConnected(int u_, int v_) {
        return spanning_trees[0]->is_connected(u_, v_);
    }

    int GetComponentsNumber() const {
        return components;
    }

    int GetMax() const {
        return mx_level;
    }
};
