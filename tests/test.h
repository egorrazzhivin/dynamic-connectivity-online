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
#include <dynamic_connectivity_online.h>

std::mt19937 generator(std::chrono::steady_clock::now().time_since_epoch().count());
std::uniform_int_distribution<int64_t> prior(0, 1e15);

// note: i used these tests to check that time complexity is adequate
// correctness of algorithm was tested on private tests

void RunRandomTest(int n, int q) {
    std::uniform_int_distribution<int> vertex(0, n - 1);
    std::uniform_int_distribution<int> typo(0, 2);
    DynamicGraph DG = DynamicGraph(n);
    std::map<std::pair<int, int>, bool> graph;
    std::vector<int> comps;
    for (int i = 0; i < q; i++) {
        int type = typo(generator);
        if (type == 0) {
            // add
            int u = vertex(generator);
            int v = vertex(generator);
            while (v == u) {
                v = vertex(generator);
            }
            if (graph.find(std::make_pair(u, v)) == graph.end()) {
                graph[std::make_pair(u, v)] = true;
                graph[std::make_pair(v, u)] = true;
                DG.AddEdge(u, v);
            }
        } else if (type == 1) {
            // erase
            int u = vertex(generator);
            int v = vertex(generator);
            while (v == u) {
                v = vertex(generator);
            }
            if (graph.find(std::make_pair(u, v)) != graph.end()) {
                graph.erase(graph.find(std::make_pair(u, v)));
                graph.erase(graph.find(std::make_pair(v, u)));
                DG.RemoveEdge(u, v);
            }
        } else {
            // get components
            comps.push_back(DG.GetComponentsNumber());
        }
    }
}


void RunFullGraphTest(int n) {
    DynamicGraph DG = DynamicGraph(n);
    int cnt = 0;
    int limit = 300000;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            DG.AddEdge(i, j);
            cnt++;
            if (cnt > limit) {
              break;
            }
        }
    }
    std::cout << "done " << '\n';
    cnt = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            DG.RemoveEdge(i, j);
            cnt++;
            if (cnt > limit) {
              break;
        }
    }
}

void QUniqueEdges(int n, int q) {
    int need = q / 2;
    std::vector<std::pair<int, int>> edges;
    std::unordered_map<std::pair<int, int>, bool, hash> graph;
    std::uniform_int_distribution<int> vertex(0, n - 1);
    DynamicGraph DG = DynamicGraph(n);
    for (int i = 0; i < need; ++i) {
        int u = vertex(generator);
        int v = vertex(generator);
        while (v == u) {
            v = vertex(generator);
        }
        if (graph.find(std::make_pair(u, v)) == graph.end()) {
            graph[std::make_pair(u, v)] = true;
            graph[std::make_pair(v, u)] = true;
            edges.emplace_back(std::make_pair(u, v));
        }
    }
    for (auto [u, v] : edges) {
        DG.AddEdge(u, v);
    }
    std::cout << "done " << '\n';
    for (auto [u, v] : edges) {
        DG.RemoveEdge(u, v);
    }
}
