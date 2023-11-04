#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

// #define DEBUG
void help_message() {}

class DisjointSet {
   public:
    DisjointSet(int n) {
        parent = new int[n];
        rank = new int[n];
        for (int i = 0; i < n; i++) {
            parent[i] = i;
            rank[i] = 0;
        }
    }
    ~DisjointSet() {
        delete[] parent;
        delete[] rank;
    }
    int find(int x) {
        if (parent[x] != x) parent[x] = find(parent[x]);
        return parent[x];
    }
    void merge(int x, int y) {
        int xRoot = find(x);
        int yRoot = find(y);
        if (xRoot == yRoot) return;
        if (rank[xRoot] < rank[yRoot])
            parent[xRoot] = yRoot;
        else if (rank[xRoot] > rank[yRoot])
            parent[yRoot] = xRoot;
        else {
            parent[yRoot] = xRoot;
            rank[xRoot]++;
        }
    }

   private:
    int* parent;
    int* rank;
};
class edge {
   public:
    int end1, end2, weight;
    edge(){};
    edge(int _end1, int _end2, int _weight) : end1(_end1), end2(_end2), weight(_weight){};

    bool operator<(const edge& e) const { return weight < e.weight; }
    bool operator>(const edge& e) const { return weight > e.weight; }
    int operator/(const int divider) const { return weight / divider; }
    friend ostream& operator<<(ostream& os, const edge& e) {
        os << e.end1 << " " << e.end2 << " " << e.weight;
        return os;
    }
};

/*****************sorting********************/
int getMax(vector<edge>& array, int n) {
    int max = abs(array[0].weight);
    for (int i = 1; i < n; i++)
        if (abs(array[i].weight) > max) max = abs(array[i].weight);
    return max;
}
void countingSort(vector<edge>& array, int size, int place, int maxVal) {
    vector<edge> output(size);
    int count[maxVal];

    for (int i = 0; i < maxVal; ++i) count[i] = 0;
    for (int i = 0; i < size; i++) count[abs(array[i] / place) % 10]++;
    for (int i = 1; i < maxVal; i++) count[i] += count[i - 1];
    for (int i = size - 1; i >= 0; i--) {
        output[count[abs(array[i] / place) % 10] - 1] = array[i];
        count[abs(array[i] / place) % 10]--;
    }
    for (int i = 0; i < size; i++) array[i] = output[i];
}
void countingSort_sign(vector<edge>& array, int size) {
    vector<edge> output(size);
    int count[2];
    for (int i = 0; i < 2; ++i) count[i] = 0;
    for (int i = 0; i < size; i++) count[int(array[i].weight > 0)]++;
    for (int i = 1; i < 2; i++) count[i] += count[i - 1];
    int negative = count[0];
    for (int i = size - 1; i >= 0; i--) {
        output[count[int(array[i].weight > 0)] - 1] = array[i];
        count[int(array[i].weight > 0)]--;
    }
    reverse(output.begin(), output.begin() + negative);
    for (int i = 0; i < size; i++) array[i] = output[i];
}

// implement radix sort that can handdle negative integer
void radixsort(vector<edge>& array, int size, int maxVal) {
    int max = getMax(array, size);
    for (int place = 1; max / place > 0; place *= 10) countingSort(array, size, place, 10);
    countingSort_sign(array, size);
}
/********************************************/

// Print an array
void printEdges(vector<edge>& array, int size) {
    int i;
    cout << "Sorted Edges: \n";
    cout << "end1 end2 weight\n";
    for (i = 0; i < size; i++) cout << array[i] << "\n";
    cout << endl;
}
// print the path matrix
void printPath(bool** path, int size) {
    cout << "All path: \n";
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) cout << path[i][j] << " ";
        cout << "\n";
    }
    cout << endl;
}
bool check_sorted(vector<edge>& array, int size) {
    for (int i = 0; i < size - 1; i++)
        if (array[i] > array[i + 1]) return false;
    return true;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        help_message();
        return 0;
    }
    //////////// read the input file /////////////

    fstream fin(argv[1]);
    char type;
    fin >> type;
    int n, m;  // vertices, edges
    fin >> n >> m;
    ////////////     main part    ////////////////
    vector<edge> edges;
    vector<edge> discardEdges;
    int end1, end2, weight;
    for (int i = 0; i < m; i++) {
        fin >> end1 >> end2 >> weight;
        edges.push_back(edge(end1, end2, weight));
    }
    radixsort(edges, edges.size(), 100);
    // printEdges(edges, edges.size());
    // cout << "check_sorted " << check_sorted(edges, edges.size()) << endl;
    int count = 0;
    int sum = 0;
    if (type == 'u') {
        DisjointSet vertices(n);
        for (int i = edges.size() - 1; i >= 0; i--) {
            if (vertices.find(edges[i].end1) != vertices.find(edges[i].end2)) {
                vertices.merge(edges[i].end1, edges[i].end2);
                count++;
            } else {
                sum += edges[i].weight;
                discardEdges.push_back(edges[i]);
            }
            if (count == n - 1) {
                i--;
                while (i >= 0) {
                    sum += edges[i].weight;
                    discardEdges.push_back(edges[i]);
                    i--;
                }
            }
        }
    } else {
        // paths[i][j] record whether there is a path from i to j
        bool** paths = new bool*[n];
        for (int i = 0; i < n; i++) {
            paths[i] = new bool[n];
            for (int j = 0; j < n; j++) paths[i][j] = false;
            paths[i][i] = true;
        }
        DisjointSet vertices(n);
        // if i, j are in the same component, vertices.find(i) == vertices.find(j)
        // and there is no need to add a negative edge between them, because it can reduce out cost
        // but if it is positive, as long as it won't form a cycle, we can add it
        for (int i = edges.size() - 1; i >= 0; i--) {
            bool add = true;
            for (int j = 0; j < n; j++) {
                if (paths[j][edges[i].end1] && paths[edges[i].end2][j]) {
                    add = false;
                    break;
                }
            }
            if (add) {
                if (edges[i].weight < 0 &&
                    vertices.find(edges[i].end1) == vertices.find(edges[i].end2)) {
                    sum += edges[i].weight;
                    discardEdges.push_back(edges[i]);
                } else {
                    paths[edges[i].end1][edges[i].end2] = true;
                    vertices.merge(edges[i].end1, edges[i].end2);
                    for (int j = 0; j < n; j++) {
                        if (paths[j][edges[i].end1] || j == edges[i].end1) {
                            for (int k = 0; k < n; k++) {
                                paths[j][k] = paths[j][k] || paths[edges[i].end2][k];
                            }
                        }
                    }
                }
            } else {
                sum += edges[i].weight;
                discardEdges.push_back(edges[i]);
            }
        }
    }
    //////////// write the output file ///////////
    fstream fout;
    fout.open(argv[2], ios::out);
    fout << sum << "\n";
    for (int i = 0; i < discardEdges.size(); i++) {
        fout << discardEdges[i] << "\n";
    }
    fout.close();
    return 0;
}
