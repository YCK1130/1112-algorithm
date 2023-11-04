#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

#define __RECURSIVE__
// #define DEBUG
void help_message() {}

template <class T>
void printMatrix(ostream& out, T** Matrix, int N, int M) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M - i; j++) {
            out << " " << Matrix[i][j];
        }
        out << "\n";
    }
    out << "\n";
    return;
}

void findChordBottomUp(ostream& out, char** MaxCut, int** maxsub, const vector<int>& chords,
                       vector<int>& result, const int& smallestLen, int N, int i, int j) {
    if (i >= j) return;
    int realj = j - i;
    if (maxsub[i][realj] == 0) return;
    if (MaxCut[i][realj] == 'N') {
        findChordBottomUp(out, MaxCut, maxsub, chords, result, smallestLen, N, i, j - 1);
        return;
    }
    if (MaxCut[i][realj] == 'I') {
        result.push_back(chords[j]);
        findChordBottomUp(out, MaxCut, maxsub, chords, result, smallestLen, N, i, chords[j] - 1);
        findChordBottomUp(out, MaxCut, maxsub, chords, result, smallestLen, N, chords[j] + 1,
                          j - 1);
        return;
    }
}
void findChordTopDown(ostream& out, char** MaxCut, int** maxsub, const vector<int>& chords,
                      vector<int>& result, const int& smallestLen, int N, int i, int j) {
    if (i >= j || i < 0 || j >= N) return;
    if (j - i < smallestLen) return;
    int realj = j - i - smallestLen;
    if (maxsub[i][realj] == 0) return;
    if (MaxCut[i][realj] == 'N') {
        findChordTopDown(out, MaxCut, maxsub, chords, result, smallestLen, N, i + 1, j);
        return;
    }
    if (MaxCut[i][realj] == 'I') {
        result.push_back(i);
        findChordTopDown(out, MaxCut, maxsub, chords, result, smallestLen, N, i + 1, chords[i] - 1);
        findChordTopDown(out, MaxCut, maxsub, chords, result, smallestLen, N, chords[i] + 1, j);
        return;
    }
}
int calculateMaxSubset_TopDown(char** maxCut, int** maxSubset, const vector<int>& chords,
                               const int& N, const int& smallestLen, const int i, const int j) {
    if (i >= j || i < 0 || j >= N) return 0;
    if (j - i < smallestLen) return 0;
    int realj = j - i - smallestLen;
    if (maxCut[i][realj] == 'I' || maxCut[i][realj] == 'N') return maxSubset[i][realj];

#ifdef DEBUG
    cerr << "entering i: " << i << ", j: " << j << endl;
#endif
    int k = chords[i];
    if ((k > i && k <= j) &&
        (calculateMaxSubset_TopDown(maxCut, maxSubset, chords, N, smallestLen, i + 1, k - 1) + 1 +
         calculateMaxSubset_TopDown(maxCut, maxSubset, chords, N, smallestLen, k + 1, j)) >
            calculateMaxSubset_TopDown(maxCut, maxSubset, chords, N, smallestLen, i + 1, j)) {
        maxSubset[i][realj] =
            calculateMaxSubset_TopDown(maxCut, maxSubset, chords, N, smallestLen, i + 1, k - 1) +
            1 + calculateMaxSubset_TopDown(maxCut, maxSubset, chords, N, smallestLen, k + 1, j);
        maxCut[i][realj] = 'I';  // In
    } else {
        maxSubset[i][realj] =
            calculateMaxSubset_TopDown(maxCut, maxSubset, chords, N, smallestLen, i + 1, j);
        maxCut[i][realj] = 'N';  // None
    }
#ifdef DEBUG
    cerr << maxCut[i][realj] << endl;
    printMatrix(cerr, maxSubset, N, N);
    cerr << "ending i: " << i << ", j: " << j << '\n' << endl;
#endif

    return maxSubset[i][realj];
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        help_message();
        return 0;
    }
    //////////// read the input file /////////////

    char buffer[200];
    fstream fin(argv[1]);

    int verticeNum;
    fin >> verticeNum;

    int end1, end2;
    vector<int> chords;
    int N = verticeNum;
    chords.resize(verticeNum);

    int smallestLen = N;
    for (size_t i = 0; i < N / 2; i++) {
        fin >> end1 >> end2;
        chords[end1] = end2;  // chords[end1] chords[end2] will form a chord
        chords[end2] = end1;
        if (abs(end1 - end2) < smallestLen) smallestLen = abs(end1 - end2);
    }
    fin.close();

    int** maxSubset = new int*[N];
    char** maxCut = new char*[N];
    for (int i = 0; i < N; i++) {
#ifndef __RECURSIVE__
        maxSubset[i] = new int[N - i];
        maxCut[i] = new char[N - i];
#else
        if (N - i - smallestLen <= 0) continue;
        maxSubset[i] = new int[N - i - smallestLen];
        maxCut[i] = new char[N - i - smallestLen];
#endif

        maxSubset[i][0] = 0;
        maxCut[i][0] = 'X';
    }
#ifdef DEBUG
    cout << "init done" << endl;
#endif
    ////////////     main part    ////////////////
    vector<int> result;
    result.reserve(max(N / 10, 128));
#ifdef __RECURSIVE__
    // Recursive (Top Down)(use the Right end's chord to divide space)
    calculateMaxSubset_TopDown(maxCut, maxSubset, chords, N, smallestLen, 0, N - 1);
#ifdef DEBUG
    cout << "cal done" << endl;
#endif
    findChordTopDown(cout, maxCut, maxSubset, chords, result, smallestLen, N, 0, N - 1);
#else
    // iteration (Bottom up)(use the left end's chord to divide space)
    for (int j = 1; j < verticeNum; j++) {
        int k = chords[j];
        for (int i = 0; i < j; i++) {
            // the chords containing the left end is in the internal
            if ((k >= i && k < j) &&
                (maxSubset[i][k - 1 - i] + 1 + maxSubset[k + 1][j - 1 - (k + 1)]) >
                    maxSubset[i][j - 1 - i]) {
                maxSubset[i][j - i] =
                    maxSubset[i][k - 1 - i] + 1 + maxSubset[k + 1][j - 1 - (k + 1)];
                maxCut[i][j - i] = 'I';  // In
            } else {
                maxSubset[i][j - i] = maxSubset[i][j - 1 - i];
                maxCut[i][j - i] = 'N';  // None
            }
        }
    }
#ifdef DEBUG
    cout << "cal done" << endl;
#endif
    findChordBottomUp(cout, maxCut, maxSubset, chords, result, smallestLen, N, 0, N - 1);
#endif
#ifdef DEBUG
    cout << "cal done" << endl;
#endif
    sort(result.begin(), result.end());
#ifdef DEBUG
    cout << "sort done" << endl;
#endif
    //////////// write the output file ///////////
    fstream fout;
    fout.open(argv[2], ios::out);
#ifndef __RECURSIVE__
    fout << maxSubset[0][N - 1] << "\n";
#else
    fout << maxSubset[0][N - 1 - smallestLen] << "\n";
#endif
    for (int i = 0; i < result.size(); i++) {
        fout << result[i] << ' ' << chords[result[i]] << "\n";
    }
    fout.close();
    return 0;
}
