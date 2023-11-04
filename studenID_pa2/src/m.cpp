#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#include "MPS_tool.h"
using namespace std;

// #define __RECURSIVE__
// #define DEBUG
void help_message() {
    cout << "usage: NTU_MPS <input_file> <output_file>" << endl;
    cout << "options:" << endl;
    cout << "   IS - Insersion Sort" << endl;
    cout << "   MS - Merge Sort" << endl;
    cout << "   QS - Quick Sort" << endl;
    cout << "   HS - Heap Sort" << endl;
}

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
void findChord(ostream& out, const vector<vector<char> >& MaxCut,
               const vector<vector<int> >& maxsub, const vector<int>& chord, vector<int>& result,
               int N, int i, int j) {
    if (i >= j) return;
    if (maxsub[i][j] == 0) return;
    if (MaxCut[i][j] == 'B') {
        findChord(out, MaxCut, maxsub, chord, result, N, i, chord[i]);
        findChord(out, MaxCut, maxsub, chord, result, N, chord[j], j);
        findChord(out, MaxCut, maxsub, chord, result, N, chord[i], chord[j]);
        return;
    }
    if (MaxCut[i][j] == 'N') {
        findChord(out, MaxCut, maxsub, chord, result, N, i + 1, j - 1);
        return;
    }
    if (MaxCut[i][j] == 'S') {
        // out << i << ' ' << j << "\n";
        result.push_back(i);
        findChord(out, MaxCut, maxsub, chord, result, N, i + 1, j - 1);
        return;
    }
    if (MaxCut[i][j] == 'L') {
        findChord(out, MaxCut, maxsub, chord, result, N, i + 1, j);
        return;
    }
    if (MaxCut[i][j] == 'R') {
        findChord(out, MaxCut, maxsub, chord, result, N, i, j - 1);
        return;
    }
}
void findChord2(ostream& out, char** MaxCut, int** maxsub, const vector<int>& chords,
                vector<int>& result, int N, int i, int j) {
    if (i >= j) return;
    int realj = j - i;
    if (maxsub[i][realj] == 0) return;
    if (MaxCut[i][realj] == 'N') {
        findChord2(out, MaxCut, maxsub, chords, result, N, i + 1, j - 1);
        return;
    }
    if (MaxCut[i][realj] == 'S') {
        result.push_back(i);
        findChord2(out, MaxCut, maxsub, chords, result, N, i + 1, j - 1);
        return;
    }
    if (MaxCut[i][realj] == 'L') {
        findChord2(out, MaxCut, maxsub, chords, result, N, chords[j], j);
        findChord2(out, MaxCut, maxsub, chords, result, N, i, chords[j] - 1);
        return;
    }
    if (MaxCut[i][realj] == 'R') {
        findChord2(out, MaxCut, maxsub, chords, result, N, i, chords[i]);
        findChord2(out, MaxCut, maxsub, chords, result, N, chords[i] + 1, j);
        return;
    }
}
int calculateMaxSubset2(char** maxCut, int** maxSubset, const vector<int>& chords, const int& N,
                        const int& smallestLen, const int i, const int j) {
    if (i >= j || j - i < smallestLen) return 0;
    if (i < 0 || j >= N) return 0;
    int realj = j - i;
    if (maxCut[i][realj] != 'X') return maxSubset[i][realj];

#ifdef DEBUG
    cerr << "entering i: " << i << ", j: " << j << endl;
#endif
    // case 1
    // both two chords containing the internal ends are outside of the internal
    if ((chords[i] < i || chords[i] > j) && (chords[j] < i || chords[j] > j)) {
        maxSubset[i][realj] =
            calculateMaxSubset2(maxCut, maxSubset, chords, N, smallestLen, i + 1, j - 1);
        maxCut[i][realj] = 'N';  // None
    }
    // case 2
    // not forms a chord
    else if (chords[i] != j) {
        maxSubset[i][realj] =
            calculateMaxSubset2(maxCut, maxSubset, chords, N, smallestLen, i + 1, j - 1);
        maxCut[i][realj] = 'N';  // None
        if ((chords[i] < j && chords[i] > i) &&
            calculateMaxSubset2(maxCut, maxSubset, chords, N, smallestLen, i, chords[i]) +
                    calculateMaxSubset2(maxCut, maxSubset, chords, N, smallestLen, chords[i] + 1,
                                        j) >=
                maxSubset[i][realj]) {
            maxSubset[i][realj] =
                maxSubset[i][chords[i] - i] + maxSubset[chords[i] + 1][j - (chords[i] + 1)];
            maxCut[i][realj] = 'R';  // Right
        }
        if ((chords[j] < j && chords[j] > i) &&
            calculateMaxSubset2(maxCut, maxSubset, chords, N, smallestLen, chords[j], j) +
                    calculateMaxSubset2(maxCut, maxSubset, chords, N, smallestLen, i,
                                        chords[j] - 1) >=
                maxSubset[i][realj]) {
            maxSubset[i][realj] =
                maxSubset[chords[j]][j - chords[j]] + maxSubset[i][chords[j] - 1 - i];
            maxCut[i][realj] = 'L';  // Left
        }
    }
    // case 4
    // the internal ends forms a chord
    else if (chords[i] == j) {
        maxSubset[i][realj] =
            calculateMaxSubset2(maxCut, maxSubset, chords, N, smallestLen, i + 1, j - 1) + 1;
        maxCut[i][realj] = 'S';  // Self
    } else {
        cerr << "something went wrong in i: " << i << ", j: " << j << endl;
        printMatrix(cerr, maxSubset, N, N);
    }
#ifdef DEBUG
    cerr << maxCut[i][realj] << endl;
    printMatrix(cerr, maxSubset, N, N);
    cerr << "ending i: " << i << ", j: " << j << '\n' << endl;
#endif

    return maxSubset[i][realj];
    // printMatrix(cout, maxCut, N, N);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        help_message();
        return 0;
    }
    //////////// read the input file /////////////

    char buffer[200];
    fstream fin(argv[1]);
    // fin.getline(buffer, 200);
    // fin.getline(buffer, 200);

    int verticeNum;
    fin >> verticeNum;

    int end1, end2;
    vector<int> chords;
    int N = verticeNum;
    chords.resize(verticeNum);
    cout << "chords size " << chords.capacity() << endl;

    int smallestLen = N;
    for (size_t i = 0; i < N / 2; i++) {
        fin >> end1 >> end2;
        // Chord chord(end1, end2);
        chords[end1] = end2;  // chords[end1] chords[end2] will form a chord
        chords[end2] = end1;
        if (abs(end1 - end2) < smallestLen) smallestLen = abs(end1 - end2);
    }
    fin.close();

    int** maxSubset = new int*[N];
    char** maxCut = new char*[N];
    for (int i = 0; i < N; i++) {
        maxSubset[i] = new int[N - i];
        maxCut[i] = new char[N - i];
        for (int j = 0; j < N - i; j++) {
            maxSubset[i][j] = 0;
            maxCut[i][j] = 'X';
        }
    }
    // vector<vector<int> > maxSubset(N, vector<int>(N));
    // vector<vector<char> > maxCut(N, vector<char>(N));
    // for (int i = 0; i < N; i++) {
    //     for (int j = i; j < N; j++) {
    //         maxSubset[i][j] = 0;
    //         maxCut[i][j] = 'X';
    //     }
    // }
    // cout << "maxSubset size " << maxSubset.capacity() << endl;
    cout << "maxSubset allocate done" << endl;
#ifdef DEBUG
    cout << "init done" << endl;
#endif
    ////////////     main part    ////////////////
#ifndef __RECURSIVE__
    for (int cut = smallestLen; cut < verticeNum; cut++) {
        for (int i = 0; i < verticeNum - cut; i++) {
            int j = i + cut;
            // case 1
            // both two chords containing the internal ends are outside of the internal
            if ((chords[i] < i || chords[i] > j) && (chords[j] < i || chords[j] > j)) {
                maxSubset[i][j - i] = maxSubset[i + 1][j - 1 - i - 1];
                maxCut[i][j - i] = 'N';  // None
            }
            // case 2
            // not form a chord
            else if (chords[i] != j) {
                maxSubset[i][j - i] = maxSubset[i + 1][j - 1 - i - 1];
                maxCut[i][j - i] = 'N';  // None
                if ((chords[i] < j && chords[i] > i) &&
                    maxSubset[i][chords[i] - i] + maxSubset[chords[i] + 1][j - 1 - chords[i]] >=
                        maxSubset[i + 1][j - 1 - i - 1]) {
                    maxSubset[i][j - i] =
                        maxSubset[i][chords[i] - i] + maxSubset[chords[i] + 1][j - 1 - chords[i]];
                    maxCut[i][j - i] = 'R';  // Right
                }
                if ((chords[j] < j && chords[j] > i) &&
                    maxSubset[chords[j]][j - chords[j]] + maxSubset[i][chords[j] - 1 - i] >=
                        maxSubset[i][j - i]) {
                    maxSubset[i][j - i] =
                        maxSubset[chords[j]][j - chords[j]] + maxSubset[i][chords[j] - 1 - i];
                    maxCut[i][j - i] = 'L';  // Left
                }

            }
            // case 3
            // the internal ends forms a chord
            else if (chords[i] == j) {
                maxSubset[i][j - i] = maxSubset[i + 1][j - 1 - i - 1] + 1;
                maxCut[i][j - i] = 'S';  // Self
            } else {
                cerr << "something went wrong in i: " << i << ", j: " << j << endl;
                printMatrix(cerr, maxSubset, N, N);
            }
        }
    }
#endif
#ifdef __RECURSIVE__
    calculateMaxSubset2(maxCut, maxSubset, chords, N, smallestLen, 0, N - 1);
#endif
    vector<int> result;
    // result.resize(maxSubset[0][N - 1]);
    cout << maxSubset[0][N - 1] << "\n";

    findChord2(cout, maxCut, maxSubset, chords, result, N, 0, N - 1);
    // cerr << "result size: " << result.size() << endl;
    sort(result.begin(), result.end());
    // cout << "-------------\n";
    // for (int i = 0; i < result.size(); i++) {
    //     cout << result[i] << ' ' << chords[result[i]] << "\n";
    // }
    // cout << "-------------\n";
    // cerr << "result size: " << result.size() << endl;
    cout << "sort done" << endl;
    //////////// write the output file ///////////
    fstream fout;
    fout.open(argv[2], ios::out);
    // printMatrix(fout, maxSubset, N, N);
    // printMatrix(fout, maxCut, N, N);
    fout << maxSubset[0][N - 1] << "\n";
    for (int i = 0; i < result.size(); i++) {
        fout << result[i] << ' ' << chords[result[i]] << "\n";
    }
    // for (int i = 0; i < chords.size(); i++) fout << i << " " << chords[i] << endl;
    fout.close();
    return 0;
}
