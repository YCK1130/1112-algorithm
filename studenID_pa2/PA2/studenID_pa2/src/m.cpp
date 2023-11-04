#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#include "MPS_tool.h"
using namespace std;
void help_message() {
    cout << "usage: NTU_MPS <input_file> <output_file>" << endl;
    cout << "options:" << endl;
    cout << "   IS - Insersion Sort" << endl;
    cout << "   MS - Merge Sort" << endl;
    cout << "   QS - Quick Sort" << endl;
    cout << "   HS - Heap Sort" << endl;
}

template <class T>
void printMatrix(ostream& out, vector<vector<T>> Matrix, int N, int M) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            out << " " << Matrix[i][j];
        }
        out << "\n";
    }
    out << "\n";
    return;
}
void findChord(ostream& out, const vector<vector<char>>& MaxCut, const vector<vector<int>>& maxsub,
               const vector<int>& chord, vector<int>& result, int N, int i, int j) {
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
int calculateMaxSubset(vector<vector<char>>& maxCut, vector<vector<int>>& maxSubset,
                       const vector<int>& chords, int N, int i, int j) {
    if (maxCut[i][j] != 'X') return maxSubset[i][j];
    // case 1
    // both two chords containing the internal ends are outside of the internal
    if ((chords[i] < i || chords[i] > j) && (chords[j] < i || chords[j] > j)) {
        maxSubset[i][j] = calculateMaxSubset(maxCut, maxSubset, chords, N, i + 1, j - 1);
        maxCut[i][j] = 'N';  // None
    }
    // case 2
    // both two chords containing the internal ends are inside of the internal
    else if ((chords[i] > i && chords[i] < j) && (chords[j] > i && chords[j] < j)) {
        // if intercect
        if (chords[i] > chords[j]) {
            if (calculateMaxSubset(maxCut, maxSubset, chords, N, i, j - 1) >=
                calculateMaxSubset(maxCut, maxSubset, chords, N, i + 1, j)) {
                if (maxSubset[i][j - 1] >=
                    calculateMaxSubset(maxCut, maxSubset, chords, N, i + 1, j - 1)) {
                    maxSubset[i][j] = maxSubset[i][j - 1];
                    maxCut[i][j] = 'R';  // Right
                } else {
                    maxSubset[i][j] = maxSubset[i + 1][j - 1];
                    maxCut[i][j] = 'N';  // None
                }
            } else {
                if (maxSubset[i + 1][j] >= maxSubset[i + 1][j - 1]) {
                    maxSubset[i][j] = maxSubset[i + 1][j];
                    maxCut[i][j] = 'L';  // Left
                } else {
                    maxSubset[i][j] = maxSubset[i + 1][j - 1];
                    maxCut[i][j] = 'N';  // None
                }
            }
        } else {
            maxSubset[i][j] =
                calculateMaxSubset(maxCut, maxSubset, chords, N, chords[j], j) +
                calculateMaxSubset(maxCut, maxSubset, chords, N, i, chords[i]) +
                calculateMaxSubset(maxCut, maxSubset, chords, N, chords[i], chords[j]);
            maxCut[i][j] = 'B';  // Both
        }
    }
    // case 3
    // the internal ends forms a chord
    else if (chords[i] == j) {
        maxSubset[i][j] = calculateMaxSubset(maxCut, maxSubset, chords, N, i + 1, j - 1) + 1;
        maxCut[i][j] = 'S';  // Self
    }
    // case 4
    // one of two chords containing the internal ends is inside of the internal
    else if (chords[i] > i && chords[i] < j) {
        if (calculateMaxSubset(maxCut, maxSubset, chords, N, i, j - 1) >=
            calculateMaxSubset(maxCut, maxSubset, chords, N, i + 1, j - 1)) {
            maxSubset[i][j] = maxSubset[i][j - 1];
            maxCut[i][j] = 'R';  // Right
        } else {
            maxSubset[i][j] = maxSubset[i + 1][j - 1];
            maxCut[i][j] = 'N';  // None
        }

    } else if (chords[j] >= i && chords[j] <= j) {
        if (calculateMaxSubset(maxCut, maxSubset, chords, N, i + 1, j) >=
            calculateMaxSubset(maxCut, maxSubset, chords, N, i + 1, j - 1)) {
            maxSubset[i][j] = maxSubset[i + 1][j];
            maxCut[i][j] = 'L';  // Left
        } else {
            maxSubset[i][j] = maxSubset[i + 1][j - 1];
            maxCut[i][j] = 'N';  // None
        }
    } else {
        cerr << "something went wrong in i: " << i << ", j: " << j << endl;
        printMatrix(cerr, maxSubset, N, N);
    }

    return maxSubset[i][j];
    // cerr << "i: " << i << ", j: " << j << endl;
    // printMatrix(cout, maxSubset, N, N);
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
    int smallestLen = N;
    for (size_t i = 0; i < N / 2; i++) {
        fin >> end1 >> end2;
        // Chord chord(end1, end2);
        chords[end1] = end2;  // chords[end1] chords[end2] will form a chord
        chords[end2] = end1;
        if (abs(end1 - end2) < smallestLen) smallestLen = abs(end1 - end2);
    }
    fin.close();

    vector<vector<int>> maxSubset(N, vector<int>(N, 0));
    vector<vector<char>> maxCut(N, vector<char>(N, 'X'));
    ////////////     main part    ////////////////
    for (int cut = smallestLen; cut < verticeNum; cut++) {
        for (int i = 0; i < verticeNum - cut; i++) {
            // case 1
            // both two chords containing the internal ends are outside of the internal
            if ((chords[i] < i || chords[i] > i + cut) &&
                (chords[i + cut] < i || chords[i + cut] > i + cut)) {
                maxSubset[i][i + cut] = maxSubset[i + 1][i + cut - 1];
                maxCut[i][i + cut] = 'N';  // None
            }
            // case 2
            // both two chords containing the internal ends are inside of the internal
            else if ((chords[i] > i && chords[i] < i + cut) &&
                     (chords[i + cut] > i && chords[i + cut] < i + cut)) {
                // if intercect
                if (chords[i] > chords[i + cut]) {
                    if (maxSubset[i][i + cut - 1] >= maxSubset[i + 1][i + cut]) {
                        if (maxSubset[i][i + cut - 1] >= maxSubset[i + 1][i + cut - 1]) {
                            maxSubset[i][i + cut] = maxSubset[i][i + cut - 1];
                            maxCut[i][i + cut] = 'R';  // Right
                        } else {
                            maxSubset[i][i + cut] = maxSubset[i + 1][i + cut - 1];
                            maxCut[i][i + cut] = 'N';  // None
                        }
                    } else {
                        if (maxSubset[i + 1][i + cut] >= maxSubset[i + 1][i + cut - 1]) {
                            maxSubset[i][i + cut] = maxSubset[i + 1][i + cut];
                            maxCut[i][i + cut] = 'L';  // Left
                        } else {
                            maxSubset[i][i + cut] = maxSubset[i + 1][i + cut - 1];
                            maxCut[i][i + cut] = 'N';  // None
                        }
                    }
                } else {
                    maxSubset[i][i + cut] = maxSubset[i][chords[i]] +
                                            maxSubset[chords[i + cut]][i + cut] +
                                            maxSubset[chords[i]][chords[i + cut]];
                    maxCut[i][i + cut] = 'B';  // Both
                }
            }
            // case 3
            // the internal ends forms a chord
            else if (chords[i] == i + cut) {
                maxSubset[i][i + cut] = maxSubset[i + 1][i + cut - 1] + 1;
                maxCut[i][i + cut] = 'S';  // Self
            }
            // case 4
            // one of two chords containing the internal ends is inside of the internal
            else if (chords[i] > i && chords[i] < i + cut) {
                if (maxSubset[i][i + cut - 1] >= maxSubset[i + 1][i + cut - 1]) {
                    maxSubset[i][i + cut] = maxSubset[i][i + cut - 1];
                    maxCut[i][i + cut] = 'R';  // Right
                } else {
                    maxSubset[i][i + cut] = maxSubset[i + 1][i + cut - 1];
                    maxCut[i][i + cut] = 'N';  // None
                }
            } else if (chords[i + cut] >= i && chords[i + cut] <= i + cut) {
                if (maxSubset[i + 1][i + cut] >= maxSubset[i + 1][i + cut - 1]) {
                    maxSubset[i][i + cut] = maxSubset[i + 1][i + cut];
                    maxCut[i][i + cut] = 'L';  // Left
                } else {
                    maxSubset[i][i + cut] = maxSubset[i + 1][i + cut - 1];
                    maxCut[i][i + cut] = 'N';  // None
                }
            } else {
                cerr << "something went wrong in i: " << i << ", j: " << i + cut << endl;
                printMatrix(cerr, maxSubset, N, N);
            }
            // cerr << "i: " << i << ", j: " << i + cut << endl;
            // printMatrix(cout, maxSubset, N, N);
            // printMatrix(cout, maxCut, N, N);
        }
    }
    calculateMaxSubset(maxCut, maxSubset, chords, N, 0, N - 1);
    vector<int> result;
    // result.resize(maxSubset[0][N - 1]);
    cout << maxSubset[0][N - 1] << "\n";

    findChord(cout, maxCut, maxSubset, chords, result, N, 0, N - 1);
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
