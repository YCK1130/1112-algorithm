// **************************************************************************
//  File       [sort_tool.cpp]
//  Author     [Yu-Hao Ho]
//  Synopsis   [The implementation of the SortTool Class]
//  Modify     [2020/9/15 Mu-Ting Wu]
// **************************************************************************

#include "sort_tool.h"

#include <iostream>

#define INT_MAX 2147483647
// Constructor
SortTool::SortTool() {}

// Insertsion sort method
void SortTool::InsertionSort(vector<int>& data) {
    // Function : Insertion sort
    // TODO : Please complete insertion sort code here

    int key;
    for (int i = 1; i < data.size(); i++) {
        key = data[i];
        int j = i - 1;
        while (j >= 0 && data[j] > key) {
            data[j + 1] = data[j];
            j--;
        }
        data[j + 1] = key;
    }
}

// Quick sort method
void SortTool::QuickSort(vector<int>& data) { QuickSortSubVector(data, 0, data.size() - 1); }
// Sort subvector (Quick sort)
void SortTool::QuickSortSubVector(vector<int>& data, int low, int high) {
    // Function : Quick sort subvector
    // TODO : Please complete QuickSortSubVector code here
    // Hint : recursively call itself
    //        Partition function is needed
    if (low >= high) return;
    int q = Partition(data, low, high);
    QuickSortSubVector(data, low, q - 1);
    QuickSortSubVector(data, q + 1, high);
}

int SortTool::Partition(vector<int>& data, int low, int high) {
    // Function : Partition the vector
    // TODO : Please complete the function
    // Hint : Textbook page 171
    int pivot = high;
    int lowest, highest;
    lowest = data[low] > data[high] ? high : low;
    highest = data[low] <= data[high] ? high : low;

    pivot = data[lowest] > data[(low + high) / 2] ? lowest : (low + high) / 2;
    pivot = data[highest] > data[pivot] ? pivot : highest;

    int pivot_value = data[pivot];
    Swap(data, high, pivot);
    pivot = high;
    int i = low - 1;
    for (int j = low; j <= high; j++) {
        if (j == pivot) continue;
        if (data[j] <= pivot_value) {
            i++;
            Swap(data, i, j);
        }
    }
    Swap(data, i + 1, pivot);
    return i + 1;
}
void SortTool::Swap(vector<int>& data, int n, int m) {
    int tmp = data[n];
    data[n] = data[m];
    data[m] = tmp;
}
// Merge sort method
void SortTool::MergeSort(vector<int>& data) { MergeSortSubVector(data, 0, data.size() - 1); }

// Sort subvector (Merge sort)
void SortTool::MergeSortSubVector(vector<int>& data, int low, int high) {
    // Function : Merge sort subvector
    // TODO : Please complete MergeSortSubVector code here
    // Hint : recursively call itself
    //        Merge function is needed
    if (low >= high) return;
    int q = (low + high) / 2;
    MergeSortSubVector(data, low, q);
    MergeSortSubVector(data, q + 1, high);
    Merge(data, low, q, q + 1, high);
}

// Merge
void SortTool::Merge(vector<int>& data, int low, int middle1, int middle2, int high) {
    // Function : Merge two sorted subvector
    // TODO : Please complete the function
    vector<int> L(data.begin() + low, data.begin() + middle1 + 1);
    vector<int> R(data.begin() + middle2, data.begin() + high + 1);
    L.push_back(INT_MAX);
    R.push_back(INT_MAX);
    int i = 0, j = 0;
    for (int k = low; k < high + 1; k++) {
        if (L[i] <= R[j]) {
            data[k] = L[i];
            i++;
        } else {
            data[k] = R[j];
            j++;
        }
    }
}

// Heap sort method
int HeapParentIndex(int myIndex) { return (myIndex - 1) / 2; }
int HeapLeftIndex(int myIndex) { return myIndex * 2 + 1; }
int HeapRightIndex(int myIndex) { return myIndex * 2 + 2; }
void SortTool::HeapSort(vector<int>& data) {
    // Build Max-Heap
    BuildMaxHeap(data);
    // 1. Swap data[0] which is max value and data[i] so that the max value will be in correct
    // location
    // 2. Do max-heapify for data[0]
    for (int i = data.size() - 1; i >= 1; i--) {
        swap(data[0], data[i]);
        heapSize--;
        MaxHeapify(data, 0);
    }
}

// Max heapify
void SortTool::MaxHeapify(vector<int>& data, int root) {
    // Function : Make tree with given root be a max-heap if both right and left sub-tree are
    // max-heap
    // TODO : Please complete max-heapify code here
    int L_Index = HeapLeftIndex(root);
    int R_Index = HeapRightIndex(root);
    int largest;
    if (L_Index < heapSize && data[L_Index] > data[root])
        largest = L_Index;
    else
        largest = root;
    if (R_Index < heapSize && data[R_Index] > data[largest]) largest = R_Index;
    if (largest != root) {
        Swap(data, root, largest);
        MaxHeapify(data, largest);
    }
}

// Build max heap
void SortTool::BuildMaxHeap(vector<int>& data) {
    heapSize = data.size();  // initialize heap size
    // Function : Make input data become a max-heap
    // TODO : Please complete BuildMaxHeap code here
    for (int i = data.size(); i >= 0; i--) {
        MaxHeapify(data, i);
    }
}
