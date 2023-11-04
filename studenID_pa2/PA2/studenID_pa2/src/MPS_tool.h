#ifndef _MPS_TOOL_
#define _MPS_TOOL_

struct Chord {
    Chord(int _start, int _end) : start(_start), end(_end){};
    int start;
    int end;
};

class Matrix2D {
   public:
    Matrix2D(){};
    Matrix2D(int _N, int _M) : N(_N), M(_M){};
    int N;
    int M;

   private:
    int **maxSubset;
};
#endif