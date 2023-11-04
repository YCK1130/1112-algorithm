#ifndef _MPS_TOOL_
#define _MPS_TOOL_

struct Chord {
    Chord(int _start, int _end) : start(_start), end(_end){};
    int start;
    int end;
};

class int_Matrix2D {
   public:
    int_Matrix2D(){};
    int_Matrix2D(int _N, int _M) : N(_N), M(_M){};
    int N;
    int M;

   private:
    int **arr;
};
#endif