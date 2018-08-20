#ifndef PUNCTURE_H_
#define PUNCTURE_H_

//#ifndef DEBUG
//#define DEBUG

using namespace std;
using namespace arma;

struct IncGenerator {
    int current;
    IncGenerator (int start) : current(start) {}
    int operator() () {return current++;}
};

vector<int> PunctureCode(umat H);
void combinations(vector<int> INDEX,int offset, int k, vector<int> comb, int must_include);
void vector_print(vector<int> vectorToPrint);
void vector_print(vector< vector<int> > vectorToPrint);


//#endif /* DEBUG */
#endif /* PUNCTURE_H_ */
