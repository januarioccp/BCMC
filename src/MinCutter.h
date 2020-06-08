#ifndef _MY_MINCUTTER_
#define _MY_MINCUTTER_
#include "DisjSet.h"
#include <vector>
#include <iostream>
#include <set>
#include <algorithm>
#include <random>
using namespace std;

class MinCutter{
private:
    double minCut;
    vector<int> S1;
    vector<int> S2;
    vector<int> G;
    DisjSet* dSet;
    vector<vector<double> > w;
    void MINIMUMCUT();
    int MINIMUMCUTPHASE();
    friend ostream& operator<<(ostream& os, const MinCutter& dt);

public:
    MinCutter(const vector<vector<double> > &w);
    MinCutter(const vector<vector<int> > &w);
    ~MinCutter();
    double getMinCut(){ return minCut;}
    pair<vector<int>, vector<int> > getPartition(){return make_pair(S1,S2);}
};

#endif