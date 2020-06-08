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
    double minCut;
    vector<int> S1;
    vector<int> S2;
    vector<int> G;
    DisjSet* dSet;
    void MINIMUMCUT(vector<vector<double> > &w);
    void MINIMUMCUTPHASE(vector<vector<double> > &w);
    friend ostream& operator<<(ostream& os, const MinCutter& dt);

public:
    MinCutter(vector<vector<double> > &w);
    MinCutter(int n);
    ~MinCutter();
    void MINIMUMCUTUPDATE(vector<vector<double> > &w);
    double getMinCut(){ return minCut;}
    pair<vector<int>, vector<int> > getPartition(){return make_pair(S1,S2);}

private:
    int n;
    
};

#endif