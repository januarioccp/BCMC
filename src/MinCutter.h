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
    vector<int> G;
    vector<int> A;
    DisjSet* partition;
    DisjSet* best_partition;
    int best_last;
    vector<vector<double> > w;
    void minimumCut(vector<vector<double> > &w, vector<int> &S);
    int MINIMUMCUTPHASE(vector<vector<double> > &w);

public:
    void updateMinCut(vector<vector<double> > &w, vector<int> &S);
    MinCutter(vector<vector<double> > &w, vector<int> &S);
    ~MinCutter();
    double getMinCut(){ return minCut;}
};

#endif