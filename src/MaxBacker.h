#ifndef _MY_MAXBACKER_
#define _MY_MAXBACKER_
#include "DisjSet.h"
#include <vector>
#include <iostream>
#include <set>
#include <algorithm>
#include <random>
using namespace std;

class MaxBacker{
private:
    double Cutmin;
    void maximumBack(const vector<vector<double> > &wf, vector<int> &S1);
public:
    void updateMaxBack(const vector<vector<double> > &wf, vector<int> &S1);
    MaxBacker(const vector<vector<double> > &wf, vector<int> &S1);
    double getMaxBack(){ return Cutmin;}
};

#endif