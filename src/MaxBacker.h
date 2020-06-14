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
    vector<int> S;
    vector<vector<double> > w;
    void maximumBack();
    friend ostream& operator<<(ostream& os,MaxBacker& mb);

public:
    vector<int> Smin;
    void updateMaxBack(const vector<vector<double> > &w);
    MaxBacker(const vector<vector<double> > &w);
    double getMaxBack(){ return Cutmin;}
    pair<vector<int>, vector<int> > getPartition();
};

#endif