#ifndef _MY_DISJSET_
#define _MY_DISJSET_
#include<vector>
using namespace std;

class DisjSet;

class DisjSet { 
    vector<int> rank;
    vector<int> parent;
    int n; 
public: 
    // Constructor to create and initialize sets of n items 
    DisjSet(int n);
     
    // Creates n single item sets 
    void makeSet();
  
    // Finds set of given item x 
    int find(int x);
  
    // Do union of two sets represented 
    // by x and y. 
    void Union(int x, int y) ;
};
#endif
