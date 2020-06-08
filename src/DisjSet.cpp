#include "DisjSet.h"

DisjSet::DisjSet(int n)
{
    rank.resize(n);
    parent.resize(n);
    this->n = n;
    makeSet();
}

void DisjSet::makeSet(){
    for (int i = 0; i < n; i++)
        parent[i] = i;
}

int DisjSet::find(int x) 
    { 
        // Finds the representative of the set 
        // that x is an element of 
        if (parent[x] != x) { 
  
            // if x is not the parent of itself 
            // Then x is not the representative of 
            // his set, 
            parent[x] = find(parent[x]); 
  
            // so we recursively call Find on its parent 
            // and move i's node directly under the 
            // representative of this set 
        } 
  
        return parent[x]; 
    } 

void DisjSet::Union(int x, int y) 
    { 
        // Find current sets of x and y 
        int xset = find(x); 
        int yset = find(y); 
  
        // If they are already in same set 
        if (xset == yset) 
            return; 
  
        // Put smaller ranked item under 
        // bigger ranked item if ranks are 
        // different 
        if (rank[xset] < rank[yset]) { 
            parent[xset] = yset; 
        } 
        else if (rank[xset] > rank[yset]) { 
            parent[yset] = xset; 
        } 
  
        // If ranks are same, then increment 
        // rank. 
        else { 
            parent[yset] = xset; 
            rank[xset] = rank[xset] + 1; 
        } 
    } 