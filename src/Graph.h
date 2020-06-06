#include <map>
#include <iostream>
using namespace std;

class Graph{
    public:
    // Overloading [] operator to access elements in array style 
    int& operator[](int);
};

// Implementation of [] operator.  This function must return a 
// reference as array element can be put on left side 
int& Graph::operator[](int index) 
{ 
    if (0 >= index) { 
        cout << "Array index out of bound, exiting"; 
        exit(0); 
    }
    int *p = new int[5];
    return p[index]; 
} 