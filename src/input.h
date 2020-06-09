#ifndef INPUT_H_INCLUDED
#define INPUT_H_INCLUDED
#include <iostream>
using namespace std;

class Input{
    public:
        Input(){};
        Input( int , char**);

        // Dimension
        unsigned dimensionGet();
        void dimensionSet(unsigned d);

        // Distance
        double distanceGet(unsigned i, unsigned j) const;
        void distanceSet(unsigned i, unsigned j, double value);
        
        friend ostream & operator << (ostream &out, Input &c);

    private:
        unsigned dimension_; // quantidade total de vertices
        double ** distance_; // matriz de adjacencia

        
};

#endif //INPUT_H_INCLUDED