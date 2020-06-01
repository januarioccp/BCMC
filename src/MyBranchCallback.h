#ifndef MY_BRANCH_CALLBACK_H
#define MY_BRANCH_CALLBACK_H

#include <ilcplex/ilocplex.h>
#include "Util.h"

/** Branch callback that keeps track of node depth. * The callback creates exactly the same branches as CPLEX would do. * Additionally, it uses the NodeData mechanism to store at each node * its depth in the branch and bound tree. */ 
class MyBranchCallback : public IloCplex::BranchCallbackI 
{
   public:

      MyBranchCallback(IloEnv env);
      
      IloCplex::CallbackI *duplicateCallback() const; 
      
      void main(); 

      static double getRootGap();

      static struct timeval tempo_start_root;

      static double getTimeAtRoot();

      static bool hasBranched();

   private:

      static double root_gap;

      static double time_at_root;

      static bool branched;

//      IloBoolVarArray y;

//      bool y_branching;

};

#endif
