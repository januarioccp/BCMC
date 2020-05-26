#ifndef NODE_INFO_H_
#define NODE_INFO_H_

#include <ilcplex/ilocplex.h>

/** A subclass of IloCplex::MIPCallbackI::NodeData that stores data related to the node **/
class NodeInfo : public IloCplex::MIPCallbackI::NodeData 
{

   //Depth of node at which this data is stored. 
   unsigned int const depth;

   //stores how many times the cutcallback have been called in the node for each of the 2 cuts
   unsigned int p3LoopCount;

   //stores the last iteration of the cut loop in which the node objective value was improved
   //This can be used, for example, to prevent the callback execution after a certain number of iterations without improvement.
   unsigned int lastGoodLoop;

   //objective function value at the node
   float obj_value;

   //determines wether or not the separating routine will run ih the node
   bool separating_p3;


  public: 

      //cplex does not store node data for the root node. It is then necessary to explicitily maintain the NodeInfo object for the root of the B&C tree.
      //This object is a static object of the class NodeInfo.
      static NodeInfo* rootData;
      //call the method below before the B&C start
      static void initRootData();
      
      NodeInfo(unsigned int idepth); 

      bool fixedVariables;

      unsigned int getDepth() const;
      unsigned int getLastGoodLoop() const;
    
      unsigned int getP3LoopCount() const;
       
      //increment the cut loop count
      void incP3LoopCount();
      
      void setLastGoodLoop(unsigned int iter);
      void setNodeObjValue(float obj_value);
      float getNodeObjValue() const;
      
      bool isSeparating();
      bool isSeparatingP3();
      void startSeparationP3();
      void stopSeparationP3();
      void stopSeparation();
};

#endif
