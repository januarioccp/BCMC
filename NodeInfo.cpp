#include "NodeInfo.h"

//metodos da classe DepthData.
//Um objeto desta classe é associado a cada nó do B&C duranto o branching. Esse objeto é recuperado no cut callback.
NodeInfo::NodeInfo(unsigned  int depth) : depth(depth)
{
   p3LoopCount = lastGoodLoop = obj_value = 0; 
   separating_p3 = true;
} 

void NodeInfo::initRootData()
{
   rootData = new NodeInfo(0);
}

unsigned int NodeInfo::getDepth() const 
{ 
   return this->depth; 
}

unsigned int NodeInfo::getLastGoodLoop() const
{
   return this->lastGoodLoop;
}

unsigned int NodeInfo::getP3LoopCount() const
{
   return this->p3LoopCount;
}

float NodeInfo::getNodeObjValue() const
{
   return this->obj_value;
}

void NodeInfo::incP3LoopCount()
{
   this->p3LoopCount++;
}

void NodeInfo::setLastGoodLoop(unsigned int iter)
{
   this->lastGoodLoop = iter;
}

void NodeInfo::setNodeObjValue(float obj_value)
{
   this->obj_value = obj_value;
}

bool NodeInfo::isSeparating()
{
   return this->separating_p3;
}

bool NodeInfo::isSeparatingP3()
{
   return this->separating_p3;
}

void NodeInfo::stopSeparationP3()
{
   this->separating_p3 = false;
}

void NodeInfo::stopSeparation()
{
   this->separating_p3 = false;
}


void NodeInfo::startSeparationP3()
{
   this->separating_p3 = true;
}

//ponteiro estático para o objeto NodeInfo que contém os dados da raiz
NodeInfo* NodeInfo::rootData = NULL;
