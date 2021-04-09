/**
 * @file MdFlexSingleNode.h
 * @author J. Körner
 * @date 07.04.2021
 */

#include "MdFlexSingleNode.h"

MdFlexibleSingleNode::Initialize(int argc, char** argv){
#if defined(AUTOPAS_INTERNODE_TUNING)
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout << "rank: " << rank << std::endl;
#endif

}

MdFlexibleSingleNode::Run(){
}

MdFlexibleSingleNode::Finalize(){
#if defined(AUTOPAS_INTERNODE_TUNING)
  MPI_Finalize();
#endif
}
