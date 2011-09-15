#include "stdinc.h"

#include "DivListModel.h"
#include "KDTreeModel.h"
#include "Simulation.h"
#include "Timer.h"
#include <iostream>

int main() {
  const size_t repeats = IF_DEBUG(true ? 1 :) 1;
  Simulation sim(repeats, true);
  Timer timer;
  std::cout << "Generating simulation. ";
  timer.print(std::cout);
  std::cout << std::endl;
#ifdef DEBUG
  sim.makeStandard(10, 200, 1000, false);
#else
  sim.makeStandard(10, 1000000, 0, false);
#endif
  
  /*
  for (int minimizeOnInsert = 1; minimizeOnInsert <= 1; ++minimizeOnInsert) {
	for (int order = 0; order <= 2; ++order) {
      bool moveDivisorToFront = (order == 1);
      bool sortOnInsert = (order == 2);
      sim.run<DivListModel<1> >
		(minimizeOnInsert, moveDivisorToFront, sortOnInsert);
      sim.run<DivListModel<0> >
		(minimizeOnInsert, moveDivisorToFront, sortOnInsert);
	}
  }
  //*/

  IF_DEBUG(sim.run<DivListModel<0> >(1, 0, 0));
  /*
  sim.run<DivListModel<0> >(1, 0, 0)
  sim.run<DivListModel<0> >(1, 1, 0);
  sim.run<DivListModel<0> >(1, 0, 1);
  sim.run<DivListModel<1> >(1, 0, 0);
  sim.run<DivListModel<1> >(1, 1, 0);
  sim.run<DivListModel<1> >(1, 0, 1);
  //*/

  //*
  //sim.run<DivListModel<0> >(1, 0, 0);
  for (int mini = 1; mini <= 1; ++mini) {
	for (int sortOnInsert = 0; sortOnInsert <= 0; ++sortOnInsert) {
  	  for (int useDivisorCache = 0; useDivisorCache <= 0; ++useDivisorCache) {
        IF_DEBUG(sim.run<KDTreeModel<1> >(10, mini,sortOnInsert,useDivisorCache, 0.01, 1));
        
        
        sim.run<KDTreeModel<0> >(10, mini,sortOnInsert,useDivisorCache, 0, 0);
        sim.run<KDTreeModel<1> >(10, mini,sortOnInsert,useDivisorCache, 0, 0);
        sim.run<KDTreeModel<1> >(10, mini,sortOnInsert,useDivisorCache, 3, 2000);
        /*sim.run<KDTreeModel<0> >(8, mini,sortOnInsert,useDivisorCache, 0, 0);
        sim.run<KDTreeModel<1> >(8, mini,sortOnInsert,useDivisorCache, 0, 0);
        sim.run<KDTreeModel<0> >(20, mini,sortOnInsert,useDivisorCache, 0, 0);
        sim.run<KDTreeModel<1> >(20, mini,sortOnInsert,useDivisorCache, 0, 0);
        sim.run<KDTreeModel<0> >(40, mini,sortOnInsert,useDivisorCache, 0, 0);
        sim.run<KDTreeModel<1> >(40, mini,sortOnInsert,useDivisorCache, 0, 0);
        sim.run<KDTreeModel<0> >(60, mini,sortOnInsert,useDivisorCache, 0, 0);
        sim.run<KDTreeModel<1> >(60, mini,sortOnInsert,useDivisorCache, 0, 0);*/

        /*
        sim.run<KDTreeModel<1> >(10, mini,sortOnInsert,useDivisorCache, 0.001, 2000);
        sim.run<KDTreeModel<1> >(10, mini,sortOnInsert,useDivisorCache, 0.001, 4000);

        sim.run<KDTreeModel<0> >(10, mini,sortOnInsert,useDivisorCache, 0.001, 6000);
        sim.run<KDTreeModel<1> >(10, mini,sortOnInsert,useDivisorCache, 0.001, 6000);

        sim.run<KDTreeModel<1> >(10, mini,sortOnInsert,useDivisorCache, 0.001, 8000);
        sim.run<KDTreeModel<1> >(10, mini,sortOnInsert,useDivisorCache, 0.001, 10000);*/
        //sim.run<KDTreeModel>(10, mini,sortOnInsert,useDivisorCache, 4, 2000);
        //sim.run<KDTreeModel>(20, mini,sortOnInsert,useDivisorCache, 0.75, 2000);
        //sim.run<KDTreeModel>(8, mini,sortOnInsert,useDivisorCache, 0.75, 2000);
      }
	}
  }
  //*/

  std::cout << "\n\n";
  sim.printData(std::cout);
}
