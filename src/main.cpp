#include "stdinc.h"

#include "DivListModel.h"
#include "KDTreeModel.h"
#include "Simulation.h"
#include <iostream>

int main() {
  const size_t repeats = IF_DEBUG(true ? 1 :) 1;
  Simulation sim(repeats, true);
  sim.makeStandard(10, 1000000, 0);
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
  */
  /*
  sim.run<DivListModel<0> >(1, 0, 0);
  sim.run<DivListModel<0> >(1, 1, 0);
  sim.run<DivListModel<0> >(1, 0, 1);
  sim.run<DivListModel<1> >(1, 0, 0);
  sim.run<DivListModel<1> >(1, 1, 0);
  sim.run<DivListModel<1> >(1, 0, 1);
  */

  //*
  for (int minimizeOnInsert = 1; minimizeOnInsert <= 1; ++minimizeOnInsert) {
	for (int sortOnInsert = 1; sortOnInsert <= 1; ++sortOnInsert) {
      //sim.run<KDTreeModel>(100000, minimizeOnInsert, sortOnInsert);
      //sim.run<KDTreeModel>(400, minimizeOnInsert, sortOnInsert);
      //sim.run<KDTreeModel>(200, minimizeOnInsert, sortOnInsert);
      //sim.run<KDTreeModel>(100, minimizeOnInsert, sortOnInsert);
      sim.run<KDTreeModel>(50, minimizeOnInsert, sortOnInsert);
	}
  }
  //*/

  std::cout << "\n\n";
  sim.printData(std::cout);
}
