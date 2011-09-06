#include "stdinc.h"

#include "DivListModel.h"
#include "Simulation.h"
#include <iostream>

int main() {
  const size_t repeats = IF_DEBUG(true ? 1 :) 10;
  Simulation sim(repeats, true);
  sim.makeStandard(10, 1000, 100000);

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

  std::cout << "\n\n";
  sim.printData(std::cout);
}
