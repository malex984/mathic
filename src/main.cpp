#include "stdinc.h"

#include "DivListModel.h"
#include "Simulation.h"
#include <iostream>

int main() {
  Simulation sim(2, true);
  sim.makeStandard(10, 1000, 100000);
  sim.run<DivListModel<0> >(0);
  sim.run<DivListModel<0> >(1);
  sim.run<DivListModel<1> >(1);
  sim.run<DivListModel<1> >(0);
  std::cout << "\n\n";
  sim.printData(std::cout);
}
