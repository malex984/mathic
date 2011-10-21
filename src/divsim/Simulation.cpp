#include "stdinc.h"
#include "Simulation.h"

#include "ColumnPrinter.h"
#include <cstdlib>
#include <algorithm>

namespace {
  std::string commafy(unsigned long long l) {
    std::stringstream out;
    out << l;
    std::string str;
    for (size_t i = 0; i < out.str().size(); ++i) {
      str += out.str()[i];
      if (i != out.str().size() - 1 && ((out.str().size() - i) % 3) == 1)
        str += ',';
    }
    return str;
  }

  void makeRandom(std::vector<int>& monomial) {
    for (size_t var = 0; var < monomial.size(); ++var)
      monomial[var] = rand() % 1000;
  }
}

void Simulation::makeStandard
  (size_t varCount, size_t inserts, size_t queries, bool findAll) {
  srand(0);

  _findAll = findAll;
  _varCount = varCount;
  _events.clear();
  for (size_t i = 0; i < inserts + queries; ++i) {
    /*Event event2;
    event2._type = StateUnknown;
    _events.push_back(event2);*/

    Event event;
    event._monomial.resize(varCount);
    makeRandom(event._monomial);
    event._type = (i <= inserts ? InsertUnknown : QueryUnknown);
    _events.push_back(event);
  }
}

void Simulation::printData(std::ostream& out) const {
  std::vector<SimData> sorted(_data);
  std::sort(sorted.begin(), sorted.end());
  out << "*** Simulation outcome for "
      << _repeats << " repeats ***" << std::endl;
  mic::ColumnPrinter pr;
  pr.addColumn(true);
  pr.addColumn(false, " ", "ms");
  pr.addColumn(false, " ", "eqs");
  for (std::vector<SimData>::const_iterator it = sorted.begin();
    it != sorted.end(); ++it) {
    pr[0] << it->_name << '\n';
    pr[1] << commafy(it->_mseconds) << '\n';
    pr[2] << commafy(it->_expQueryCount) << '\n';
  }
  pr.print(out);
}

void Simulation::SimData::print(std::ostream& out) {
  out << _name
    << " " << commafy(_mseconds) << " ms"
    << " " << commafy(_expQueryCount) << " eqs"
    << '\n';
}

bool Simulation::SimData::operator<(const SimData& sd) const {
  return _mseconds < sd._mseconds;
}
