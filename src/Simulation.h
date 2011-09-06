#ifndef SIMULATION_GUARD
#define SIMULATION_GUARD

#include "Timer.h"
#include <vector>
#include <iostream>
#include <cstdlib>

class Simulation {
 public:
  Simulation(size_t repeats, bool printPartialData):
   _repeats(repeats), _printPartialData(printPartialData), _simType("none") {}

  void makeStandard(size_t varCount, size_t inserts, size_t queries);

  template<class DivFinder>
  void run();
  template<class DivFinder, class Param1>
  void run(const Param1& param1);
  template<class DivFinder, class Param1, class Param2>
  void run(const Param1& param1, const Param2& param2);
  template<class DivFinder, class Param1, class Param2, class Param3>
  void run(const Param1& param1, const Param2& param2, const Param3& param3);

  void printData(std::ostream& out) const;

 private:
  struct SimData {
    bool operator<(const SimData& sd) const;
    void print(std::ostream& out);

    std::string _name;
    unsigned long _mseconds;
    unsigned long _expQueryCount;
  };

  template<class DivFinder>
  void run(DivFinder& finder);

  enum EventType {
    Insertion,
    QueryNoDivisor,
    QueryHasDivisor,
    QueryUnknown
  };
  struct Event {
    EventType _type;
	std::vector<int> _monomial;
  };

  std::vector<Event> _events;
  std::vector<SimData> _data;
  size_t _varCount;
  size_t _repeats;
  bool _printPartialData;
  std::string _simType;
};

template<class DivFinder>
void Simulation::run() {
  DivFinder finder(_varCount);
  run(finder);
}

template<class DivFinder, class Param1>
void Simulation::run(const Param1& param1) {
  DivFinder finder(_varCount, param1);
  run(finder);
}

template<class DivFinder, class Param1, class Param2>
void Simulation::run(const Param1& param1, const Param2& param2) {
  DivFinder finder(_varCount, param1, param2);
  run(finder);
}

template<class DivFinder, class Param1, class Param2, class Param3>
void Simulation::run
(const Param1& param1, const Param2& param2, const Param3& param3) {
  DivFinder finder(_varCount, param1, param2, param3);
  run(finder);
}

template<class DivFinder>
void Simulation::run(DivFinder& finder) {
  Timer timer;
  for (size_t step = 0; step < _repeats; ++step) {
	for (size_t i = 0; i < _events.size(); ++i) {
	  Event& e = _events[i];
	  if (e._type == Insertion)
		finder.insert(e._monomial);
	  else {
		typename DivFinder::const_iterator it = finder.findDivisor(e._monomial);
		if (it == finder.end()) {
		  if (e._type == QueryHasDivisor) {
			std::cerr << "Divisor finder \"" << finder.getName()
					  << "\" failed to find divisor." << std::endl;
			std::exit(1);
		  }
		  e._type = QueryNoDivisor;
		} else {
#ifdef DEBUG
		  for (size_t var = 0; var < _varCount; ++var) {
			ASSERT((*it)[var] <= e._monomial[var]);
		  }
#endif
		  if (e._type == QueryNoDivisor) {
			std::cerr << "Divisor finder \"" << finder.getName() <<
			  "\" found incorrect divisor." << std::endl;
			std::exit(1);
		  }
		  e._type = QueryHasDivisor;
		}
	  }
	}
  }

  SimData data;
  data._mseconds = (unsigned long)timer.getMilliseconds();
  data._name = finder.getName();
  data._expQueryCount = finder.getExpQueryCount();
  _data.push_back(data);
  if (_printPartialData)
	data.print(std::cerr);
}

#endif
