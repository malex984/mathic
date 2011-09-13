#ifndef SIMULATION_GUARD
#define SIMULATION_GUARD

#include "Monomial.h"
#include "Timer.h"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <algorithm>

class Simulation {
 public:
  Simulation(size_t repeats, bool printPartialData):
   _repeats(repeats), _printPartialData(printPartialData), _simType("none") {}

  void makeStandard(size_t varCount, size_t inserts, size_t queries, bool findAll);

  template<class DivFinder>
  void run();
  template<class DivFinder, class Param1>
  void run(const Param1& param1);
  template<class DivFinder, class Param1, class Param2>
  void run(const Param1& param1, const Param2& param2);
  template<class DivFinder, class Param1, class Param2, class Param3>
  void run(const Param1& param1, const Param2& param2, const Param3& param3);
  template<class DivFinder, class P1, class P2, class P3, class P4>
  void run(const P1& p1, const P2& p2, const P3& param3, const P4& param4);
  template<class DivFinder, class P1, class P2, class P3, class P4, class P5>
  void run(const P1& p1, const P2& p2, const P3& param3, const P4& param4,
    const P5& p5);
  template<class DivFinder, class P1, class P2, class P3, class P4, class P5,
  class P6>
  void run(const P1& p1, const P2& p2, const P3& param3, const P4& param4,
    const P5& p5, const P6& p6);

  void printData(std::ostream& out) const;

 private:
  struct SimData {
    bool operator<(const SimData& sd) const;
    void print(std::ostream& out);

    std::string _name;
    unsigned long _mseconds;
    unsigned long long _expQueryCount;
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
#ifdef DEBUG
    std::vector<Monomial> _allDivisors;
#else
    size_t _divisorCount;
#endif
  };
  class DivisorStore;

  bool _findAll;
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

template<class DivFinder, class P1, class P2, class P3, class P4>
void Simulation::run
(const P1& param1, const P2& param2, const P3& param3, const P4& param4) {
  DivFinder finder(_varCount, param1, param2, param3, param4);
  run(finder);
}

template<class DivFinder, class P1, class P2, class P3, class P4, class P5>
void Simulation::run
(const P1& p1, const P2& p2, const P3& p3, const P4& p4, const P5& p5) {
  DivFinder finder(_varCount, p1, p2, p3, p4, p5);
  run(finder);
}

template<class DivFinder, class P1, class P2, class P3, class P4,
class P5, class P6>
void Simulation::run
(const P1& p1, const P2& p2, const P3& p3, const P4& p4, const P5& p5,
const P6& p6) {
  DivFinder finder(_varCount, p1, p2, p3, p4, p5, p6);
  run(finder);
}

class Simulation::DivisorStore {
public:
  DivisorStore() {clear();}
  void clear() {
#ifdef DEBUG
    _divisors.clear();
#else
    _divisorCount = 0;
#endif
  }

  void push_back(const Monomial& divisor) {
#ifdef DEBUG
    _divisors.push_back(divisor);
#else
    ++_divisorCount;
#endif
  }

  template<class Finder>
  void check(Event& e, const Finder& finder ) {
#ifdef DEBUG
    for (size_t d = 0; d < _divisors.size(); ++d) {
      for (size_t var = 0; var < e._monomial.size(); ++var) {
        ASSERT(_divisors[d][var] <= e._monomial[var]);
      }
    }
    std::sort(_divisors.begin(), _divisors.end());
#endif

    if (e._type == QueryUnknown) {          
      bool noDivisors;
#ifdef DEBUG
      e._allDivisors.clear();
      for (size_t i = 0; i < _divisors.size(); ++i)
        e._allDivisors.push_back(_divisors[i]);
      noDivisors = _divisors.empty();
#else
      e._divisorCount = _divisorCount;
      noDivisors = _divisorCount == 0;
#endif
      e._type = noDivisors ? QueryNoDivisor : QueryHasDivisor;
    } else {
#ifdef DEBUG
      for (size_t i = 0; i < _divisors.size(); ++i)
        ASSERT(_divisors[i] == e._allDivisors[i]);
#else
      if (_divisorCount != e._divisorCount) {
        std::cerr << "Divisor finder \"" << finder.getName() <<
          "\" found incorrect number of divisors." << std::endl;
        std::exit(1);
      }
#endif
    }
  }

private:
#ifdef DEBUG
  std::vector<Monomial> _divisors;
#else
  size_t _divisorCount;
#endif
};

template<class DivFinder>
void Simulation::run(DivFinder& finder) {
  Timer timer;
  std::vector<Monomial> divisors;
  for (size_t step = 0; step < _repeats; ++step) {
	for (size_t i = 0; i < _events.size(); ++i) {
	  Event& e = _events[i];
	  if (e._type == Insertion)
		finder.insert(e._monomial);
	  else if (!_findAll) {
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
      } else {
        ASSERT(_findAll);
        divisors.clear();
        DivisorStore store;
        const_cast<const DivFinder&>(finder) // to test const interface
          .findAllDivisors(e._monomial, store);
        store.check(e, finder);
        continue;



#ifdef DEBUG
        for (size_t d = 0; d < divisors.size(); ++d) {
          for (size_t var = 0; var < _varCount; ++var) {
            ASSERT(divisors[d][var] <= e._monomial[var]);
          }
        }
        std::sort(divisors.begin(), divisors.end());
#endif
        
        if (e._type == QueryUnknown) {          
#ifdef DEBUG
          e._allDivisors = divisors;
#else
          e._divisorCount = divisors.size();
#endif
          e._type = divisors.empty() ? QueryNoDivisor : QueryHasDivisor;
        } else {
#ifdef DEBUG
          ASSERT(divisors == e._allDivisors);
#else
          if (divisors.size() != e._divisorCount) {
            std::cerr << "Divisor finder \"" << finder.getName() <<
              "\" found incorrect number of divisors." << std::endl;
            std::exit(1);
          }
#endif
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
  std::cout << finder.size() << std::endl;
}

#endif
