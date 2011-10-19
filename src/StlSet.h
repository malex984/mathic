#ifndef MATHIC_STL_SET_GUARD
#define MATHIC_STL_SET_GUARD

#include "stdinc.h"
#include <set>
#include <string>
#include <ostream>

namespace mathic {
  template<class C>
	class StlSet {
  public:
	typedef C Configuration;
	typedef typename Configuration::Entry Entry;

  StlSet(const Configuration& configuration):
    _conf(configuration), _set(Cmp(&_conf)) {}
	Configuration& getConfiguration() {return _conf;}
	const Configuration& getConfiguration() const {return _conf;}

	std::string getName() const {return "stlset";}
	void push(Entry entry) {_set.insert(entry);}
	void push(const Entry* begin, const Entry* end) {
	  for (; begin != end; ++begin)
		push(*begin);
	}
	Entry pop() {
	  Entry top = *_set.begin();
	  _set.erase(_set.begin());
	  return top;
	}
	Entry top() {
	  return *_set.begin();
	}
	bool empty() const {return _set.empty();}
	void print(std::ostream& out) const {
	  out << getName() << ":\n";
	  for (typename std::multiset<Entry, Cmp>::const_iterator it = _set.begin();
		   it != _set.end(); ++it)
		out << ' ' << *it;
	  out << '\n';
	}

  private:
	Configuration _conf;
	struct Cmp {
    Cmp(const Configuration* conf): _conf(conf) {}
	  bool operator()(const Entry& a, const Entry& b) {
		return _conf->cmpLessThan(_conf->compare(b, a));
	  }
	  const Configuration* _conf; // should be &, but gcc complains
	};
	std::multiset<Entry, Cmp> _set;
  };
}

#endif
