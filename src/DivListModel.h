#ifndef DIV_ARRAY_MODEL_GUARD
#define DIV_ARRAY_MODEL_GUARD

#include "DivList.h"
#include <string>
#include <vector>

/** Helper class for DivListModel. */
class DivListModelConfiguration {
 public:
  typedef int Exponent;
  typedef std::vector<Exponent> Monomial;
  typedef Monomial Entry;

  DivListModelConfiguration(size_t varCount, bool sortOnInsert):
  _varCount(varCount), _sortOnInsert(sortOnInsert), _expQueryCount(0) {}

  size_t getVarCount() const {return _varCount;}
  bool getSortOnInsert() const {return _sortOnInsert;}

  Exponent getExponent(const Monomial& monomial, size_t var) const {
    ++_expQueryCount;
    ASSERT(var < monomial.size());
    return monomial[var];
  }

  bool divides(const Monomial& a, const Monomial& b) const {
    for (size_t var = 0; var < getVarCount(); ++var)
      if (getExponent(b, var) < getExponent(a, var))
        return false;
    return true;
  }

  bool isLessThan(const Monomial& a, const Monomial& b) const {
    for (size_t var = 0; var < getVarCount(); ++var) {
      if (getExponent(a, var) < getExponent(b, var))
        return true;
      if (getExponent(b, var) < getExponent(a, var))
        return false;
	}
    return false;
  }

  unsigned long getExpQueryCount() const {return _expQueryCount;}

 private:
  const size_t _varCount;
  const bool _sortOnInsert;
  mutable unsigned long _expQueryCount;
};

template<bool UseLinkedList>
class DivListModel;

/** An instantiation of the capabilities of DivList. */
template<bool ULL>
class DivListModel {
 private:
  typedef DivListModelConfiguration C;
  typedef DivList<C, ULL> Finder;
 public:
  typedef typename Finder::iterator iterator;
  typedef typename Finder::const_iterator const_iterator;
  typedef typename Finder::Monomial Monomial;
  typedef typename Finder::Entry Entry;

 DivListModel(size_t varCount,
			  bool minimizeOnInsert,
			  bool moveDivisorToFront,
			  bool sortOnInsert):
  _finder(C(varCount, sortOnInsert)),
  _minimizeOnInsert(minimizeOnInsert),
  _moveDivisorToFront(moveDivisorToFront) {
    ASSERT(!sortOnInsert || !moveDivisorToFront);
  }

  void insert(const Entry& entry);
  iterator findDivisor(const Monomial& monomial) {
    iterator it = _finder.findDivisor(monomial);
    if (_moveDivisorToFront && it != _finder.end()) {
      _finder.moveToFront(it);
	  it = _finder.begin();
	}
    return it;
  }
  const_iterator findDivisor(const Monomial& monomial) const {
    return const_cast<DivListModel<ULL>&>(*this).findDivisor(monomial);
  }
  std::string getName() const;

  iterator begin() {return _finder.begin();}
  const_iterator begin() const {return _finder.begin();}
  iterator end() {return _finder.end();}
  const_iterator end() const {return _finder.end();}

  unsigned long getExpQueryCount() const {
    return _finder.getConfiguration().getExpQueryCount();
  }

 private:
  Finder _finder;
  bool _minimizeOnInsert;
  bool _moveDivisorToFront;
};

template<bool ULL>
inline void DivListModel<ULL>::insert(const Entry& entry) {
  if (!_minimizeOnInsert) {
    _finder.insert(entry);
    return;
  }
  if (findDivisor(entry) != _finder.end())
    return;
  bool hasMultiples = _finder.removeMultiples(entry);
  _finder.insert(entry);
  if (_moveDivisorToFront && hasMultiples) {
    iterator it = _finder.end();
    _finder.moveToFront(--it);
  }
} 

template<bool ULL>
inline std::string DivListModel<ULL>::getName() const {
  return _finder.getName() +
    (_minimizeOnInsert ? " remin" : " nomin") +
    (_moveDivisorToFront ? " toFront" : "");
}

#endif
