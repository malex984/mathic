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

  DivListModelConfiguration(size_t varCount):
   _varCount(varCount), _expQueryCount(0) {}

  size_t getVarCount() const {return _varCount;}

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

  unsigned long getExpQueryCount() const {return _expQueryCount;}

 private:
  const size_t _varCount;
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

 DivListModel(size_t varCount, bool minimizeOnInsert):
  _finder(C(varCount)), _minimizeOnInsert(minimizeOnInsert) {}

  void insert(const Entry& entry);
  iterator findDivisor(const Monomial& monomial) {
    return _finder.findDivisor(monomial);
  }
  const_iterator findDivisor(const Monomial& monomial) const {
    return _finder.findDivisor(monomial);
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
};

template<bool ULL>
inline void DivListModel<ULL>::insert(const Entry& entry) {
  if (_minimizeOnInsert)
    _finder.insertReminimize(entry);
  else
    _finder.insert(entry);
} 

template<bool ULL>
inline std::string DivListModel<ULL>::getName() const {
  return _finder.getName() +
    (_minimizeOnInsert ? " remin" : "");
}

#endif
