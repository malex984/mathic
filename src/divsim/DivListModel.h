#ifndef DIV_ARRAY_MODEL_GUARD
#define DIV_ARRAY_MODEL_GUARD

#include "DivList.h"
#include "Monomial.h"
#include <string>
#include <vector>

/** Helper class for DivListModel. */
template<bool UseLinkedList, bool UseDivMask>
class DivListModelConfiguration;

template<bool ULL, bool UDM>
class DivListModelConfiguration {
public:
  typedef int Exponent;
  typedef ::Monomial Monomial;
  typedef Monomial Entry;

  DivListModelConfiguration
    (size_t varCount,
    bool sortOnInsert,
     double rebuildRatio,
     size_t minRebuild):
  _varCount(varCount),
  _sortOnInsert(sortOnInsert),
  _useAutomaticRebuild((rebuildRatio > 0.0 || minRebuild > 0) && UDM),
  _rebuildRatio(rebuildRatio),
  _minRebuild(minRebuild),
  _expQueryCount(0) {}

  static const bool UseLinkedList = ULL;
  static const bool UseDivMask = UDM;

  bool getDoAutomaticRebuilds() const {return _useAutomaticRebuild;}
  double getRebuildRatio() const {return _rebuildRatio;}
  size_t getRebuildMin() const {return _minRebuild;}
  bool getSortOnInsert() const {return _sortOnInsert;}


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

  bool isLessThan(const Monomial& a, const Monomial& b) const {
    for (size_t var = 0; var < getVarCount(); ++var) {
      if (getExponent(a, var) < getExponent(b, var))
        return true;
      if (getExponent(b, var) < getExponent(a, var))
        return false;
    }
    return false;
  }

  unsigned long long getExpQueryCount() const {return _expQueryCount;}

 private:
  const size_t _varCount;
  const bool _sortOnInsert;
  const bool _useAutomaticRebuild;
  const double _rebuildRatio;
  const size_t _minRebuild;
  mutable unsigned long long _expQueryCount;
};

template<bool UseLinkedList, bool UseDivMask>
class DivListModel;

/** An instantiation of the capabilities of DivList. */
template<bool ULL, bool UDM>
class DivListModel {
 private:
  typedef DivListModelConfiguration<ULL, UDM> C;
  typedef mathic::DivList<C> Finder;
 public:
  typedef typename Finder::iterator iterator;
  typedef typename Finder::const_iterator const_iterator;
  typedef typename Finder::Monomial Monomial;
  typedef typename Finder::Entry Entry;

  DivListModel(size_t varCount,
              bool minimizeOnInsert,
              bool moveDivisorToFront,
              bool sortOnInsert,
              double rebuildRatio,
              size_t minRebuild):
  _finder(C(varCount, sortOnInsert, rebuildRatio, minRebuild)),
  _minimizeOnInsert(minimizeOnInsert),
  _moveDivisorToFront(moveDivisorToFront) {
    ASSERT(!sortOnInsert || !moveDivisorToFront);
  }

  void insert(const Entry& entry);
  template<class MultipleOutput>
  void insert(const Entry& entry, MultipleOutput& removed);

  Entry* findDivisor(const Monomial& monomial) {
    return _finder.findDivisor(monomial);
  }
  const Entry* findDivisor(const Monomial& monomial) const {
    return const_cast<DivListModel<ULL, UDM>&>(*this).findDivisor(monomial);
  }

  template<class DO>
  void findAllDivisors(const Monomial& monomial, DO& out) {
    _finder.findAllDivisors(monomial, out);
  }
  template<class DO>
  void findAllDivisors(const Monomial& monomial, DO& out) const {
    _finder.findAllDivisors(monomial, out);
  }
  template<class EO>
  void forAll(EO& out) {
    _finder.forAll(out);
  }
  template<class EO>
  void forAll(EO& out) const {
    _finder.forAll(out);
  }

  std::string getName() const;

  size_t size() const {return _finder.size();}

  unsigned long long getExpQueryCount() const {
    return _finder.getConfiguration().getExpQueryCount();
  }

 private:
  Finder _finder;
  const bool _minimizeOnInsert;
  const bool _moveDivisorToFront;
};

template<bool ULL, bool UDM>
inline void DivListModel<ULL, UDM>::insert(const Entry& entry) {
  if (!_minimizeOnInsert) {
    _finder.insert(entry);
    return;
  }
  if (findDivisor(entry) != 0)
    return;
  _finder.removeMultiples(entry);
  _finder.insert(entry);
}

template<bool ULL, bool UDM>
template<class MO>
inline void DivListModel<ULL, UDM>::insert(const Entry& entry, MO& out) {
  if (!_minimizeOnInsert) {
    _finder.insert(entry);
    return;
  }
  if (findDivisor(entry) != 0)
    return;
  _finder.removeMultiples(entry, out);
  _finder.insert(entry);
}

template<bool ULL, bool UDM>
inline std::string DivListModel<ULL, UDM>::getName() const {
  return _finder.getName() +
    (_minimizeOnInsert ? " remin" : " nomin") +
    (_moveDivisorToFront ? " toFront" : "");
}

#endif
