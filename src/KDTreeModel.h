#ifndef K_D_TREE_MODEL_GUARD
#define K_D_TREE_MODEL_GUARD

#include "Monomial.h"
#include "KDTree.h"
#include <string>
#include <vector>

template<bool UseDivMask>
class KDTreeModelConfiguration;

/** Helper class for KDTreeModel. */
template<bool UDM>
class KDTreeModelConfiguration {
 public:
  typedef int Exponent;
  typedef ::Monomial Monomial;
  typedef Monomial Entry;

  KDTreeModelConfiguration
    (size_t varCount,
     size_t leafSize,
     bool sortOnInsert,
     bool useDivisorCache,
     double rebuildRatio,
     size_t minRebuild):
   _varCount(varCount),
   _leafSize(leafSize),
   _sortOnInsert(sortOnInsert),
   _useDivisorCache(useDivisorCache),
   _useAutomaticRebuild(rebuildRatio > 0.0),
   _rebuildRatio(rebuildRatio),
   _minRebuild(minRebuild),
   _expQueryCount(0) {
     ASSERT(rebuildRatio >= 0);
   }

  size_t getVarCount() const {return _varCount;}
  bool getSortOnInsert() const {return _sortOnInsert;}

  Exponent getExponent(const Monomial& monomial, size_t var) const {
    ++_expQueryCount;
    ASSERT(var < monomial.size());
    return monomial[var];
  }

  NO_PINLINE bool divides(const Monomial& a, const Monomial& b) const {
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

  size_t getLeafSize() const {return _leafSize;}
  bool getUseDivisorCache() const {return _useDivisorCache;}
  bool getDoAutomaticRebuilds() const {return _useAutomaticRebuild;}
  double getRebuildRatio() const {return _rebuildRatio;}
  size_t getRebuildMin() const {return _minRebuild;}

  static const bool UseDivMask = UDM;

  unsigned long long getExpQueryCount() const {return _expQueryCount;}

  class Comparer {
  public:
  Comparer(const KDTreeModelConfiguration& conf): _conf(conf) {}
	bool operator()(const Entry& a, const Entry& b) const {
	  return _conf.isLessThan(a, b);
	}
  private:
	const KDTreeModelConfiguration& _conf;
  };
  Comparer getComparer() const {return *this;}

 private:
  const size_t _varCount;
  const size_t _leafSize;
  const bool _sortOnInsert;
  const bool _useDivisorCache;
  const bool _useAutomaticRebuild;
  const double _rebuildRatio;
  const size_t _minRebuild;
  mutable unsigned long long _expQueryCount;
};

/** An instantiation of the capabilities of KDTree. */
template<bool UseDivMask = true>
class KDTreeModel {
 private:
  typedef KDTreeModelConfiguration<UseDivMask> C;
  typedef KDTree<C> Finder;
 public:
  typedef typename Finder::iterator iterator;
  typedef typename Finder::const_iterator const_iterator;
  typedef typename Finder::Monomial Monomial;
  typedef typename Finder::Entry Entry;

 KDTreeModel(size_t varCount,
			 size_t leafSize,
			 bool minimizeOnInsert,
			 bool sortOnInsert,
             bool useDivisorCache,
             double rebuildRatio,
             size_t minRebuild):
  _finder(C(varCount, leafSize, sortOnInsert, useDivisorCache, rebuildRatio, minRebuild)),
  _minimizeOnInsert(minimizeOnInsert) {}

  void insert(const Entry& entry);
  iterator findDivisor(const Monomial& monomial) {
    return _finder.findDivisor(monomial);
  }
  const_iterator findDivisor(const Monomial& monomial) const {
    return _finder.findDivisor(monomial);
  }
  std::string getName() const;

  template<class DO>
  void findAllDivisors(const Monomial& monomial, DO& out) {
    _finder.findAllDivisors(monomial, out);
  }
  template<class DO>
  void findAllDivisors(const Monomial& monomial, DO& out) const {
    _finder.findAllDivisors(monomial, out);
  }

  iterator begin() {return _finder.begin();}
  const_iterator begin() const {return _finder.begin();}
  iterator end() {return _finder.end();}
  const_iterator end() const {return _finder.end();}
  size_t size() const {return _finder.size();}

  unsigned long long getExpQueryCount() const {
    return _finder.getConfiguration().getExpQueryCount();
  }

  class Comparer;

 private:
  Finder _finder;
  bool _minimizeOnInsert;
};

template<bool UDM>
inline void KDTreeModel<UDM>::insert(const Entry& entry) {
  if (!_minimizeOnInsert) {
    _finder.insert(entry);
    return;
  }
  if (findDivisor(entry) != _finder.end())
    return;
  _finder.removeMultiples(entry);
  _finder.insert(entry);
}

template<bool UDM>
inline std::string KDTreeModel<UDM>::getName() const {
  return _finder.getName() +
    (_minimizeOnInsert ? " remin" : " nomin");
}

#endif
