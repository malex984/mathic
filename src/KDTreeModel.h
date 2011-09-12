#ifndef K_D_TREE_MODEL_GUARD
#define K_D_TREE_MODEL_GUARD

#include "KDTree.h"
#include <string>
#include <vector>

/** Helper class for KDTreeModel. */
class KDTreeModelConfiguration {
 public:
  typedef int Exponent;
  typedef std::vector<Exponent> Monomial;
  typedef Monomial Entry;

  KDTreeModelConfiguration
    (size_t varCount,
     size_t leafSize,
      bool sortOnInsert,
      bool useDivisorCache):
   _varCount(varCount),
   _leafSize(leafSize),
   _sortOnInsert(sortOnInsert),
   _useDivisorCache(useDivisorCache),
   _expQueryCount(0) {}

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

  size_t getLeafSize() const {return _leafSize;}
  bool getUseDivisorCache() const {return _useDivisorCache;}

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
  mutable unsigned long long _expQueryCount;
};

/** An instantiation of the capabilities of KDTree. */
class KDTreeModel {
 private:
  typedef KDTreeModelConfiguration C;
  typedef KDTree<C> Finder;
 public:
  typedef Finder::iterator iterator;
  typedef Finder::const_iterator const_iterator;
  typedef Finder::Monomial Monomial;
  typedef Finder::Entry Entry;

 KDTreeModel(size_t varCount,
			 size_t leafSize,
			 bool minimizeOnInsert,
			 bool sortOnInsert,
             bool useDivisorCache):
  _finder(C(varCount, leafSize, sortOnInsert, useDivisorCache)),
  _minimizeOnInsert(minimizeOnInsert) {}

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
  size_t size() const {return _finder.size();}

  unsigned long long getExpQueryCount() const {
    return _finder.getConfiguration().getExpQueryCount();
  }

  class Comparer;

 private:
  Finder _finder;
  bool _minimizeOnInsert;
};

inline void KDTreeModel::insert(const Entry& entry) {
  if (!_minimizeOnInsert) {
    _finder.insert(entry);
    return;
  }
  if (findDivisor(entry) != _finder.end())
    return;
  _finder.removeMultiples(entry);
  _finder.insert(entry);
}

inline std::string KDTreeModel::getName() const {
  return _finder.getName() +
    (_minimizeOnInsert ? " remin" : " nomin");
}

#endif
