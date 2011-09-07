#ifndef K_D_TREE_GUARD
#define K_D_TREE_GUARD

#include "KDTreeLeaf.h"
#include "Arena.h"
#include <list>
#include <string>
#include <algorithm>
#include <iterator>
#include <sstream>

/** An object that supports queries for divisors of a monomial using
 a KD Tree (K Dimensional Tree). See DivFinder for more documentation.

 Extra fields for Configuration:

 * bool getSortOnInsert() const
  Keep the monomials sorted to speed up queries.
*/
template<class Configuration>
class KDTree;

template<class C>
class KDTree {
 private:
  typedef KDTreeLeaf<C> Leaf;
  typedef typename Leaf::iterator LeafIt;
  typedef std::list<Leaf> Tree;
  typedef typename Tree::iterator TreeIt;
  typedef typename Tree::const_iterator CTreeIt;

 public:
  typedef typename C::Monomial Monomial;
  typedef typename C::Entry Entry;
  typedef typename C::Exponent Exponent;

  class iterator;
  class const_iterator;

  KDTree(const C& configuration);

  bool removeMultiples(const Monomial& monomial) {
    bool changed = false;
    for (TreeIt it = _leaves.begin(); it != _leaves.end(); ++it)
      if (it->removeMultiples(monomial, _conf))
        changed = true;
    return changed;
  }

  void insert(const Entry& entry) {
    if (_leaves.empty() || _leaves.back().size() >= _conf.getLeafSize()) {
      _leaves.resize(_leaves.size() + 1);
      _leaves.back().reset(_arena, _conf);
    }
    _leaves.back().insert(entry, _conf);
  }

  iterator findDivisor(const Monomial& monomial) {
    for (TreeIt it = _leaves.begin(); it != _leaves.end(); ++it) {
      LeafIt leafIt = it->findDivisor(monomial, _conf);
      if (leafIt != it->end())
        return iterator(it, _leaves.end(), leafIt);
    }
    return end();
  }
  const_iterator findDivisor(const Monomial& monomial) const {
    return const_cast<KDTree<C>&>(*this).findDivisor(monomial);
  }
  
  iterator begin() {return iterator(_leaves.begin(), _leaves.end());}
  const_iterator begin() const {return const_cast<KDTree<C>&>(*this).begin();}
  iterator end() {return iterator(_leaves.end());}
  const_iterator end() const {return const_cast<KDTree<C>&>(*this).end();}
  size_t size() const;

  std::string getName() const;

  C& getConfiguration() {return _conf;}
  const C& getConfiguration() const {return _conf;}

 private:
  template<class MacroIt, class MicroIt>
  class IterHelper;

  Arena _arena;
  C _conf;
  Tree _leaves;
};

template<class C>
KDTree<C>::KDTree(const C& configuration): _conf(configuration) {
}

template<class C>
size_t KDTree<C>::size() const {
  size_t sum = 0;
  for (CTreeIt it = _leaves.begin(); it != _leaves.end(); ++it)
    sum += it->size();
  ASSERT(sum == static_cast<size_t>(std::distance(begin(), end())));
  return sum;
}

template<class C>
std::string KDTree<C>::getName() const {
  std::stringstream out;
  out << "KDTree-" << _conf.getLeafSize()
    << (_conf.getSortOnInsert() ? " sort" : "");
  return out.str();
}

template<class C>
template<class MacroIt, class MicroIt>
class KDTree<C>::IterHelper {
 public:
  IterHelper() {} // singular iterator
  IterHelper(const IterHelper& it):
   _macro(it._macro), _macroEnd(it._macroEnd), _micro(it._micro) {
    avoidEndOfMacro();
  }
  IterHelper(MacroIt macroEnd):
   _macro(macroEnd), _macroEnd(macroEnd) {
    avoidEndOfMacro();
  }
  IterHelper(const MacroIt& macro, const MacroIt& macroEnd):
   _macro(macro), _macroEnd(macroEnd) {
    if (macro != macroEnd)
      _micro = macro->begin();
    avoidEndOfMacro();
   }
  IterHelper
    (const MacroIt& macro, const MacroIt& macroEnd, const MicroIt& micro):
   _macro(macro), _macroEnd(macroEnd), _micro(micro) {}

  void assign(const IterHelper& it) {
	_macro = it._macro;
	_macroEnd = it._macroEnd;
    _micro = it._micro;
  }

  void increment() {
	ASSERT(_macro != _macroEnd);
	++_micro;
    if (_micro == _macro->end()) {
	  ++_macro;
	  if (_macro != _macroEnd)
		_micro = _macro->begin();
      avoidEndOfMacro();
	}
  }

  void decrement() {
	while (_macro == _macroEnd || _micro == _macro->begin()) {
	  --_macro;
	  _micro = _macro->end();
	}
	--_micro;
	return *this;
  }

  template<class T>
  bool equals(const T& it) const {
	if (_macro == _macroEnd)
	  return _macro == it._macro;
	else if (it._macro == it._macroEnd)
	  return false;
	else
	  return _micro == it._micro && _macro == it._macro;
  }

  const MicroIt& get() const {return _micro;}

 private:
  void avoidEndOfMacro() {
    if (_macro == _macroEnd)
      return;
    while (_micro == _macro->end()) {
      ++_macro;
      if (_macro == _macroEnd)
        return;
      _micro = _macro->begin();
    }
  }

  MacroIt _macro;
  MacroIt _macroEnd;
  MicroIt _micro;
};

template<class C>
class KDTree<C>::iterator :
public std::iterator<std::bidirectional_iterator_tag, Entry> {
  friend class KDTree;
 private:
  typedef typename Tree::const_iterator MacroIt;
  typedef typename Leaf::const_iterator MicroIt;
 public:
  iterator() {} // singular iterator
  iterator(const iterator& it): _it(it._it) {}
  iterator& operator=(const iterator& it) {_it.assign(it._it);}

  iterator& operator++() {_it.increment(); return *this;}
  iterator operator++(int) {iterator tmp = *this; operator++(); return tmp;}
  iterator& operator--() {_it.decrement(); return *this;}
  iterator operator--(int) {iterator tmp = *this; operator--(); return tmp;}

  bool operator==(const iterator& it) const {return _it.equals(it._it);}
  bool operator!=(const iterator& it) const {return !(*this == it);}
  bool operator==(const const_iterator& it) const {return _it.equals(it._it);}
  bool operator!=(const const_iterator& it) const {return !(*this == it);}
  Entry& operator*() const {return *_it.get();}
  Entry* operator->() const {return &*_it.get();}

 protected:
  iterator
    (const MacroIt& macro, const MacroIt& macroEnd, const MicroIt& micro):
   _it(macro, macroEnd, micro) {}
  iterator(const MacroIt& macro, const MacroIt& macroEnd):
   _it(macro, macroEnd) {}
  iterator(const MacroIt& macroEnd): _it(macroEnd) {}

 private:
  IterHelper<MacroIt, MicroIt> _it;
};

template<class C>
class KDTree<C>::const_iterator :
public std::iterator<std::bidirectional_iterator_tag, const Entry> {
  friend class KDTree;
 private:
  typedef typename Tree::const_iterator MacroIt;
  typedef typename Leaf::const_iterator MicroIt;
 public:
  const_iterator() {} // singular iterator
  const_iterator(const typename KDTree<C>::iterator& it): _it(it._it) {}
  const_iterator(const const_iterator& it): _it(it._it) {}
  const_iterator& operator=(const const_iterator& it) {_it.assign(it._it);}

  const_iterator& operator++() {_it.increment(); return *this;}
  const_iterator operator++(int) {
	const_iterator tmp = *this;
	operator++();
	return tmp;
  }
  const_iterator& operator--() {_it.decrement(); return *this;}
  const_iterator operator--(int) {
	const_iterator tmp = *this;
	operator--();
	return tmp;
  }
  
  bool operator==(const iterator& it) const {return _it.equals(it._it);}
  bool operator!=(const iterator& it) const {return !(*this == it);}
  bool operator==(const const_iterator& it) const {return _it.equals(it._it);}
  bool operator!=(const const_iterator& it) const {return !(*this == it);}
  const Entry& operator*() const {return *_it.get();}
  const Entry* operator->() const {return &*_it.get();}

 protected:
  const_iterator
    (const MacroIt& macro, const MacroIt& macroEnd, const MicroIt& micro):
   _it(macro, macroEnd, micro) {}
  const_iterator(const MacroIt& macro, const MacroIt& macroEnd):
   _it(macro, macroEnd) {}
  const_iterator(const MacroIt& macroEnd): _it(macroEnd) {}

 private:
  IterHelper<MacroIt, MicroIt> _it;
};  

#endif
