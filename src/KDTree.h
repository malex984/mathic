#ifndef K_D_TREE_GUARD
#define K_D_TREE_GUARD

#include "KDTreeWalker.h"
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
  typedef KDTreeNode<C> Node;
  typedef KDTreeInterior<C> Interior;
  typedef KDTreeLeaf<C> Leaf;

  typedef typename Leaf::iterator LeafIt;
  typedef std::list<Leaf> Tree;
  typedef typename Tree::iterator TreeIt;
  typedef typename Tree::const_iterator CTreeIt;

 public:
  typedef typename C::Monomial Monomial;
  typedef typename C::Entry Entry;
  typedef typename C::Exponent Exponent;

  //class iterator;
  class iterator2;
  class const_iterator2;
  typedef iterator2 iterator;
  typedef const_iterator2 const_iterator;

  typedef std::reverse_iterator<iterator> reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

  KDTree(const C& configuration);

  bool removeMultiples(const Monomial& monomial);
  void insert(const Entry& entry);
  iterator findDivisor(const Monomial& monomial);
  const_iterator findDivisor(const Monomial& monomial) const {
    return const_cast<KDTree<C>&>(*this).findDivisor(monomial);
  }
  /*
  iterator begin() {return iterator(_leaves.begin(), _leaves.end());}
  const_iterator begin() const {return const_cast<KDTree<C>&>(*this).begin();}
  iterator end() {return iterator(_leaves.end());}
  const_iterator end() const {return const_cast<KDTree<C>&>(*this).end();}
  /*/
  iterator begin() {return iterator::makeBegin(_root);}
  const_iterator begin() const {return const_cast<KDTree<C>&>(*this).begin();}
  iterator end() {return iterator::makeEnd(_root);}
  const_iterator end() const {return const_cast<KDTree<C>&>(*this).end();}
  // */

  reverse_iterator rbegin() {return reverse_iterator(end());}
  const_reverse_iterator rbegin() const {return const_reverse_iterator(end());}
  reverse_iterator rend() {return reverse_iterator(begin());}
  const_reverse_iterator rend() const {return const_reverse_iterator(begin());}

  size_t size() const;

  std::string getName() const;

  C& getConfiguration() {return _conf;}
  const C& getConfiguration() const {return _conf;}

 private:
   KDTree(const KDTree<C>&); // unavailable
   void operator=(const KDTree<C>&); // unavailable

  template<class MacroIt, class MicroIt>
  class IterHelper;
  template<class TreeIt>
  class IterHelper2;

  Arena _arena;
  C _conf;
  Tree _leaves;
  KDTreeNode<C>* _root;
};

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
// ********************************

template<class C>
template<class LeafIt>
class KDTree<C>::IterHelper2 {
  typedef IterHelper2<C> Helper2;
  typedef KDTreeWalker<C> Walker;
  typedef KDTreeNode<C> Node;
 public:
  IterHelper2() {} // singular iterator
  IterHelper2(const IterHelper2<C>& it):
   _walker(it.getWalker()), _leafIt(it.getLeafIterator()) {}
  template<class T>
  IterHelper2(const T& it):
   _walker(it.getWalker()), _leafIt(it.getLeafIterator()) {}
  IterHelper2(Leaf& leaf, LeafIt it): _walker(Walker::makeAt(&leaf, leaf.getParent())), _leafIt(it) {}
  IterHelper2(const Walker& walker, LeafIt it): _walker(walker), _leafIt(it) {}
  static Helper2 makeBegin(Node* root) {
    Walker walker = Walker::makeAtFirstNonEmptyLeaf(root);
    LeafIt it = walker.atEnd() ? 0 : walker.asLeaf().begin();
    return Helper2(walker, it);
  }
  static Helper2 makeEnd(Node* root) {return Helper2(Walker::makeAtEnd(root), 0);}

  void assign(const IterHelper2& it) {
    _walker = it._walker;
    _leafIt = it._leafIt;
  }

  void increment() {
    ASSERT(!_walker.atEnd());
    ASSERT(_walker.atLeaf());
    ASSERT(!_walker.asLeaf().empty());
    ASSERT(_leafIt != _walker.asLeaf().end());
	++_leafIt;
    if (_leafIt == _walker.asLeaf().end()) {
      _walker.nextNonEmptyLeaf();
      ASSERT(_walker.atEnd() || !_walker.asLeaf().empty());
      _leafIt = _walker.atEnd() ? 0 : _walker.asLeaf().begin();
	}
    ASSERT((_walker.atEnd() && _leafIt == 0) ||
      _leafIt != _walker.asLeaf().end());
  }

  void decrement() {
    ASSERT(_walker.atEnd() || _walker.atLeaf());
    ASSERT((_walker.atEnd() && _leafIt == 0) ||
      _leafIt != _walker.asLeaf().end());
    if (!_walker.atEnd() && _leafIt != _walker.asLeaf().begin())
      --_leafIt;
    else {
      _walker.toPrevNonEmptyLeaf();
      ASSERT(!_walker.asLeaf().empty());
      _leafIt = _walker.asLeaf().end();
      --_leafIt;
    }
  }

  template<class T>
  bool equals(const T& it) const {
    return _walker == it.getWalker() && _leafIt == it.getLeafIterator();
  }

  const LeafIt& get() const {return _leafIt;}
  const LeafIt& getLeafIterator() const {return _leafIt;}
  const Walker& getWalker() const {return _walker;}

 private:
  Walker _walker;
  LeafIt _leafIt;
};

template<class C>
class KDTree<C>::iterator2 :
public std::iterator<std::bidirectional_iterator_tag, Entry> {
  friend class KDTree;
  typedef typename Leaf::iterator LeafIt;
  typedef IterHelper2<LeafIt> Helper;
 public:
  iterator2() {} // singular iterator
  iterator2(const iterator& it): _it(it._it) {}
  iterator2& operator=(const iterator& it) {_it.assign(it._it);}

  iterator2& operator++() {_it.increment(); return *this;}
  iterator2 operator++(int) {iterator tmp = *this; operator++(); return tmp;}
  iterator2& operator--() {_it.decrement(); return *this;}
  iterator2 operator--(int) {iterator tmp = *this; operator--(); return tmp;}

  bool operator==(const iterator2& it) const {return _it.equals(it._it);}
  bool operator!=(const iterator2& it) const {return !(*this == it);}
  bool operator==(const const_iterator& it) const {return _it.equals(it._it);}
  bool operator!=(const const_iterator& it) const {return !(*this == it);}
  Entry& operator*() const {return *_it.get();}
  Entry* operator->() const {return &*_it.get();}

 protected:
  iterator2(Leaf& leaf, LeafIt it): _it(leaf, it) {}
  static iterator2 makeBegin(Node* root) {
    return iterator2(Helper::makeBegin(root));
  }
  static iterator2 makeEnd(Node* root) {
    return iterator2(Helper::makeEnd(root));
  }

 private:
  iterator2(const Helper& helper): _it(helper) {}

  Helper _it;
};

template<class C>
class KDTree<C>::const_iterator2 :
public std::iterator<std::bidirectional_iterator_tag, Entry> {
  friend class KDTree;
  typedef typename Leaf::const_iterator LeafIt;
  typedef IterHelper2<LeafIt> Helper;
 public:
  const_iterator2() {} // singular iterator
  const_iterator2(const typename KDTree<C>::iterator2& it): _it(it._it) {}
  const_iterator2(const const_iterator2& it): _it(it._it) {}
  const_iterator2& operator=(const const_iterator2& it) {_it.assign(it._it);}

  const_iterator2& operator++() {_it.increment(); return *this;}
  const_iterator2 operator++(int) {const_iterator tmp = *this; operator++(); return tmp;}
  const_iterator2& operator--() {_it.decrement(); return *this;}
  const_iterator2 operator--(int) {const_iterator tmp = *this; operator--(); return tmp;}

  bool operator==(const iterator2& it) const {return _it.equals(it._it);}
  bool operator!=(const iterator2& it) const {return !(*this == it);}
  bool operator==(const const_iterator& it) const {return _it.equals(it._it);}
  bool operator!=(const const_iterator& it) const {return !(*this == it);}
  Entry& operator*() const {return *_it.get();}
  Entry* operator->() const {return &*_it.get();}

 protected:
  const_iterator2(Leaf& leaf, LeafIt it): _it(leaf, it) {}
  static const_iterator2 makeBegin(Node* root) {
    return const_iterator2(Helper::makeBegin(root));
  }
  static const_iterator2 makeEnd(Node* root) {
    return const_iterator2(Helper::makeEnd(root));
  }

 private:
  const_iterator2(const Helper& helper): _it(helper) {}

  Helper _it;
};

/*
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
*/
template<class C>
KDTree<C>::KDTree(const C& configuration): _conf(configuration) {
  _leaves.resize(1);
  _leaves.back().reset(_arena, _conf);
  _root = &_leaves.back();
}

template<class C>
size_t KDTree<C>::size() const {
  size_t sum = 0;
  for (CTreeIt it = _leaves.begin(); it != _leaves.end(); ++it)
    sum += it->size();
  ASSERT(sum == static_cast<size_t>(std::distance(begin(), end())));
  ASSERT(sum == static_cast<size_t>(std::distance(rbegin(), rend())));
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
bool KDTree<C>::removeMultiples(const Monomial& monomial) {
  bool changed = false;
  for (TreeIt it = _leaves.begin(); it != _leaves.end(); ++it)
    if (it->removeMultiples(monomial, _conf))
      changed = true;
  return changed;
}

template<class C>
void KDTree<C>::insert(const Entry& entry) {
  if (_leaves.back().size() >= _conf.getLeafSize()) {
    Leaf& toSplit = _leaves.back();
    _leaves.resize(_leaves.size() + 1);
    _leaves.back().reset(_arena, _conf);
    KDTreeInterior<C>* node = _arena.allocObject<KDTreeInterior<C> >();
    node->reset(toSplit, _leaves.back(), _arena, _conf);
    if (&toSplit == _root)
      _root = node;
  }
  _leaves.back().insert(entry, _conf);
}

template<class C>
typename KDTree<C>::iterator KDTree<C>::findDivisor(const Monomial& monomial) {
  for (TreeIt it = _leaves.begin(); it != _leaves.end(); ++it) {
    LeafIt leafIt = it->findDivisor(monomial, _conf);
    if (leafIt != it->end())
      return iterator(*it, leafIt);
  }
  return end();
}

#endif
