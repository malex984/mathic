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
#include <vector>

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
  typedef KDTreeWalker<C> Walker;

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
  typedef std::reverse_iterator<iterator> reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

  KDTree(const C& configuration);

  bool removeMultiples(const Monomial& monomial);
  void insert(const Entry& entry);
  iterator findDivisor(const Monomial& monomial);
  const_iterator findDivisor(const Monomial& monomial) const {
    return const_cast<KDTree<C>&>(*this).findDivisor(monomial);
  }

  iterator begin() {return iterator::makeBegin(_root);}
  const_iterator begin() const {return const_cast<KDTree<C>&>(*this).begin();}
  iterator end() {return iterator::makeEnd(_root);}
  const_iterator end() const {return const_cast<KDTree<C>&>(*this).end();}

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

  template<class TreeIt>
  class IterHelper;

#ifdef DEBUG
  bool debugIsValid() const;
#endif

  std::vector<Node*> _tmp;
  Arena _arena;
  C _conf;
  Node* _root;
};

template<class C>
template<class LeafIt>
class KDTree<C>::IterHelper {
  typedef IterHelper<C> Helper;
  typedef KDTreeWalker<C> Walker;
  typedef KDTreeNode<C> Node;
 public:
  IterHelper() {} // singular iterator
  IterHelper(const IterHelper<C>& it):
   _walker(it.getWalker()), _leafIt(it.getLeafIterator()) {}
  template<class T>
  IterHelper(const T& it):
   _walker(it.getWalker()), _leafIt(it.getLeafIterator()) {}
  IterHelper(Leaf& leaf, LeafIt it): _walker(Walker::makeAt(&leaf, leaf.getParent())), _leafIt(it) {}
  IterHelper(const Walker& walker, LeafIt it): _walker(walker), _leafIt(it) {}
  static Helper makeBegin(Node* root) {
    Walker walker = Walker::makeAtFirstNonEmptyLeaf(root);
    LeafIt it = walker.atEnd() ? 0 : walker.asLeaf().begin();
    return Helper(walker, it);
  }
  static Helper makeEnd(Node* root) {return Helper(Walker::makeAtEnd(root), 0);}

  void assign(const IterHelper& it) {
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
class KDTree<C>::iterator :
public std::iterator<std::bidirectional_iterator_tag, Entry> {
  friend class KDTree;
  typedef typename Leaf::iterator LeafIt;
  typedef IterHelper<LeafIt> Helper;
 public:
  iterator() {} // singular iterator
  iterator(const iterator& it): _it(it._it) {}
  iterator& operator=(const iterator& it) {_it.assign(it._it); return *this;}

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
  iterator(Leaf& leaf, LeafIt it): _it(leaf, it) {}
  static iterator makeBegin(Node* root) {
    return iterator(Helper::makeBegin(root));
  }
  static iterator makeEnd(Node* root) {
    return iterator(Helper::makeEnd(root));
  }

 private:
  iterator(const Helper& helper): _it(helper) {}

  Helper _it;
};

template<class C>
class KDTree<C>::const_iterator :
public std::iterator<std::bidirectional_iterator_tag, Entry> {
  friend class KDTree;
  typedef typename Leaf::const_iterator LeafIt;
  typedef IterHelper<LeafIt> Helper;
  typedef typename KDTree<C>::iterator iterator;
 public:
  const_iterator() {} // singular iterator
  const_iterator(const typename KDTree<C>::iterator& it): _it(it._it) {}
  const_iterator(const const_iterator& it): _it(it._it) {}
  const_iterator& operator=(const const_iterator& it) {_it.assign(it._it); return *this;}

  const_iterator& operator++() {_it.increment(); return *this;}
  const_iterator operator++(int) {const_iterator tmp = *this; operator++(); return tmp;}
  const_iterator& operator--() {_it.decrement(); return *this;}
  const_iterator operator--(int) {const_iterator tmp = *this; operator--(); return tmp;}

  bool operator==(const iterator& it) const {return _it.equals(it._it);}
  bool operator!=(const iterator& it) const {return !(*this == it);}
  bool operator==(const const_iterator& it) const {return _it.equals(it._it);}
  bool operator!=(const const_iterator& it) const {return !(*this == it);}
  Entry& operator*() const {return *_it.get();}
  Entry* operator->() const {return &*_it.get();}

 protected:
  const_iterator(Leaf& leaf, LeafIt it): _it(leaf, it) {}
  static const_iterator makeBegin(Node* root) {
    return const_iterator(Helper::makeBegin(root));
  }
  static const_iterator makeEnd(Node* root) {
    return const_iterator(Helper::makeEnd(root));
  }

 private:
  const_iterator(const Helper& helper): _it(helper) {}

  Helper _it;
};

template<class C>
KDTree<C>::KDTree(const C& configuration): _conf(configuration) {
  ASSERT(_conf.getLeafSize() >= 2);
  _root = new (_arena.allocObjectNoCon<Leaf>()) Leaf(_arena, _conf);
}

template<class C>
size_t KDTree<C>::size() const {
  size_t sum = 0;
  for (Walker walker(_root); !walker.atEnd(); walker.next())
    if (walker.atLeaf())
      sum += walker.asLeaf().size();
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
NO_PINLINE bool KDTree<C>::removeMultiples(const Monomial& monomial) {
  ASSERT(_tmp.empty());
  bool changed = false;
  Node* node = _root;
  while (true) {
    while (node->isInterior()) {
      Interior& interior = node->asInterior();
      if (!(interior.getExponent() <
        _conf.getExponent(monomial, interior.getVar())))
        _tmp.push_back(&interior.getEqualOrLess());
      node = &interior.getStrictlyGreater();
    }
    ASSERT(node->isLeaf());
    if (node->asLeaf().removeMultiples(monomial, _conf))
      changed = true;
    if (_tmp.empty())
      break;
    node = _tmp.back();
    _tmp.pop_back();
  }
  ASSERT(debugIsValid());
  ASSERT(_tmp.empty());
  return changed;
}

template<class C>
NO_PINLINE void KDTree<C>::insert(const Entry& entry) {
  Node* node = _root;
  while (node->isInterior())
    node = &node->asInterior().getChildFor(entry, _conf);
  Leaf* leaf = &node->asLeaf();

  ASSERT(leaf->size() <= _conf.getLeafSize());
  if (leaf->size() == _conf.getLeafSize()) {
    Interior& interior = leaf->split(_arena, _conf);
    if (leaf == _root)
      _root = &interior;
    leaf = &interior.getChildFor(entry, _conf).asLeaf();
  }
  ASSERT(leaf->size() < _conf.getLeafSize());
  leaf->insert(entry, _conf);

  ASSERT(debugIsValid());
}

template<class C>
NO_PINLINE typename KDTree<C>::iterator KDTree<C>::findDivisor(const Monomial& monomial) {  
  ASSERT(_tmp.empty());
  Node* node = _root;
  while (true) {
    while (node->isInterior()) {
      Interior& interior = node->asInterior();
      if (interior.getExponent() <
        _conf.getExponent(monomial, interior.getVar()))
        _tmp.push_back(&interior.getStrictlyGreater());
      node = &interior.getEqualOrLess();
    }
    ASSERT(node->isLeaf());
    Leaf& leaf = node->asLeaf();
    LeafIt leafIt = leaf.findDivisor(monomial, _conf);
    if (leafIt != leaf.end()) {
      _tmp.clear();
      return iterator(leaf, leafIt);
    }
    if (_tmp.empty())
      break;
    node = _tmp.back();
    _tmp.pop_back();
  }
  ASSERT(_tmp.empty());
  return end();
}

#ifdef DEBUG
template<class C>
bool KDTree<C>::debugIsValid() const {
  ASSERT(_tmp.empty());
  Walker walker(_root);
  if (walker.atEnd())
    return true;
  ASSERT(walker.getNode()->getParent() == 0);
  for (; !walker.atEnd(); walker.next()) {
    if (walker.atLeaf()) {
      Leaf& leaf = walker.asLeaf();
      Walker ancestor(walker);
      // Check all interior nodes above leaf have each monomial in the
      // leaf in the correct child.
      while (!ancestor.atRoot()) {
        Node* child = ancestor.getNode();
        ancestor.toParent();
        ASSERT(!ancestor.atEnd());
        Interior& interior = ancestor.asInterior();
        if (child == &interior.getEqualOrLess()) {
          typename Leaf::const_iterator it = leaf.begin();
          for (; it != leaf.end(); ++it) {
            ASSERT(!(interior.getExponent() <
              _conf.getExponent(*it, interior.getVar())));
          }
        } else {
          ASSERT(child == &ancestor.getStrictlyGreater());
          typename Leaf::const_iterator it = leaf.begin();
          for (; it != leaf.end(); ++it) {
            ASSERT(interior.getExponent() <
              _conf.getExponent(*it, interior.getVar()));
          }
        }
      }
    } else {
      Interior& interior = walker.asInterior();
      ASSERT(walker.getEqualOrLess().getParent() == &interior);
      ASSERT(walker.getStrictlyGreater().getParent() == &interior);
      ASSERT(interior.getVar() < _conf.getVarCount());
    }
  }
  return true;
}
#endif

#endif
