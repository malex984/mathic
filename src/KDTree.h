#ifndef K_D_TREE_GUARD
#define K_D_TREE_GUARD

#include "KDTreeWalker.h"
#include "KDTreeLeaf.h"
#include "ArenaVector.h"
#include "Arena.h"
#include "DivMask.h"
#include <list>
#include <string>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <vector>

/** An object that supports queries for divisors of a monomial using
 a KD Tree (K Dimensional Tree). See DivFinder.h for more documentation.

 Extra fields for Configuration:

 * bool getSortOnInsert() const
  Return true to keep the monomials in leaves sorted to speed up queries.

 * size_t getLeafSize() const
  Return the fixed maximal size of a leaf.

 * static const bool UseDivMask
  Set to true to use div masks to speed up queries. This must be a
  static const data member - not a function.

 * size_t getUseAutomaticRebuild() const
 * double getRebuildRatio() const
 * size_t getRebuildMin() const
  If getUseAutomaticRebuild() returns true, the tree will call rebuild
  on itself after a total of max(size() * getRebuildRatio(), getRebuildMin())
  entry insertions and removals have occurred.
*/
template<class Configuration>
class KDTree;

template<class C>
class KDTree {
public:
  typedef C Configuration;
  static const bool UseDivMask = C::UseDivMask;
  typedef typename C::Monomial Monomial;
  typedef typename C::Entry Entry;
  typedef typename C::Exponent Exponent;

private:
  typedef typename DivMask::Extender<Entry, C::UseDivMask> ExtEntry;
  typedef typename DivMask::Extender<const Monomial&,C::UseDivMask> ExtMonoRef;
  typedef typename DivMask::Calculator<C> DivMaskCalculator;

  typedef KDTreeNode<C, ExtEntry> Node;
  typedef KDTreeInterior<C, ExtEntry> Interior;
  typedef KDTreeLeaf<C, ExtEntry> Leaf;
  typedef KDTreeWalker<C, ExtEntry> Walker;

  typedef typename Leaf::iterator LeafIt;

public:
  /** Iterator for enumerating all entries. */
  class iterator;

  /** Const iterator for enumerating all entries. */
  class const_iterator;

  /** Reverse iterator for enumerating all entries. */
  typedef std::reverse_iterator<iterator> reverse_iterator;

  /** Reverse const iterator for enumerating all entries. */
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

  /** Constructs an object with the given configuration. The configuration
   is copied into the object, so a reference to the passed-in object is
   not kept. */
  KDTree(const C& configuration);
  ~KDTree();

  /** Removes all multiples of monomial. A duplicate of monomial counts
   as a multiple. Returns true if any multiples were removed. */
  bool removeMultiples(const Monomial& monomial);

  /** Inserts entry into the data structure. Does NOT remove multiples
   of entry and entry is inserted even if it is a multiple of another
   entry. */
  void insert(const Entry& entry);

  /** Inserts the entries in the range [begin, end) into the data
   structure. Does NOT remove multiples of entry and entry is inserted
   even if it is a multiple of another entry.
   
   The elements in the range [begin, end) may be rearranged by this
   function, so the range must be mutable. If that is not acceptable,
   call the one element insert method for each element. */
  template<class Iter>
  void insert(Iter begin, Iter end);

  /** Returns the position of a divisor of monomial. Returns end() if no
   entries divide monomial. */
  iterator findDivisor(const Monomial& monomial);

  /** Calls output.push_back(entry) for each entry that divides monomial
   So a std::vector will work, but any other class with a push_back method
   accepting an iterator will also work. */
  template<class DivisorOutput>
  void findAllDivisors(const Monomial& monomial, DivisorOutput& out);

  /** Calls output.push_back(entry) for each entry that divides monomial
   So a std::vector will work, but any other class with a push_back method
   accepting an iterator will also work. */
  template<class DivisorOutput>
  void findAllDivisors(const Monomial& monomial, DivisorOutput& out) const;

  /** Returns the position of a divisor of monomial. Returns end() if no
   entries divide monomial. */
  const_iterator findDivisor(const Monomial& monomial) const {
    return const_cast<KDTree<C>&>(*this).findDivisor(monomial);
  }

  /** Removes all entries. Does not reset the configuration object. */
  void clear();

  /** Rebuilds the data structure. */
  void rebuild();

  iterator begin() {return iterator::makeBegin(_root);}
  const_iterator begin() const {return const_cast<KDTree<C>&>(*this).begin();}
  iterator end() {return iterator::makeEnd(_root);}
  const_iterator end() const {return const_cast<KDTree<C>&>(*this).end();}

  reverse_iterator rbegin() {return reverse_iterator(end());}
  const_reverse_iterator rbegin() const {return const_reverse_iterator(end());}
  reverse_iterator rend() {return reverse_iterator(begin());}
  const_reverse_iterator rend() const {return const_reverse_iterator(begin());}

  /** Returns whether there are entries. Not O(1) so don't call it in a loop.
   @todo Make this faster than size(). */
  bool empty() const {return size() == 0;}

  /** Returns the number of entries. Not O(1) so don't call it in a loop. */
  size_t size() const;

  /** Returns a string that describes the data structure. */
  std::string getName() const;

  /** Returns a reference to this object's configuration object. */
  C& getConfiguration() {return _conf;}

  /** Returns a reference to this object's configuration object. */
  const C& getConfiguration() const {return _conf;}

 private:
  KDTree(const KDTree<C>&); // unavailable
  void operator=(const KDTree<C>&); // unavailable

  void resetNumberOfChangesTillRebuild();
  void reportChanges(size_t changesMadeCount);

  /** Encapsulates common code between iterator and const_iterator. */
  template<class TreeIt>
  class IterHelper;

  /** Transfers push_back from entry to const entry. */
  template<class DO>
  class ConstDivisorOutput {
  public:
    ConstDivisorOutput(DO& out): _out(out) {}
    void push_back(Entry& entry) {
      const Entry& constEntry = entry;
      _out.push_back(constEntry);
    }
  private:
    DO& _out;
  };

  template<class Iter>
  struct InsertTodo {
    Iter begin;
    Iter end;
    Interior* parent;
  };

#ifdef DEBUG
  bool debugIsValid() const;
#endif

  std::vector<Node*> _tmp;
  Arena _arena;
  C _conf;
  Node* _root;
  iterator _divisorCache; /// The divisor in the previous query. Can be end().
  size_t _changesTillRebuild; /// Update using reportChanges().
  DivMaskCalculator _divMaskCalculator;
};

template<class C>
template<class LeafIt>
class KDTree<C>::IterHelper {
  typedef IterHelper<C> Helper;
  typedef KDTreeWalker<C, ExtEntry> Walker;
  typedef KDTreeNode<C, ExtEntry> Node;
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
  Entry& operator*() const {return _it.get()->get();}
  Entry* operator->() const {return &_it.get()->get();}

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
public std::iterator<std::bidirectional_iterator_tag, const Entry> {
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
  Entry& operator*() const {return _it.get()->get();}
  Entry* operator->() const {return &_it.get()->get();}

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
KDTree<C>::KDTree(const C& configuration):
_conf(configuration),
_divMaskCalculator(configuration) {
  ASSERT(_conf.getLeafSize() >= 2);
  _root = new (_arena.allocObjectNoCon<Leaf>()) Leaf(_arena, _conf);
  if (_conf.getUseDivisorCache())
    _divisorCache = end();
  resetNumberOfChangesTillRebuild();
  ASSERT(debugIsValid());
}

template<class C>
KDTree<C>::~KDTree() {
  clear();
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
  out << "KDTree leaf:" << _conf.getLeafSize();
  if (UseDivMask && _conf.getDoAutomaticRebuilds()) {
    out << " autob:" << _conf.getRebuildRatio()
      << '/' << _conf.getRebuildMin();
  }
  out << (C::UseDivMask && !C::UseTreeDivMask ? " dmask" : "")
    << (C::UseTreeDivMask ? " tree-dmask" : "")
    << (_conf.getSortOnInsert() ? " sort" : "")
    << (_conf.getUseDivisorCache() ? " cache" : "");
  return out.str();
}

/// @todo: return number instead of bool
template<class C>
NO_PINLINE bool KDTree<C>::removeMultiples(const Monomial& monomial) {
  ExtMonoRef extMonomial(monomial, _divMaskCalculator, _conf);

  ASSERT(_tmp.empty());
  size_t removedCount = 0;
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
    removedCount += node->asLeaf().removeMultiples(extMonomial, _conf);
    if (_tmp.empty())
      break;
    node = _tmp.back();
    _tmp.pop_back();
  }
  ASSERT(debugIsValid());
  ASSERT(_tmp.empty());
  if (_conf.getUseDivisorCache() && removedCount > 0)
    _divisorCache = end();
  reportChanges(removedCount);
  return removedCount > 0;
}

template<class C>
NO_PINLINE void KDTree<C>::insert(const Entry& entry) {
  ExtEntry extEntry(entry, _divMaskCalculator, _conf);

  Node* node = _root;
  while (node->isInterior()) {
    node->updateToLowerBound(extEntry);
    node = &node->asInterior().getChildFor(extEntry, _conf);
  }
  Leaf* leaf = &node->asLeaf();

  ASSERT(leaf->size() <= _conf.getLeafSize());
  if (leaf->size() == _conf.getLeafSize()) {
    Interior& interior = leaf->split(_arena, _conf);
    interior.updateToLowerBound(extEntry);
    if (leaf == _root)
      _root = &interior;
    leaf = &interior.getChildFor(extEntry, _conf).asLeaf();
  }
  ASSERT(leaf->size() < _conf.getLeafSize());
  leaf->insert(extEntry, _conf);

  if (_conf.getUseDivisorCache())
    _divisorCache = end();
  ASSERT(debugIsValid());
  reportChanges(1);
}

/// @todo: this function is too big and it knows too much about the details
/// inside nodes. Also, it allocates a std::vector every time.
template<class C>
template<class Iter>
NO_PINLINE void KDTree<C>::insert(Iter insertBegin, Iter insertEnd) {
  if (!empty()) {
    for (; insertBegin != insertEnd; ++insertBegin)
      insert(*insertBegin);
    return;
  } else if (insertBegin == insertEnd)
    return;

  _arena.freeAll();
  _root = 0;
  _divMaskCalculator.rebuild(insertBegin, insertEnd, _conf);
  if (_conf.getUseDivisorCache())
    _divisorCache = this->end();

  typedef InsertTodo<Iter> Task;
  typedef std::vector<Task> TaskCont;
  TaskCont todo;

  Interior* parent = 0;
  bool isEqualOrLessChild = false;
  while (true) {
    Node* node = 0;
    const size_t insertCount = std::distance(insertBegin, insertEnd);
    const bool isLeaf = (insertCount <= _conf.getLeafSize());
    if (isLeaf)
      node = Leaf::makeLeafCopy(parent, insertBegin,
      insertEnd, _arena, _divMaskCalculator, _conf);
    else {
      std::pair<Interior*, Iter> p =
        Node::preSplit(parent, insertBegin, insertEnd, _arena, _conf);
      ASSERT(p.second != insertBegin && p.second != insertEnd);
      // push strictly-greater on todo
      Task task;
      task.begin = p.second;
      task.end = insertEnd;
      task.parent = p.first;
      todo.push_back(task);
      // set up equal-or-less
      insertEnd = p.second;
      node = p.first;
    }

    if (parent == 0)
      _root = node;
    else if (isEqualOrLessChild)
      parent->setEqualOrLess(node);
    else
      parent->setStrictlyGreater(node);

    if (isLeaf) {
      // grab next item from todo
      if (todo.empty())
        break;
      Task task = todo.back();
      todo.pop_back();
      insertBegin = task.begin;
      insertEnd = task.end;
      parent = task.parent;
      // only strictly-greater goes on todo
      isEqualOrLessChild = false;
    } else {
      // continue with equal-or-less as next item
      parent = &node->asInterior();
      isEqualOrLessChild = true;
    }
  }
  if (_conf.getUseDivisorCache())
    _divisorCache = this->end();

  // set div masks
  if (C::UseTreeDivMask) {
    for (Walker walker(_root); !walker.atEnd(); walker.next())
      if (walker.atInterior()) {
        walker.getNode()->updateToLowerBound(walker.asInterior().getEqualOrLess());
        walker.getNode()->updateToLowerBound(walker.asInterior().getStrictlyGreater());
      }
  }

  ASSERT(debugIsValid());
  // range insert into empty container is equivalent to a rebuild.
  resetNumberOfChangesTillRebuild();
}

template<class C>
NO_PINLINE typename KDTree<C>::iterator KDTree<C>::findDivisor(const Monomial& monomial) {
  if (_conf.getUseDivisorCache() &&
    _divisorCache != end() && _conf.divides(*_divisorCache, monomial))
    return _divisorCache;

  ExtMonoRef extMonomial(monomial, _divMaskCalculator, _conf);

  ASSERT(debugIsValid());
  ASSERT(_tmp.empty());
  Node* node = _root;
  while (true) {
    while (node->isInterior()) {
      if (C::UseTreeDivMask &&
        !node->getDivMask().canDivide(extMonomial.getDivMask()))
        goto next;

      Interior& interior = node->asInterior();
      if (interior.getExponent() <
        _conf.getExponent(monomial, interior.getVar()))
        _tmp.push_back(&interior.getStrictlyGreater());
      node = &interior.getEqualOrLess();
    }

    {
      ASSERT(node->isLeaf());
      Leaf& leaf = node->asLeaf();
      LeafIt leafIt = leaf.findDivisor(extMonomial, _conf);
      if (leafIt != leaf.end()) {
        ASSERT(_conf.divides(leafIt->get(), extMonomial.get()));
        _tmp.clear();
        return _divisorCache = iterator(leaf, leafIt);
      }
    }
  next:
    if (_tmp.empty())
      break;
    node = _tmp.back();
    _tmp.pop_back();
  }
  ASSERT(_tmp.empty());
  return end();
}

template<class C>
template<class DO>
void KDTree<C>::findAllDivisors(const Monomial& monomial, DO& output) {
  ExtMonoRef extMonomial(monomial, _divMaskCalculator, _conf);

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
    leaf.findAllDivisors(extMonomial, output, _conf);
    if (_tmp.empty())
      break;
    node = _tmp.back();
    _tmp.pop_back();
  }
  ASSERT(_tmp.empty());
}

template<class C>
template<class DO>
void KDTree<C>::findAllDivisors(const Monomial& monomial, DO& output) const {
  ConstDivisorOutput<DO> constOutput(output);
  const_cast<KDTree<C>&>(*this).findAllDivisors(monomial, constOutput);
}

template<class C>
void KDTree<C>::clear() {
  for (Walker walker(_root); !walker.atEnd(); walker.next())
    if (walker.atLeaf())
      walker.asLeaf().clear();
  _arena.freeAll();
  _root = new (_arena.allocObjectNoCon<Leaf>()) Leaf(_arena, _conf);
  resetNumberOfChangesTillRebuild();
  _divMaskCalculator.rebuildDefault(_conf);
  if (_conf.getUseDivisorCache())
    _divisorCache = end();
}

template<class C>
void KDTree<C>::rebuild() {
  const size_t totalSize = size();
  typedef ArenaVector<Entry, true> TmpContainer;
  TmpContainer tmpCopy(Arena::getArena(), totalSize);
  for (Walker walker(_root); !walker.atEnd(); walker.next()) {
    if (walker.atLeaf()) {
      Leaf& leaf = walker.asLeaf();
      typename Leaf::const_iterator stop = leaf.end();
      for (typename Leaf::const_iterator it = leaf.begin(); it != stop; ++it)
        tmpCopy.push_back(it->get());
      leaf.clear();
    }
  }
  _arena.freeAll();
  _root = new (_arena.allocObjectNoCon<Leaf>()) Leaf(_arena, _conf);
  resetNumberOfChangesTillRebuild();
  // range insert into empty container IS a rebuild.
  insert(tmpCopy.begin(), tmpCopy.end());
}

template<class C>
void KDTree<C>::resetNumberOfChangesTillRebuild() {
  if (!_conf.getDoAutomaticRebuilds())
    return;
  ASSERT(_conf.getRebuildRatio() > 0);
  _changesTillRebuild = std::max
    (static_cast<size_t>(size() * _conf.getRebuildRatio()),
    _conf.getRebuildMin());
}

template<class C>
void KDTree<C>::reportChanges(size_t changesMadeCount) {
  // note how negative value/overflow of _changesTillRebuild cannot
  // happen this way.
  if (!_conf.getDoAutomaticRebuilds())
    return;
  if (_changesTillRebuild > changesMadeCount)
    _changesTillRebuild -= changesMadeCount;
  else
    rebuild();
}

#ifdef DEBUG
template<class C>
bool KDTree<C>::debugIsValid() const {
  ASSERT(_tmp.empty());
  ASSERT(!_conf.getDoAutomaticRebuilds() || _conf.getRebuildRatio() > 0);
  Walker walker(_root);
  if (walker.atEnd())
    return true;
  ASSERT(walker.getNode()->getParent() == 0);
  for (; !walker.atEnd(); walker.next()) {
    if (walker.atLeaf()) {
      Leaf& leaf = walker.asLeaf();
      typedef typename Leaf::const_iterator LeafCIter;
      if (C::UseTreeDivMask) {
        for (LeafCIter it = leaf.begin(); it != leaf.end(); ++it) {
          ASSERT(leaf.getDivMask().canDivide(it->getDivMask()));
        }
      }

      Walker ancestor(walker);
      // Check all interior nodes above leaf have each monomial in the
      // leaf in the correct child. Also check div masks.
      while (!ancestor.atRoot()) {
        Node* child = ancestor.getNode();
        ancestor.toParent();
        ASSERT(!ancestor.atEnd());
        Interior& interior = ancestor.asInterior();
        if (child == &interior.getEqualOrLess()) {
          typename Leaf::const_iterator it = leaf.begin();
          for (; it != leaf.end(); ++it) {
            ASSERT(!(interior.getExponent() <
              _conf.getExponent(it->get(), interior.getVar())));
            if (C::UseTreeDivMask) {
              ASSERT(interior.getDivMask().canDivide(it->getDivMask()));
            }
          }
        } else {
          ASSERT(child == &ancestor.getStrictlyGreater());
          typename Leaf::const_iterator it = leaf.begin();
          for (; it != leaf.end(); ++it) {
            ASSERT(interior.getExponent() <
              _conf.getExponent(it->get(), interior.getVar()));
            if (C::UseTreeDivMask) {
              ASSERT(interior.getDivMask().canDivide(it->getDivMask()));
            }
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
