#ifndef K_D_TREE_LEAF_GUARD
#define K_D_TREE_LEAF_GUARD

#include "DivMask.h"
#include "Arena.h"
#include <algorithm>

/** A helper class for KDTree. A node in the tree. The ExtEntry
 comes from the KdTree. */
template<class Configuration, class ExtEntry>
class KDTreeNode;

/** A helper class for KDTree. An interior node in the tree. The ExtEntry
 comes from the KdTree. @todo: rename to KDTreeInternal. */
template<class Configuration, class ExtEntry>
class KDTreeInterior;

/** A helper class for KDTree. Represents a leaf in the tree. Leaves
 hold the monomials. The Configuration is as for KdTree. The ExtEntry
 comes from the KdTree. */
template<class Configuration, class ExtEntry>
class KDTreeLeaf;

namespace KDTreeHelpers {
  template<class C>
  class Comparer {
  public:
  Comparer(const C& conf): _conf(conf) {}
    template<class A, class B>
	bool operator()(const A& a, const B& b) const {
	  return _conf.isLessThan(a.get(), b.get());
	}
  private:
	const C& _conf;
  };
}

template<class EE, bool UseMask>
class Masked {
public:
  void resetDivMask() {ASSERT(false);}
  DivMask getDivMask() const {ASSERT(false); return DivMask();}
  template<class T>
  void updateMask(const T& entry) {}
};
template<class EE>
class Masked<EE, true> {
public:
  Masked() {resetDivMask();}
  void resetDivMask() {_mask = DivMask::getMaxMask();}
  DivMask& getDivMask() {return _mask;}
  const DivMask& getDivMask() const {return _mask;}
  template<class T>
  void updateMask(const T& t) {_mask.combineAnd(t.getDivMask());}

private:
  DivMask _mask;
};

template<class C, class EE>
class KDTreeNode : public Masked<EE, C::UseTreeDivMask> {
  typedef KDTreeNode<C, EE> Node;
  typedef KDTreeLeaf<C, EE> Leaf;
  typedef KDTreeInterior<C, EE> Interior;
  friend class KDTreeInterior<C, EE>;
  friend class KDTreeLeaf<C, EE>;
public:
  bool isLeaf() const {return _isLeaf;}
  const Leaf& asLeaf() const {
    ASSERT(isLeaf());
    return static_cast<const Leaf&>(*this);
  }
  Leaf& asLeaf() {
    ASSERT(isLeaf());
    return static_cast<Leaf&>(*this);
  }

  bool isInterior() const {return !isLeaf();}
  const Interior& asInterior() const {
    ASSERT(isInterior());
    return static_cast<Interior&>(*this);
  }
  Interior& asInterior() {
    ASSERT(isInterior());
    return static_cast<Interior&>(*this);
  }

  bool isEqualOrLessChild() const {
    return hasParent() && &getParent()->getEqualOrLess() == this;
  }
  bool isStrictlyGreaterChild() const {
    return hasParent() && &getParent()->getStrictlyGreater() == this;
  }

  bool hasParent() const {return _parent != 0;}
  Interior* getParent() {return _parent;}
  const Interior* getParent() const {return _parent;}

  /** Partitions [begin, end) into two parts. The
   returned node has the information about the split, while the returned
   iterator it is such that the equal-or-less part of the partition
   is [begin, it) and the striclty-greater part is [it, end). */
  template<class Iter>
  static std::pair<Interior*, Iter> preSplit
    (Interior* parent, Iter begin, Iter end, Arena& arena, const C& conf);

protected:
  KDTreeNode(bool leaf, Interior* parent):
    _isLeaf(leaf), _parent(parent) {}
  void setParent(Interior* parent) {_parent = parent;}

private:
  class SplitEqualOrLess;
  const bool _isLeaf;
  Interior* _parent;
};

template<class C, class EE>
class KDTreeInterior : public KDTreeNode<C, EE> {
public:
  typedef typename C::Exponent Exponent;
  typedef KDTreeInterior<C, EE> Interior;
  typedef KDTreeLeaf<C, EE> Leaf;
  typedef KDTreeNode<C, EE> Node;

  using Node::getParent;

  KDTreeInterior
    (Interior* parent,
     Node& equalOrLess,
     Node& strictlyGreater,
     size_t var,
     Exponent exponent):
   Node(false, parent),
   _equalOrLess(&equalOrLess),
   _strictlyGreater(&strictlyGreater),
   _var(var),
   _exponent(exponent) {
  }
  KDTreeInterior
    (Interior* parent,
     size_t var,
     Exponent exponent):
   Node(false, parent),
   _equalOrLess(0),
   _strictlyGreater(0),
   _var(var),
   _exponent(exponent) {
  }
  size_t getVar() const {return _var;}
  Exponent getExponent() const {return _exponent;}

  Node& getEqualOrLess() {return *_equalOrLess;}
  const Node& getEqualOrLess() const {return *_equalOrLess;}
  void setEqualOrLess(Node* equalOrLess) {_equalOrLess = equalOrLess;}

  Node& getStrictlyGreater() {return *_strictlyGreater;}
  const Node& getStrictlyGreater() const {return *_strictlyGreater;}
  void setStrictlyGreater(Node* strictlyGreater) {
    _strictlyGreater = strictlyGreater;
  }

  Node& getChildFor(const EE& entry, const C& conf) {
    if (getExponent() < conf.getExponent(entry.get(), getVar()))
      return getStrictlyGreater();
    else
      return getEqualOrLess();
  }

private:
  KDTreeInterior(const Interior&); // unavailable
  void operator=(const Interior&); // unavailable

  Node* _equalOrLess;
  Node* _strictlyGreater;
  size_t _var;
  Exponent _exponent;
};

template<class C, class EE>
class KDTreeLeaf : public KDTreeNode<C, EE> {
  typedef KDTreeInterior<C, EE> Interior;
  typedef KDTreeLeaf<C, EE> Leaf;
  typedef KDTreeNode<C, EE> Node;
 public:
  typedef typename C::Monomial Monomial;
  typedef EE* iterator;
  typedef const EE* const_iterator;
  typedef const EE& const_reference;
  typedef EE value_type;

  KDTreeLeaf(Arena& arena, const C& conf);
  KDTreeLeaf(Arena& arena, size_t capacity);

  void clear();

  typedef DivMask::Calculator<C> DivMaskCalculator;

  /** Copies [begin, end) into the new leaf. */
  template<class Iter>
  static Leaf* makeLeafCopy
    (Interior* parent, Iter begin, Iter end, Arena& arena,
    const DivMaskCalculator& calc, const C& conf);

  iterator begin() {return _begin;}
  const_iterator begin() const {return _begin;}
  iterator end() {return _end;}
  const_iterator end() const {return _end;}

  bool empty() const {return _begin == _end;}
  size_t size() const {return std::distance(_begin, _end);}
  EE& front() {ASSERT(!empty()); return *_begin;}
  const EE& front() const {ASSERT(!empty()); return *_begin;}
  EE& back() {ASSERT(!empty()); return *(_end - 1);}
  const EE& back() const {ASSERT(!empty()); return *(_end - 1);}

  void push_back(const EE& entry);
  void pop_back();
  void insert(iterator it, const EE& entry);

  void insert(const EE& entry, const C& conf);

  /** Returns how many were removed. */
  template<class EM>
  size_t removeMultiples(const EM& monomial, const C& conf);

  template<class EM>
  iterator findDivisor(const EM& extMonomial, const C& conf);
  template<class EM>
  const_iterator findDivisor(const EM& extMonomial, const C& conf) const {
    return const_cast<Leaf&>(*this).findDivisor(extMonomial, conf);
  }

  template<class EM, class DO>
  void findAllDivisors(const EM& extMonomial, DO& out, const C& conf);

  Interior& split(Arena& arena, const C& conf);

 private:
  KDTreeLeaf(const KDTreeLeaf& t); // unavailable
  void operator=(const KDTreeLeaf&); // unavailable

  iterator _begin;
  iterator _end;
#ifdef DEBUG
  const size_t _capacityDebug;
  const bool _sortOnInsertDebug;
#endif
};

template<class C, class EE>
KDTreeLeaf<C, EE>::KDTreeLeaf(Arena& arena, const C& conf):
Node(true, 0)
#ifdef DEBUG
,_capacityDebug(conf.getLeafSize())
,_sortOnInsertDebug(conf.getSortOnInsert())
#endif
{
  _begin = arena.allocArrayNoCon<EE>(conf.getLeafSize()).first;
  _end = _begin;
}

template<class C, class EE>
KDTreeLeaf<C,EE>::KDTreeLeaf(Arena& arena, size_t capacity):
Node(true, 0)
#ifdef DEBUG
,_capacityDebug(capacity)
,_sortOnInsertDebug(false)
#endif
{
  _begin = arena.allocArrayNoCon<EE>(capacity).first;
  _end = _begin;
}

template<class C, class EE>
void KDTreeLeaf<C, EE>::push_back(const EE& entry) {
  ASSERT(size() < _capacityDebug);
  new (_end) EE(entry);
  updateMask(entry);
  ++_end;
}

template<class C, class EE>
void KDTreeLeaf<C, EE>::pop_back() {
  ASSERT(!empty());
  --_end;
  _end->~EE();
}

template<class C, class EE>
void KDTreeLeaf<C, EE>::insert(iterator it, const EE& entry) {
  ASSERT(size() < _capacityDebug);
  if (it == end()) {
    push_back(entry);
    return;
  }
  push_back(back());
  iterator moveTo = end();
  for (--moveTo; moveTo != it; --moveTo)
    *moveTo = *(moveTo - 1);
  *it = entry;
}

template<class C, class EE>
void KDTreeLeaf<C, EE>::insert(const EE& entry, const C& conf) {
  ASSERT(size() < _capacityDebug);
  if (!conf.getSortOnInsert())
    push_back(entry);
  else {
    iterator it = std::upper_bound(begin(), end(), entry,
      KDTreeHelpers::Comparer<C>(conf));
	insert(it, entry);
  }
}

template<class C, class EE>
void KDTreeLeaf<C, EE>::clear() {while (!empty())
  while (!empty())
    pop_back();
}

template<class C, class EE>
class KDTreeNode<C, EE>::SplitEqualOrLess {
public:
  typedef typename C::Exponent Exponent;
  typedef typename C::Entry Entry;
  SplitEqualOrLess(size_t var, const Exponent& exp, const C& conf):
  _var(var), _exp(exp), _conf(conf) {}
  bool operator()(const Entry& entry) const {
    return !(_exp < _conf.getExponent(entry, _var));
  }
private:
  size_t _var;
  const Exponent& _exp;
  const C& _conf;
};

template<class C, class EE>
template<class Iter>
std::pair<KDTreeInterior<C, EE>*, Iter> KDTreeNode<C, EE>::preSplit
(Interior* parent, Iter begin, Iter end, Arena& arena, const C& conf) {
  ASSERT(begin != end);

  size_t var = (parent != 0 ? parent->getVar() : static_cast<size_t>(-1));
  while (true) {
    var = (var + 1) % conf.getVarCount();

    typename C::Exponent min = conf.getExponent(*begin, var);
    typename C::Exponent max = conf.getExponent(*begin, var);
    for (Iter it = begin; it != end; ++it) {
      min = std::min(min, conf.getExponent(*it, var));
      max = std::max(max, conf.getExponent(*it, var));
    }
    // todo: avoid infinite loop if all duplicates
    if (min == max)
      continue;
    // this formula for avg avoids overflow
    typename C::Exponent exp = min + (max - min) / 2;
    Interior* interior = new (arena.allocObjectNoCon<Interior>())
      Interior(parent, var, exp);
    SplitEqualOrLess cmp(var, exp, conf);
    Iter middle = std::partition(begin, end, cmp);
    return std::make_pair(interior, middle);
  }
}

template<class C, class EE>
template<class Iter>
KDTreeLeaf<C, EE>* KDTreeLeaf<C, EE>::makeLeafCopy 
(Interior* parent, Iter begin, Iter end,
Arena& arena, const DivMaskCalculator& calc, const C& conf) {
  ASSERT(static_cast<size_t>(std::distance(begin, end)) <= conf.getLeafSize());
  Leaf* leaf = new (arena.allocObjectNoCon<Leaf>()) Leaf(arena, conf);
  leaf->_parent = parent;
  // cannot directly copy as memory is not constructed.
  for (; begin != end; ++begin)
    leaf->push_back(EE(*begin, calc, conf));
  if (conf.getSortOnInsert())
    std::sort(leaf->begin(), leaf->end(), KDTreeHelpers::Comparer<C>(conf));
  return leaf;
}

template<class C, class EE>
template<class EM>
NO_PINLINE size_t KDTreeLeaf<C, EE>::removeMultiples(const EM& monomial, const C& conf) {
  iterator it = begin();
  iterator oldEnd = end();
  for (; it != oldEnd; ++it)
	if (monomial.divides(*it, conf))
	  break;
  if (it == oldEnd)
	return 0;
  iterator newEnd = it;
  for (++it; it != oldEnd; ++it) {
	if (!monomial.divides(*it, conf)) {
	  *newEnd = *it;
	  ++newEnd;
	}
  }
  // cannot just adjust _end as the superfluous
  // entries at the end need to be destructed.
  const size_t newSize = std::distance(begin(), newEnd);
  const size_t removedCount = size() - newSize;
  ASSERT(newSize < size());
  do {
    pop_back();
  } while (newSize < size());
  ASSERT(size() == newSize);
  return removedCount;
}

template<class C, class EE>
template<class EM>
NO_PINLINE typename KDTreeLeaf<C, EE>::iterator
KDTreeLeaf<C, EE>::findDivisor(const EM& extMonomial, const C& conf) {
  if (!conf.getSortOnInsert()) {
	const iterator stop = end();
	for (iterator it = begin(); it != stop; ++it)
	  if (it->divides(extMonomial, conf))
		return it;
	return stop;
  } else {
    iterator rangeEnd =
      std::upper_bound(begin(), end(), extMonomial,
      KDTreeHelpers::Comparer<C>(conf));
    iterator it = begin();
    for (; it != rangeEnd; ++it)
	  if (it->divides(extMonomial, conf))
	    return it;
	return end();
  }
}

template<class C, class EE>
template<class EM, class DO>
NO_PINLINE void KDTreeLeaf<C, EE>::
findAllDivisors(const EM& extMonomial, DO& out, const C& conf) {
  if (!conf.getSortOnInsert()) {
	const iterator stop = end();
	for (iterator it = begin(); it != stop; ++it)
	  if (it->divides(extMonomial, conf))
        out.push_back(it->get());
  } else {
    iterator rangeEnd =
      std::upper_bound(begin(), end(), extMonomial,
      KDTreeHelpers::Comparer<C>(conf));
    iterator it = begin();
    for (; it != rangeEnd; ++it)
	  if (it->divides(extMonomial, conf))
        out.push_back(it->get());
  }  
}

template<class C, class EE>
struct ExpOrder {
  typedef typename C::Entry Entry;
  ExpOrder(size_t var, const C& conf): _var(var), _conf(conf) {}
  bool operator()(const EE& a, const EE& b) const {
    return _conf.getExponent(a.get(), _var) < _conf.getExponent(b.get(), _var);
  }
private:
  const size_t _var;
  const C& _conf;
};

template<class C, class EE>
NO_PINLINE KDTreeInterior<C, EE>&
KDTreeLeaf<C, EE>::split(Arena& arena, const C& conf) {
  ASSERT(conf.getVarCount() > 0);
  ASSERT(size() >= 2);
  // ASSERT not all equal
  Leaf& other = *new (arena.allocObjectNoCon<Leaf>()) Leaf(arena, conf);

  typename C::Exponent exp;
  size_t var = Node::hasParent() ?
    Node::getParent()->getVar() : static_cast<size_t>(-1);
  while (true) {
    var = (var + 1) % conf.getVarCount();
    
    if (1) {
      typename C::Exponent min = conf.getExponent(front().get(), var);
      typename C::Exponent max = conf.getExponent(front().get(), var);
      for (iterator it = begin(); it != end(); ++it) {
        min = std::min(min, conf.getExponent(it->get(), var));
        max = std::max(max, conf.getExponent(it->get(), var));
      }
      if (min == max && size() > 1) {
        // todo: avoid infinite loop if all equal
        continue;
      }
      exp = min + (max - min) / 2; // this formula for avg avoids overflow

      iterator newEnd = begin();
      for (iterator it = begin(); it != end(); ++it) {
        if (exp < conf.getExponent(it->get(), var))
          other.push_back(*it);
        else {
          if (it != newEnd)
            *newEnd = *it;
          ++newEnd;
        }
      }
      while (newEnd != end())
        pop_back();
    } else {
      iterator middle = begin() + size() / 2;
      ExpOrder<C, EE> order(var, conf);

      std::nth_element(begin(), middle, end(), order);
      if (middle != end()) {
        exp = conf.getExponent(middle->get(), var);
        while (middle != end() && conf.getExponent(middle->get(), var) == exp)
          ++middle;
      }
      if (middle == end() && size() > 1) {
        // todoL avoid infinite loop if all equal.
        continue; // bad split, use another variable
      }
      ASSERT(middle != end());
      ASSERT(exp != conf.getExponent(middle->get(), var));

#ifdef DEBUG
      for (iterator it = begin(); it != middle; ++it) {
        ASSERT(!(exp < conf.getExponent(it->get(), var)));
      }
      for (iterator it = middle; it != end(); ++it) {
        ASSERT(!(conf.getExponent(it->get(), var) < exp));
      }
#endif
      // nth_element does not guarantee where equal elements go,
      // so we cannot just copy [middle, end()).
      iterator newEnd = begin();
      for (iterator it = begin(); it != end(); ++it) {
        if (exp < conf.getExponent(it->get(), var))
          other.push_back(*it);
        else {
          if (it != newEnd)
            *newEnd = *it;
          ++newEnd;
        }
      }
      while (newEnd != end())
        pop_back();

    }
    ASSERT(other.size() < conf.getLeafSize());
    ASSERT(size() < conf.getLeafSize());
    break;
  }

  Interior& interior = *new (arena.allocObjectNoCon<Interior>())
    Interior(Node::getParent(), *this, other, var, exp);

  if (Node::hasParent()) {
    if (Node::isEqualOrLessChild())
      Node::getParent()->setEqualOrLess(&interior);
    else {
      ASSERT(Node::isStrictlyGreaterChild());
      Node::getParent()->setStrictlyGreater(&interior);
    }
  }
  setParent(&interior);
  other.setParent(&interior);
  if (C::UseTreeDivMask) {
    KDTreeNode<C,EE>::resetDivMask();
    for (const_iterator it = begin(); it != end(); ++it)
      updateMask(*it);
    interior.updateMask(*this);
    interior.updateMask(other);
  }

  return interior;
}


#endif
