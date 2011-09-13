#ifndef K_D_TREE_LEAF_GUARD
#define K_D_TREE_LEAF_GUARD

#include "Arena.h"
#include <algorithm>

/** A helper class for KDTree. A node in the tree. */
template<class Configuration>
class KDTreeNode;

/** A helper class for KDTree. An interior node in the tree. */
template<class Configuration>
class KDTreeInterior;

/** A helper class for KDTree. Represents a leaf in the tree. Leaves
 hold the monomials. The Configuration is as for KdTree. */
template<class Configuration>
class KDTreeLeaf;

template<class C>
class KDTreeNode {
  typedef KDTreeNode<C> Node;
  typedef KDTreeLeaf<C> Leaf;
  typedef KDTreeInterior<C> Interior;
  friend class KDTreeInterior<C>;
  friend class KDTreeLeaf<C>;
public:
  virtual ~KDTreeNode() {}

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

  const bool _isLeaf; // todo: check performance of virtual replacement
  Interior* _parent;
};

template<class C>
class KDTreeInterior : public KDTreeNode<C> {
public:
  typedef typename C::Exponent Exponent;
  typedef KDTreeInterior<C> Interior;
  typedef KDTreeLeaf<C> Leaf;
  typedef KDTreeNode<C> Node;

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
  ~KDTreeInterior();
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

  Node& getChildFor(const typename C::Entry& entry, const C& conf) {
    if (getExponent() < conf.getExponent(entry, getVar()))
      return getStrictlyGreater();
    else
      return getEqualOrLess();
  }

private:
  KDTreeInterior(const Interior&); // unavailable
  void operator=(const Interior&); // unavailable

  KDTreeNode<C>* _equalOrLess;
  KDTreeNode<C>* _strictlyGreater;
  size_t _var;
  Exponent _exponent;
};

template<class C>
class KDTreeLeaf : public KDTreeNode<C> {
  typedef std::vector<typename C::Entry> List;
  typedef KDTreeInterior<C> Interior;
  typedef KDTreeLeaf<C> Leaf;
  typedef KDTreeNode<C> Node;
 public:
  typedef typename C::Monomial Monomial;
  typedef typename C::Entry Entry;
  typedef Entry* iterator;
  typedef const Entry* const_iterator;
  typedef const Entry& const_reference;
  typedef Entry value_type;

  KDTreeLeaf(Arena& arena, const C& conf);
  KDTreeLeaf(Arena& arena, size_t capacity);
  ~KDTreeLeaf();

  void clear();

  /** Copies [begin, end) into the new leaf. */
  template<class Iter>
  static Leaf* makeLeafCopy
    (Interior* parent, Iter begin, Iter end, Arena& arena, const C& conf);

  iterator begin() {return _begin;}
  const_iterator begin() const {return _begin;}
  iterator end() {return _end;}
  const_iterator end() const {return _end;}

  bool empty() const {return _begin == _end;}
  size_t size() const {return std::distance(_begin, _end);}
  Entry& front() {ASSERT(!empty()); return *_begin;}
  const Entry& front() const {ASSERT(!empty()); return *_begin;}
  Entry& back() {ASSERT(!empty()); return *(_end - 1);}
  const Entry& back() const {ASSERT(!empty()); return *(_end - 1);}

  void push_back(const Entry& entry);
  void pop_back();
  void insert(iterator it, const Entry& entry);

  void insert(const Entry& entry, const C& conf);

  /** Returns how many were removed. */
  size_t removeMultiples(const Monomial& monomial, const C& conf);

  iterator findDivisor(const Monomial& monomial, const C& conf);
  const_iterator findDivisor(const Monomial& monomial, const C& conf) const {
    return const_cast<KDTreeLeaf<C>&>(*this).findDivisor(monomial, conf);
  }

  template<class DO>
  void findAllDivisors(const Monomial& monomial, DO& out, const C& conf);

  Interior& split(Arena& arena, const C& conf);

 private:
  KDTreeLeaf(const KDTreeLeaf& t); // unavailable
  void operator=(const KDTreeLeaf&); // unavailable

  iterator _begin;
  iterator _end;
  IF_DEBUG(const size_t _capacityDebug;)
  IF_DEBUG(const bool _sortOnInsertDebug;)
};

template<class C>
KDTreeInterior<C>::~KDTreeInterior() {
  if (_equalOrLess != 0)
    _equalOrLess->~KDTreeNode();
  if (_strictlyGreater != 0)
    _strictlyGreater->~KDTreeNode();
}

template<class C>
KDTreeLeaf<C>::KDTreeLeaf(Arena& arena, const C& conf):
Node(true, 0)
#ifdef DEBUG
,_capacityDebug(conf.getLeafSize())
,_sortOnInsertDebug(conf.getSortOnInsert())
#endif
{
  _begin = arena.allocArrayNoCon<Entry>(conf.getLeafSize()).first;
  _end = _begin;
}

template<class C>
KDTreeLeaf<C>::KDTreeLeaf(Arena& arena, size_t capacity):
Node(true, 0)
#ifdef DEBUG
,_capacityDebug(capacity)
,_sortOnInsertDebug(false)
#endif
{
  _begin = arena.allocArrayNoCon<Entry>(capacity).first;
  _end = _begin;
}

template<class C>
KDTreeLeaf<C>::~KDTreeLeaf() {
  while (!empty())
    pop_back();
}

template<class C>
void KDTreeLeaf<C>::push_back(const Entry& entry) {
  ASSERT(size() < _capacityDebug);
  new (_end) Entry(entry);
  ++_end;
}

template<class C>
void KDTreeLeaf<C>::pop_back() {
  ASSERT(!empty());
  --_end;
  _end->~Entry();
}

template<class C>
void KDTreeLeaf<C>::insert(iterator it, const Entry& entry) {
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

template<class C>
void KDTreeLeaf<C>::insert(const Entry& entry, const C& conf) {
  ASSERT(size() < _capacityDebug);
  if (!conf.getSortOnInsert())
    push_back(entry);
  else {
    iterator it = std::upper_bound(begin(), end(), entry, conf.getComparer());
	insert(it, entry);
  }
}

template<class C>
void KDTreeLeaf<C>::clear() {while (!empty())
  while (!empty())
    pop_back();
}

template<class C>
class KDTreeNode<C>::SplitEqualOrLess {
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

template<class C>
template<class Iter>
std::pair<KDTreeInterior<C>*, Iter> KDTreeNode<C>::preSplit
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

template<class C>
template<class Iter>
KDTreeLeaf<C>* KDTreeLeaf<C>::makeLeafCopy 
(Interior* parent, Iter begin, Iter end, Arena& arena, const C& conf) {
  ASSERT(static_cast<size_t>(distance(begin, end)) <= conf.getLeafSize());
  KDTreeLeaf<C>* leaf = new (arena.allocObjectNoCon<KDTreeLeaf<C> >())
    KDTreeLeaf<C>(arena, conf);
  leaf->_parent = parent;
  // cannot directly copy as memory is not constructed.
  for (; begin != end; ++begin)
    leaf->push_back(*begin);
  if (conf.getSortOnInsert())
    std::sort(leaf->begin(), leaf->end(), conf.getComparer());    
  return leaf;
}

template<class C>
size_t KDTreeLeaf<C>::removeMultiples(const Monomial& monomial, const C& conf) {
  iterator it = begin();
  iterator oldEnd = end();
  for (; it != oldEnd; ++it)
	if (conf.divides(monomial, *it))
	  break;
  if (it == oldEnd)
	return 0;
  iterator newEnd = it;
  for (++it; it != oldEnd; ++it) {
	if (!conf.divides(monomial, *it)) {
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

template<class C>
NO_PINLINE typename KDTreeLeaf<C>::iterator
KDTreeLeaf<C>::findDivisor(const Monomial& monomial, const C& conf) {
  if (!conf.getSortOnInsert()) {
	const iterator stop = end();
	for (iterator it = begin(); it != stop; ++it)
	  if (conf.divides(*it, monomial))
		return it;
	return stop;
  } else {
    iterator rangeEnd =
      std::upper_bound(begin(), end(), monomial, conf.getComparer());
    iterator it = begin();
    for (; it != rangeEnd; ++it)
      if (conf.divides(*it, monomial))
	    return it;
	return end();
  }
}

template<class C>
template<class DO>
NO_PINLINE void KDTreeLeaf<C>::
findAllDivisors(const Monomial& monomial, DO& out, const C& conf) {
  if (!conf.getSortOnInsert()) {
	const iterator stop = end();
	for (iterator it = begin(); it != stop; ++it)
	  if (conf.divides(*it, monomial))
        out.push_back(*it);
  } else {
    iterator rangeEnd =
      std::upper_bound(begin(), end(), monomial, conf.getComparer());
    iterator it = begin();
    for (; it != rangeEnd; ++it)
      if (conf.divides(*it, monomial))
        out.push_back(*it);
  }  
}

template<class C>
struct ExpOrder {
  typedef typename C::Entry Entry;
  ExpOrder(size_t var, const C& conf): _var(var), _conf(conf) {}
  bool operator()(const Entry& a, const Entry& b) const {
    return _conf.getExponent(a, _var) < _conf.getExponent(b, _var);
  }
private:
  const size_t _var;
  const C& _conf;
};

template<class C>
NO_PINLINE KDTreeInterior<C>& KDTreeLeaf<C>::split(Arena& arena, const C& conf) {
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
      typename C::Exponent min = conf.getExponent(front(), var);
      typename C::Exponent max = conf.getExponent(front(), var);
      for (iterator it = begin(); it != end(); ++it) {
        min = std::min(min, conf.getExponent(*it, var));
        max = std::max(max, conf.getExponent(*it, var));
      }
      if (min == max && size() > 1) {
        if (back() == front()) {
          // no good way to split a leaf of all duplicates other than to
          // detect and remove them.
          while (back() == front() && size() > 1)
            pop_back();
        }
        continue;
      }
      exp = min + (max - min) / 2; // this formula for avg avoids overflow

      iterator newEnd = begin();
      for (iterator it = begin(); it != end(); ++it) {
        if (exp < conf.getExponent(*it, var))
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
      ExpOrder<C> order(var, conf);

      std::nth_element(begin(), middle, end(), order);
      if (middle != end()) {
        exp = conf.getExponent(*middle, var);
        while (middle != end() && conf.getExponent(*middle, var) == exp)
          ++middle;
      }
      if (middle == end() && size() > 1) {
        if (back() == front()) {
          // no good way to split a leaf of all duplicates other than to
          // detect and remove them.
          while (back() == front())
            pop_back();
        }
        continue; // bad split, use another variable
      }
      ASSERT(middle != end());
      ASSERT(exp != conf.getExponent(*middle, var));

#ifdef DEBUG
      for (iterator it = begin(); it != middle; ++it) {
        ASSERT(!(exp < conf.getExponent(*it, var)));
      }
      for (iterator it = middle; it != end(); ++it) {
        ASSERT(!(conf.getExponent(*it, var) < exp));
      }
#endif
      // nth_element does not guarantee where equal elements go,
      // so we cannot just copy [middle, end()).
      iterator newEnd = begin();
      for (iterator it = begin(); it != end(); ++it) {
        if (exp < conf.getExponent(*it, var))
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

  return interior;
}


#endif
