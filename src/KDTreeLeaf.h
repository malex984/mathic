#ifndef K_D_TREE_LEAF_GUARD
#define K_D_TREE_LEAF_GUARD

#include "Arena.h"

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

protected:
  KDTreeNode(bool leaf, Interior* parent):
    _isLeaf(leaf), _parent(parent) {}
  void setParent(Interior* parent) {_parent = parent;}

private:
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

  KDTreeLeaf(Arena& arena, const C& conf);
  ~KDTreeLeaf();

  iterator begin() {return _begin;}
  const_iterator begin() const {return _begin;}
  iterator end() {return _end;}
  const_iterator end() const {return _end;}

  bool empty() const {return _begin == _end;}
  size_t size() const {return _end - _begin;}
  Entry& front() {ASSERT(!empty()); return *_begin;}
  const Entry& front() const {ASSERT(!empty()); return *_begin;}
  Entry& back() {ASSERT(!empty()); return *(_end - 1);}
  const Entry& back() const {ASSERT(!empty()); return *(_end - 1);}

  void push_back(const Entry& entry);
  void pop_back();
  void insert(iterator it, const Entry& entry);

  void insert(const Entry& entry, const C& conf);
  bool removeMultiples(const Monomial& monomial, const C& conf);

  iterator findDivisor(const Monomial& monomial, const C& conf);
  const_iterator findDivisor(const Monomial& monomial, const C& conf) const {
    return const_cast<KDTreeLeaf<C>&>(*this).findDivisor(monomial, conf);
  }

  Interior& split(Arena& arena, const C& conf);

 private:
  KDTreeLeaf(const KDTreeLeaf& t); // unavailable
  void operator=(const KDTreeLeaf&) {ASSERT(false);} // unavailable

  iterator _begin;
  iterator _end;
  IF_DEBUG(size_t _capacityDebug;)
};

template<class C>
KDTreeLeaf<C>::KDTreeLeaf(Arena& arena, const C& conf):
 Node(true, 0), _begin(0), _end(0) {
  IF_DEBUG(_capacityDebug = 0);

  IF_DEBUG(_capacityDebug = conf.getLeafSize());
  _begin = arena.allocArrayNoCon<Entry>(conf.getLeafSize()).first;
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
  for (iterator moveTo = end() - 1; moveTo != it; --moveTo)
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
bool KDTreeLeaf<C>::removeMultiples(const Monomial& monomial, const C& conf) {
  iterator it = begin();
  iterator oldEnd = end();
  for (; it != oldEnd; ++it)
	if (conf.divides(monomial, *it))
	  break;
  if (it == oldEnd)
	return false;
  iterator newEnd = it;
  for (++it; it != oldEnd; ++it) {
	if (!conf.divides(monomial, *it)) {
	  *newEnd = *it;
	  ++newEnd;
	}
  }
  const size_t newSize = newEnd - begin();
  ASSERT(newSize < size());
  do {
    pop_back();
  } while (newSize < size());
  ASSERT(size() == newSize);
  return true;
}

template<class C>
typename KDTreeLeaf<C>::iterator
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
KDTreeInterior<C>& KDTreeLeaf<C>::split(Arena& arena, const C& conf) {
  ASSERT(conf.getVarCount() > 0);
  ASSERT(size() >= 2);
  // ASSERT not all equal
    Leaf& other = *new (arena.allocObjectNoCon<Leaf>()) Leaf(arena, conf);

    size_t var;
    typename C::Exponent exp;
    while (true) {
      var = rand() % conf.getVarCount();
      typename C::Exponent min = conf.getExponent(front(), var);
      typename C::Exponent max = conf.getExponent(front(), var);
      for (iterator it = begin(); it != end(); ++it) {
        min = std::min(min, conf.getExponent(*it, var));
        max = std::max(max, conf.getExponent(*it, var));
      }
      if (min == max) {
        if (back() == front())
          pop_back();
        else
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
