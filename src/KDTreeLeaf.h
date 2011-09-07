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

  Interior* getParent() {return _parent;}
  const Interior* getParent() const {return _parent;}

protected:
  KDTreeNode(bool isLeaf, Interior* parent):
    _isLeaf(isLeaf), _parent(parent) {}
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

  KDTreeInterior(): Node(false, 0) {}
  void reset(Leaf& toSplit, Leaf& other, Arena& arena, const C& conf) {
    ASSERT(conf.getVarCount() > 0);
    //_var = (parent == 0 ? 0 : parent->getVar() % conf.getVarCount());

    _equalOrLess = &toSplit;
    _strictlyGreater = &other;

    setParent(toSplit.getParent());
    if (getParent() != 0) {
      if (&getParent()->getEqualOrLess() == &toSplit)
        getParent()->setEqualOrLess(this);
      else {
        ASSERT(&getParent()->getStrictlyGreater() == &toSplit);
        getParent()->setStrictlyGreater(this);
      }
    }
    toSplit.setParent(this);
    other.reset(arena, conf);
    other.setParent(this);
    size_t toMove = toSplit.size() / 2;
    while (other.size() < toMove) {
      other.push_back(toSplit.back());
      toSplit.pop_back();
    }
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

private:
  KDTreeInterior(const Interior&); // unavailable
  void operator=(const Interior&); // unavailable

  size_t _var;
  Exponent _exponent;
  KDTreeNode<C>* _equalOrLess;
  KDTreeNode<C>* _strictlyGreater;
};

template<class C>
class KDTreeLeaf : public KDTreeNode<C> {
  typedef std::vector<typename C::Entry> List;
 public:
  typedef typename C::Monomial Monomial;
  typedef typename C::Entry Entry;
  typedef Entry* iterator;
  typedef const Entry* const_iterator;
  typedef KDTreeNode<C> Node;

  KDTreeLeaf();
  ~KDTreeLeaf();
  void reset(Arena& arena, const C& conf);

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

  KDTreeLeaf(const KDTreeLeaf& t): Node(t), _begin(0), _end(0) {
    IF_DEBUG(_capacityDebug = 0);
    ASSERT(t._begin == 0);
    ASSERT(t._end == 0);
    ASSERT(t._capacityDebug == 0);
  } // unavailable
  void operator=(const KDTreeLeaf&) {ASSERT(false);} // unavailable

 private:
  iterator _begin;
  iterator _end;
  IF_DEBUG(size_t _capacityDebug;)
};

template<class C>
KDTreeLeaf<C>::KDTreeLeaf(): Node(true, 0), _begin(0), _end(0) {
  IF_DEBUG(_capacityDebug = 0);
}

template<class C>
KDTreeLeaf<C>::~KDTreeLeaf() {
  while (!empty())
    pop_back();
}

template<class C>
void KDTreeLeaf<C>::reset(Arena& arena, const C& conf) {
  IF_DEBUG(_capacityDebug = conf.getLeafSize());
  _begin = arena.allocArrayNoCon<Entry>(conf.getLeafSize()).first;
  _end = _begin;
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

#endif
