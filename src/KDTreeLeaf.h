#ifndef K_D_TREE_LEAF_GUARD
#define K_D_TREE_LEAF_GUARD

#include "Arena.h"

/** A helper class for KDTree. Represents a leaf in the tree. Leaves
 hold the actual monomials. The Configuration is as for KdTree. */
template<class Configuration>
class KDTreeLeaf;

template<class C>
class KDTreeLeaf {
  typedef std::vector<typename C::Entry> List;
 public:
  typedef typename C::Monomial Monomial;
  typedef typename C::Entry Entry;
  typedef Entry* iterator;
  typedef const Entry* const_iterator;

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

  KDTreeLeaf(const KDTreeLeaf& t): _begin(0), _end(0) {
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
KDTreeLeaf<C>::KDTreeLeaf():
 _begin(0), _end(0) {
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
