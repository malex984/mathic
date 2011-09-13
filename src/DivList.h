#ifndef DIV_ARRAY_GUARD
#define DIV_ARRAY_GUARD

#include <vector>
#include <string>
#include <list>
#include <algorithm>

/** An object that supports queries for divisors of a monomial using
 an array of monomials. See DivFinder for more documentation.

 Extra fields for Configuration:

 * bool getSortOnInsert() const
  Keep the monomials sorted to speed up queries.
*/
template<class Configuration, bool UseLinkedList = false>
class DivList;

namespace DivListHelper {
  // Implementation details for DivList.

  template<bool B, class E>
  struct ListImpl;

  template<class Entry>
  struct ListImpl<false, Entry> {
    typedef std::vector<Entry> Impl;
  };
  template<class Entry>
  struct ListImpl<true, Entry> {
    typedef std::list<Entry> Impl;
  };
}

template<class C, bool ULL>
class DivList {
  typedef typename DivListHelper::ListImpl<ULL, typename C::Entry>::Impl
    List;

 public:
  static const bool UseLinkedList = ULL;

  typedef typename C::Monomial Monomial;
  typedef typename C::Entry Entry;
  typedef typename C::Exponent Exponent;

  typedef typename List::const_iterator const_iterator;
  typedef typename List::iterator iterator;

  DivList(const C& configuration): _conf(configuration) {}

  bool removeMultiples(const Monomial& monomial);
  void insert(const Entry& entry);
  iterator findDivisor(const Monomial& monomial);
  const_iterator findDivisor(const Monomial& monomial) const;

  template<class DO>
  void findAllDivisors(const Monomial& monomial, DO& out);
  template<class DO>
  void findAllDivisors(const Monomial& monomial, DO& out) const;

  iterator begin() {return _list.begin();}
  const_iterator begin() const {return _list.begin();}
  iterator end() {return _list.end();}
  const_iterator end() const {return _list.end();}
  size_t size() const {return _list.size();}

  std::string getName() const;

  C& getConfiguration() {return _conf;}
  const C& getConfiguration() const {return _conf;}

  void moveToFront(iterator pos);

 private:
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

  List _list;
  C _conf;
};

namespace DivListHelper {
  template<class C, class E, class M>
  bool removeMultiples(C& conf, std::vector<E>& list, const M& monomial) {
    typedef typename std::vector<E>::iterator iterator;
	iterator it = list.begin();
	iterator oldEnd = list.end();
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
	const size_t newSize = newEnd - list.begin();
    ASSERT(newSize < list.size());
	list.resize(newSize);
    return true;
  }

  template<class C, class E, class M>
  bool removeMultiples(C& conf, std::list<E>& list, const M& monomial) {
    typedef typename std::list<E>::iterator iterator;
	iterator it = list.begin();
	iterator oldEnd = list.end();
    bool removedSome = false;
    while (it != oldEnd) {
	  if (conf.divides(monomial, *it)) {
		removedSome = true;
		it = list.erase(it);
	  } else
		++it;
	}
	return removedSome;
  }

  template<class E, class It>
  void moveToFront(std::vector<E>& list, It pos) {
    E valueToMove = *pos;
    It begin = list.begin();
    while (pos != begin) {
      It prev = pos;
	  --pos;
      *prev = *pos;
	}
    list.front() = valueToMove;
  }

  template<class E, class It>
  void moveToFront(std::list<E>& list, It pos) {
    list.splice(list.begin(), list, pos);
  }

  template<class C, class E>
  typename std::list<E>::iterator
  insertSort(C& conf, std::list<E>& list, const E& entry) {
    typedef typename std::list<E>::iterator iterator;
    iterator end = list.end();
    iterator it = list.begin();
    for (; it != end; ++it)
	  if (conf.isLessThan(entry, *it))
		break;
	return list.insert(it, entry);
  }

  template<class C, class E>
  typename std::vector<E>::iterator
  insertSort(C& conf, std::vector<E>& list, const E& entry) {
    typedef typename std::vector<E>::iterator iterator;
    iterator it = std::upper_bound(list.begin(), list.end(), entry,
      conf.getComparer());
	return list.insert(it, entry);
  }

  template<class C, class E, class M>
  typename std::vector<E>::iterator  
  findDivisorSorted(C& conf, std::vector<E>& list, const M& monomial) {
    typedef typename std::vector<E>::iterator iterator;
    iterator rangeEnd =
      std::upper_bound(list.begin(), list.end(), monomial,
        conf.getComparer());
    iterator it = list.begin();
    for (; it != rangeEnd; ++it)
      if (conf.divides(*it, monomial))
	    return it;
	return list.end();
  }

  template<class C, class E, class M>
  typename std::list<E>::iterator  
  findDivisorSorted(C& conf, std::list<E>& list, const M& monomial) {
    typedef typename std::list<E>::iterator iterator;
    iterator end = list.end();
    iterator it = list.begin();
    size_t count = 0;
    for (; it != end; ++it) {
	  ++count;
      if (count == 35) {
		count = 0;
        if (conf.isLessThan(monomial, *it))
		  break;
	  }
      if (conf.divides(*it, monomial))
	    return it;
    }
	return end;
  }

  template<class C, class E, class M, class DO>
  void findAllDivisorsSorted
   (C& conf, std::vector<E>& list, const M& monomial, DO& out) {
    typedef typename std::vector<E>::iterator iterator;
    iterator rangeEnd =
      std::upper_bound(list.begin(), list.end(), monomial, conf.getComparer());
    iterator it = list.begin();
    for (; it != rangeEnd; ++it)
      if (conf.divides(*it, monomial))
        out.push_back(*it);
  }

  template<class C, class E, class M, class O>
  void findAllDivisorsSorted
   (C& conf, std::list<E>& list, const M& monomial, O& out) {
    typedef typename std::list<E>::iterator iterator;
    iterator end = list.end();
    iterator it = list.begin();
    size_t count = 0;
    for (; it != end; ++it) {
	  ++count;
      if (count == 35) {
		count = 0;
        if (conf.isLessThan(monomial, *it))
		  break;
	  }
      if (conf.divides(*it, monomial))
        out.push_back(*it);
    }
  }
}

template<class C, bool ULL>
void DivList<C, ULL>::insert(const Entry& entry) {
  if (!_conf.getSortOnInsert())
    _list.push_back(entry);
  else
	DivListHelper::insertSort(_conf, _list, entry);
}

template<class C, bool ULL>
bool DivList<C, ULL>::removeMultiples(const Monomial& monomial) {
  return DivListHelper::removeMultiples(_conf, _list, monomial);
}

template<class C, bool ULL>
typename DivList<C, ULL>::iterator
DivList<C, ULL>::findDivisor(const Monomial& monomial) {
  if (!_conf.getSortOnInsert()) {
	const iterator stop = end();
	for (iterator it = begin(); it != stop; ++it)
	  if (_conf.divides(*it, monomial))
		return it;
	return stop;
  } else
	return DivListHelper::findDivisorSorted(_conf, _list, monomial);
}

template<class C, bool ULL>
typename DivList<C, ULL>::const_iterator
DivList<C, ULL>::findDivisor(const Monomial& monomial) const {
  return const_cast<DivList<C>&>(*this).findDivisor(monomial);
}

template<class C, bool ULL>
template<class DO>
void DivList<C, ULL>::
findAllDivisors(const Monomial& monomial, DO& out) {
  if (!_conf.getSortOnInsert()) {
	const iterator stop = end();
	for (iterator it = begin(); it != stop; ++it)
	  if (_conf.divides(*it, monomial))
        out.push_back(*it);
  } else
    DivListHelper::findAllDivisorsSorted(_conf, _list, monomial, out);
}

template<class C, bool ULL>
template<class DO>
void DivList<C, ULL>::findAllDivisors(const Monomial& monomial, DO& output) const {
  ConstDivisorOutput<DO> constOutput(output);
  const_cast<DivList<C, ULL>&>(*this).findAllDivisors(monomial, constOutput);
}

template<class C, bool ULL>
std::string DivList<C, ULL>::getName() const {
  return std::string("DivList") +
	(_conf.getSortOnInsert() ? " sort" : "") +
    (ULL ? " linked" : " array");
}

template<class C, bool ULL>
void DivList<C, ULL>::moveToFront(iterator pos) {
  DivListHelper::moveToFront(_list, pos);
}

#endif
