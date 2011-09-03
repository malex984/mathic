#ifndef DIV_ARRAY_GUARD
#define DIV_ARRAY_GUARD

#include <vector>
#include <string>
#include <list>

/** An object that supports queries for divisors of a monomial using
 an array of monomials. See DivFinder for more documentation.
*/
template<class Configuration, bool UseLinkedList = false>
class DivList;

// implementation details for DivList.
namespace DivListHelper {
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

  void removeMultiples(const Monomial& monomial);
  void insert(const Entry& entry);
  void insertReminimize(const Entry& entry);
  iterator findDivisor(const Monomial& monomial);
  const_iterator findDivisor(const Monomial& monomial) const;

  iterator begin() {return _list.begin();}
  const_iterator begin() const {return _list.begin();}
  iterator end() {return _list.end();}
  const_iterator end() const {return _list.end();}

  std::string getName() const;

  C& getConfiguration() {return _conf;}
  const C& getConfiguration() const {return _conf;}

 private:
  List _list;
  C _conf;
};

namespace DivListHelper {
  template<class C, class E, class M>
  void removeMultiples(C& conf, std::vector<E>& list, const M& monomial) {
    typedef typename std::vector<E>::iterator iterator;
	iterator it = list.begin();
	iterator oldEnd = list.end();
	for (; it != oldEnd; ++it)
	  if (conf.divides(monomial, *it))
		break;
	if (it == oldEnd)
	  return;
	iterator newEnd = it;
	for (++it; it != oldEnd; ++it) {
	  if (!conf.divides(monomial, *it)) {
		*newEnd = *it;
		++newEnd;
	  }
	}
	list.resize(std::distance(list.begin(), newEnd));
  }

  template<class C, class E, class M>
  void removeMultiples(C& conf, std::list<E>& list, const M& monomial) {
    typedef typename std::list<E>::iterator iterator;
	iterator it = list.begin();
	iterator oldEnd = list.end();
	for (; it != oldEnd; ++it)
	  if (conf.divides(monomial, *it))
		break;
	if (it == oldEnd)
	  return;
	iterator newEnd = it;
	for (++it; it != oldEnd; ++it) {
	  if (!conf.divides(monomial, *it)) {
		*newEnd = *it;
		++newEnd;
	  }
	}
	list.resize(std::distance(list.begin(), newEnd));
  }
}

template<class C, bool ULL>
  void DivList<C, ULL>::insert(const Entry& entry) {
  _list.push_back(entry);
}

template<class C, bool ULL>
  void DivList<C, ULL>::insertReminimize(const Entry& entry) {
  removeMultiples(entry);
  insert(entry);
}

template<class C, bool ULL>
void DivList<C, ULL>::removeMultiples(const Monomial& monomial) {
  DivListHelper::removeMultiples(_conf, _list, monomial);
}

template<class C, bool ULL>
typename DivList<C, ULL>::iterator
DivList<C, ULL>::findDivisor(const Monomial& monomial) {
  const iterator stop = end();
  for (iterator it = begin(); it != stop; ++it)
    if (_conf.divides(*it, monomial))
      return it;
  return stop;
}

template<class C, bool ULL>
typename DivList<C, ULL>::const_iterator
DivList<C, ULL>::findDivisor(const Monomial& monomial) const {
  return const_cast<DivList<C>&>(*this).findDivisor(monomial);
}

template<class C, bool ULL>
std::string DivList<C, ULL>::getName() const {
  return std::string("DivList") +
    (ULL ? " linked" : " array");
}

#endif
