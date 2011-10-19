#ifndef MATHIC_HEAP_GUARD
#define MATHIC_HEAP_GUARD

#include "stdinc.h"
#include "ComTree.h"
#include <vector>
#include <ostream>
#include <string>

namespace mathic {
  /** A heap priority queue.

	  Configuration serves the same role as for Geobucket. It must have these
	  fields that work as for Geobucket.

	  * A type Entry
	  * A type CompareResult
	  * A const or static method: CompareResult compare(Entry, Entry)
	  * A const or static method: bool cmpLessThan(CompareResult)
	  * A static const bool supportDeduplication
	  * A static or const method: bool cmpEqual(CompareResult)

	  It also has these additional fields:

	  * A static const bool fastIndex
	  If this field is true, then a faster way of calculating indexes is used.
	  This requires sizeof(Entry) to be a power of two! This can be achieved
	  by adding padding to Entry, but this class does not do that for you.
  */
  template<class C>
	class Heap {
  public:
	typedef C Configuration;
	typedef typename Configuration::Entry Entry;

  Heap(const Configuration& configuration): _conf(configuration) {}
	Configuration& getConfiguration() {return _conf;}
	const Configuration& getConfiguration() const {return _conf;}

	std::string getName() const;
	void push(Entry entry);
	void push(const Entry* begin, const Entry* end);
	Entry pop();
	Entry top() const {return _tree[Node()];}
	bool empty() const {return _tree.empty();}
	void print(std::ostream& out) const;

	void decreaseTop(Entry newEntry);

  private:
	typedef ComTree<Entry, Configuration::fastIndex> Tree;
	typedef typename Tree::Node Node;

	Node moveHoleDown(Node hole);
	void moveValueUp(Node pos, Entry value);

#ifdef MATHIC_DEBUG
	bool isValid() const;
#endif

	Tree _tree;
	Configuration _conf;
  };

  template<class C>
	std::string Heap<C>::getName() const {
	return std::string("heap(") +
	  (C::fastIndex ?  "fi" : "si") +
	  (C::supportDeduplication ? " dedup" : "") +
	  ')';
  }

  template<class C>
	void Heap<C>::push(Entry entry) {
	_tree.pushBack(entry);
	moveValueUp(_tree.lastLeaf(), entry);
	MATHIC_ASSERT(isValid());
  }

  template<class C>
	void Heap<C>::push(const Entry* begin, const Entry* end) {
	for (; begin != end; ++begin)
	  push(*begin);
  }

  template<class C>
	void Heap<C>::decreaseTop(Entry newEntry) {
	moveValueUp(moveHoleDown(Node()), newEntry);
	MATHIC_ASSERT(isValid());
  }

  template<class C>
	typename Heap<C>::Entry Heap<C>::pop() {
	Entry top = _tree[Node()];
	Entry movedValue = _tree[_tree.lastLeaf()];
	_tree.popBack();
	if (!_tree.empty())
	  moveValueUp(moveHoleDown(Node()), movedValue);
	return top;
	MATHIC_ASSERT(isValid());
  }

  template<class C>
	void Heap<C>::print(std::ostream& out) const {
	out << getName() << ": {" << _tree << "}\n";
  }

  template<class C>
	typename Heap<C>::Node Heap<C>::moveHoleDown(Node hole) {
	const Node firstWithout2Children = _tree.lastLeaf().next().parent();
	while (hole < firstWithout2Children) {
	  // can assume hole has two children here
	  Node child = hole.left();
	  Node sibling = child.next();
	  if (_conf.cmpLessThan(_conf.compare(_tree[child], _tree[sibling])))
		child = sibling;
	  _tree[hole] = _tree[child];
	  hole = child;
	}
	// if we are at a node that has a single left child
	if (hole == firstWithout2Children && _tree.lastLeaf().isLeft()) {
	  Node child = hole.left();
	  _tree[hole] = _tree[child];
	  hole = child;
	}
	return hole;
  }

  template<class C>
	void Heap<C>::moveValueUp(Node pos, Entry value) {
	const Node origPos = pos;
  again:
	while (!pos.isRoot()) {
	  const Node up = pos.parent();
	  typename C::CompareResult cmp = _conf.compare(_tree[up], value);
	  if (_conf.cmpLessThan(cmp)) {
		_tree[pos] = _tree[up];
		pos = up;
	  } else if (C::supportDeduplication && _conf.cmpEqual(cmp)) {
		_tree[up] = _conf.deduplicate(_tree[up], value);
		if (pos != origPos) {
		  // move elements back
		  Entry tmp = _tree[origPos];
		  for (Node p = origPos.parent(); p != pos; p = p.parent())
			std::swap(tmp, _tree[p]);
		}
		pos = origPos;
		value = _tree[_tree.lastLeaf()];
		_tree.popBack();
		if (origPos == _tree.lastLeaf().next()) {
		  MATHIC_ASSERT(isValid());
		  return;
		}
		goto again;
	  } else
		break;
	}
	_tree[pos] = value;
	MATHIC_ASSERT(isValid());
  }

#ifdef MATHIC_DEBUG
  template<class C>
	bool Heap<C>::isValid() const {
	MATHIC_ASSERT(_tree.isValid());
	for (Node i = Node().next(); i <= _tree.lastLeaf(); ++i) {
	  MATHIC_ASSERT(!_conf.cmpLessThan(_conf.compare(_tree[i.parent()], _tree[i])));
	}
	return true;
  }
#endif
}

#endif
