#ifndef K_D_TREE_WALKER_GUARD
#define K_D_TREE_WALKER_GUARD

#include "KDTreeLeaf.h"

/** A helper class for KDTree. Encapsulates navigation in the tree.
 ExtEntry is from the kd tree. */
template<class Configuration, class ExtEntry>
class KDTreeWalker;

template<class C, class EE>
class KDTreeWalker {
public:
  typedef KDTreeNode<C, EE> Node;
  typedef KDTreeLeaf<C, EE> Leaf;
  typedef KDTreeInterior<C, EE> Interior;
  typedef KDTreeWalker<C, EE> Walker;

  KDTreeWalker(Node* root = 0): _node(root), _prev(0) {}

  static Walker makeAtFirstNonEmptyLeaf(Node* root) {
    Walker walker(root);
    walker.toNonEmptyLeaf();
    return walker;
  }

  static Walker makeAtEnd(Node* root) {
    return Walker(0, root);
  }

  static Walker makeAt(Node* node, Node* prev) {
    ASSERT(node == 0 ||
      prev == node->getParent() ||
      prev == &node->asInterior().getEqualOrLess() ||
      prev == &node->asInterior().getStrictlyGreater());
    return Walker(node, prev);
  }

  bool atBegin() const {return _prev == 0;}
  bool atEnd() const {return _node == 0;}
  bool atRoot() const {ASSERT(!atEnd()); return _node->getParent() == 0;}

  Node* getNode() {return _node;}
  const Node* getNode() const {return _node;}

  bool atLeaf() const {return _node->isLeaf();}
  bool atInterior() const {return _node->isInterior();}
  Leaf& asLeaf() {
    ASSERT(atLeaf());
    return _node->asLeaf();
  }
  const Leaf& asLeaf() const {
    ASSERT(atLeaf());
    return _node->asLeaf();
  }
  Interior& asInterior() {
    ASSERT(atInterior());
    return _node->asInterior();
  }
  const Interior& asInterior() const {
    ASSERT(atInterior());
    return _node->asInterior();
  }

  Node* getParent() {
    ASSERT(!atEnd());
    return _node->getParent();
  }
  void toParent() {
    ASSERT(!atEnd());
    _prev = _node;
    _node = _node->getParent();
  }

  Node& getEqualOrLess() {
    ASSERT(_node != 0);
    ASSERT(atInterior());
    return asInterior().getEqualOrLess();
  }
  const Node& getEqualOrLess() const {
    return const_cast<Walker&>(*this).getEqualOrLess();
  }
  void toEqualOrLess() {
    ASSERT(_node != 0);
    ASSERT(atInterior());
    _prev = _node;
    _node = &asInterior().getEqualOrLess();
    ASSERT(_node->getParent() == _prev);
  }

  Node& getStrictlyGreater() {
    ASSERT(_node != 0);
    ASSERT(atInterior());
    return asInterior().getStrictlyGreater();
  }
  const Node& getStrictlyGreater() const {
    return const_cast<Walker&>(*this).getStrictlyGreater();
  }
  void toStrictlyGreater() {
    ASSERT(_node != 0);
    ASSERT(atInterior());
    _prev = _node;
    _node = &asInterior().getStrictlyGreater();
    ASSERT(_node->getParent() == _prev);
  }

  /** Moves to next node in depth first traversal. */
  void next() {
    ASSERT(!atEnd());
    ASSERT(_prev == getParent() || !atLeaf());
    if (_prev == getParent()) {
      if (atLeaf())
        toParent();
      else
        toEqualOrLess();
    } else if (_prev == &getEqualOrLess())
      toStrictlyGreater();
    else {
      ASSERT(_prev == &getStrictlyGreater());
      toParent();
    }

    ASSERT(_node == 0 || _prev == getParent() || !atLeaf());
    ASSERT((_node == 0 && _prev->getParent() == 0) ||
      _prev == _node->getParent() ||
      _prev == &_node->asInterior().getEqualOrLess() ||
      _prev == &_node->asInterior().getStrictlyGreater());

  }
  /** Returns true if next() will move to the strictlyGreater child. */
  bool nextIsStrictlyGreater() const {
    ASSERT(!atEnd());
    return atInterior() && _prev == &getEqualOrLess();
  }

  /** Moves to next non empty leaf in depth first traversal. */
  void nextNonEmptyLeaf() {
    ASSERT(!atEnd());
    do {
      next();
    } while (!atEnd() && (!atLeaf() || asLeaf().empty()));
  }

  /** Moves to next non empty leaf in depth first traversal or
   stays in same position if already at a leaf. */
  void toNonEmptyLeaf() {
    while (!atEnd() && (!atLeaf() || asLeaf().empty()))
      next();
  }

  /** Moves to previous node in depth first traversal. */
  void prev() {
    ASSERT(!atBegin());
    if (atEnd()) {
      _node = _prev;
      _prev = atLeaf() ? 0 : &getStrictlyGreater();
    } else if (_prev == getParent()) {
      if (&_prev->asInterior().getEqualOrLess() == _node) {
        _node = _prev;
        _prev = getParent();
      } else {
        ASSERT(&_prev->asInterior().getStrictlyGreater() == _node);
        _node = _prev;
        _prev = &getEqualOrLess();
      }
    } else if (_prev->isLeaf()) {
      std::swap(_node, _prev);
    } else {
      ASSERT(_prev == &getEqualOrLess() ||
        _prev == &getStrictlyGreater());
      _node = _prev;
      _prev = &getStrictlyGreater();
    }
  }

  /** Moves to previous non empty leaf in depth first traversal. */
  void toPrevNonEmptyLeaf() {
    ASSERT(!atBegin());
    do {
      prev();
    } while (!atBegin() && (!atLeaf() || asLeaf().empty()));
  }

  bool operator==(const Walker& walker) const {
    return _node == walker._node && _prev == walker._prev;
  }
  Walker& operator=(const Walker& walker) {
    _node = walker._node;
    _prev = walker._prev;
    return *this;
  }

private:
  KDTreeWalker(Node* node, Node* previous): _node(node), _prev(previous) {}

  Node* _node;
  Node* _prev;
};

#endif
