#ifndef MATHIC_BINARY_K_D_TREE_GUARD
#define MATHIC_BINARY_K_D_TREE_GUARD

#include "stdinc.h"
#include "DivMask.h"
#include "KDEntryArray.h"
#include "memtailor/memtailor.h"

namespace mathic {
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

  template<class C, class EE>
  class KDTreeNode {
    typedef KDTreeNode<C, EE> Node;
    typedef KDTreeLeaf<C, EE> Leaf;
    typedef KDTreeInterior<C, EE> Interior;
    friend class KDTreeInterior<C, EE>;
    friend class KDTreeLeaf<C, EE>;
  public:
    bool isLeaf() const {return _isLeaf;}
    const Leaf& asLeaf() const {
      MATHIC_ASSERT(isLeaf());
      return static_cast<const Leaf&>(*this);
    }
    Leaf& asLeaf() {
      MATHIC_ASSERT(isLeaf());
      return static_cast<Leaf&>(*this);
    }

    bool isInterior() const {return !isLeaf();}
    const Interior& asInterior() const {
      MATHIC_ASSERT(isInterior());
      return static_cast<Interior&>(*this);
    }
    Interior& asInterior() {
      MATHIC_ASSERT(isInterior());
      return static_cast<Interior&>(*this);
    }

    /** Partitions [begin, end) into two parts. The
        returned node has the information about the split, while the returned
        iterator it is such that the equal-or-less part of the partition
        is [begin, it) and the strictly-greater part is [it, end). */
    template<class Iter>
      static std::pair<Interior*, Iter> preSplit
      (size_t var, Iter begin, Iter end, memt::Arena& arena, const C& conf);

  protected:
  KDTreeNode(bool leaf): _isLeaf(leaf) {}

  private:
    class SplitEqualOrLess;
    const bool _isLeaf;
  };

  template<class C, class EE>
    class KDTreeInterior : public KDTreeNode<C, EE>, public mathic::DivMask::HasDivMask<C::UseTreeDivMask> {
  public:
    typedef typename C::Exponent Exponent;
    typedef KDTreeInterior<C, EE> Interior;
    typedef KDTreeLeaf<C, EE> Leaf;
    typedef KDTreeNode<C, EE> Node;

    KDTreeInterior
      (Node& equalOrLess,
       Node& strictlyGreater,
       size_t var,
       Exponent exponent):
    Node(false),
      _equalOrLess(&equalOrLess),
      _strictlyGreater(&strictlyGreater),
      _var(var),
      _exponent(exponent) {
      }
    KDTreeInterior
      (size_t var,
       Exponent exponent):
    Node(false),
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

    using mathic::DivMask::HasDivMask<C::UseTreeDivMask>::updateToLowerBound;
    void updateToLowerBound(Node& node) {
      if (!C::UseTreeDivMask)
        return;
      if (node.isLeaf())
        mathic::DivMask::HasDivMask<C::UseTreeDivMask>::updateToLowerBound(node.asLeaf().entries());
      else
        mathic::DivMask::HasDivMask<C::UseTreeDivMask>::updateToLowerBound(node.asInterior());
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

    typedef DivMask::Calculator<C> DivMaskCalculator;

    KDTreeLeaf(memt::Arena& arena, const C& conf);

    /** Copies [begin, end) into the new leaf. */
    template<class Iter>
    KDTreeLeaf(Iter begin, Iter end, memt::Arena& arena,
      const DivMaskCalculator& calc, const C& conf);

    EntryArray<C, EE>& entries() {return _entries;}
    const EntryArray<C, EE>& entries() const {return _entries;}

    Interior& split(Interior* parent, memt::Arena& arena, const C& conf);

  private:
    KDTreeLeaf(const KDTreeLeaf& t); // unavailable
    void operator=(const KDTreeLeaf&); // unavailable

    EntryArray<C, EE> _entries;
  };

  template<class C>
  class BinaryKDTree {
  public:
    typedef typename C::Monomial Monomial;
    typedef typename C::Entry Entry;
    typedef typename C::Exponent Exponent;
    typedef typename DivMask::Extender<Entry, C::UseDivMask> ExtEntry;
    typedef typename DivMask::Extender<const Monomial&,C::UseDivMask> ExtMonoRef;
    typedef typename DivMask::Calculator<C> DivMaskCalculator;
    typedef KDTreeNode<C, ExtEntry> Node;
    typedef KDTreeInterior<C, ExtEntry> Interior;
    typedef KDTreeLeaf<C, ExtEntry> Leaf;
    typedef typename Leaf::iterator LeafIt;
    typedef C Configuration;
    static const bool UseDivMask = C::UseDivMask;

  public:
    BinaryKDTree(const C& configuration);
    ~BinaryKDTree();

    template<class MultipleOutput>
    size_t removeMultiples(const ExtMonoRef& monomial, MultipleOutput& out);

    void insert(const ExtEntry& entry);

    template<class Iter>
    void insert(Iter begin, Iter end, const DivMaskCalculator& calc);

    Entry* findDivisor(const ExtMonoRef& monomial);

    template<class DivisorOutput>
    void findAllDivisors(const ExtMonoRef& monomial, DivisorOutput& out);

    template<class EntryOutput>
    void forAll(EntryOutput& out);

    void clear();

    size_t getMemoryUsage() const;

    C& getConfiguration() {return _conf;}

#ifdef MATHIC_DEBUG
    bool debugIsValid() const;
#endif

  private:
    BinaryKDTree(const BinaryKDTree<C>&); // unavailable
    void operator=(const BinaryKDTree<C>&); // unavailable
  
    template<class Iter>
    struct InsertTodo {
      Iter begin;
      Iter end;
      Interior* parent;
    };

    memt::Arena _arena; // Everything permanent allocated from here.
    C _conf; // User supplied configuration.
    mutable std::vector<Node*> _tmp; // For navigating the tree.
    Node* _root; // Root of the tree.
  };

  template<class C>
  BinaryKDTree<C>::BinaryKDTree(const C& configuration):
  _conf(configuration) {
    MATHIC_ASSERT(_conf.getLeafSize() >= 2);
    _root = new (_arena.allocObjectNoCon<Leaf>()) Leaf(_arena, _conf);
    MATHIC_ASSERT(debugIsValid());
  }

  template<class C>
  BinaryKDTree<C>::~BinaryKDTree() {
    clear();
  }

  template<class C>
  template<class MO>
  size_t BinaryKDTree<C>::removeMultiples(const ExtMonoRef& extMonomial, MO& out) {

    MATHIC_ASSERT(_tmp.empty());
    size_t removedCount = 0;
    Node* node = _root;
    while (true) {
      while (node->isInterior()) {
        Interior& interior = node->asInterior();
        if (!(interior.getExponent() <
              _conf.getExponent(extMonomial.get(), interior.getVar())))
          _tmp.push_back(&interior.getEqualOrLess());
        node = &interior.getStrictlyGreater();
      }
      MATHIC_ASSERT(node->isLeaf());
      removedCount += node->asLeaf().entries().removeMultiples(extMonomial, out, _conf);
      if (_tmp.empty())
        break;
      node = _tmp.back();
      _tmp.pop_back();
    }
    MATHIC_ASSERT(_tmp.empty());
    MATHIC_ASSERT(debugIsValid());
    return removedCount;
  }

  template<class C>
  void BinaryKDTree<C>::insert(const ExtEntry& extEntry) {
    Interior* parent = 0;
    Node* node = _root;
    while (node->isInterior()) {
      parent = &node->asInterior();
      parent->updateToLowerBound(extEntry);
      node = &parent->getChildFor(extEntry, _conf);
    }
    Leaf* leaf = &node->asLeaf();

    MATHIC_ASSERT(leaf->entries().size() <= _conf.getLeafSize());
    if (leaf->entries().size() == _conf.getLeafSize()) {
      Interior& interior = leaf->split(parent, _arena, _conf);
      interior.updateToLowerBound(extEntry);
      if (parent == 0) {
        ASSERT(leaf == _root);
        _root = &interior;
      }
      leaf = &interior.getChildFor(extEntry, _conf).asLeaf();
    }
    MATHIC_ASSERT(leaf->entries().size() < _conf.getLeafSize());
    leaf->entries().insert(extEntry, _conf);
    MATHIC_ASSERT(debugIsValid());
  }

  /// @todo: this function is too big and it knows too much about the details
  /// inside nodes. Also, it allocates a std::vector every time.
  template<class C>
  template<class Iter>
  void BinaryKDTree<C>::insert(Iter insertBegin, Iter insertEnd, const DivMaskCalculator& calc) {
      // @todo: rename to rebuild or something like that

    _arena.freeAllAllocs();
    _root = 0;

    typedef InsertTodo<Iter> Task;
    typedef std::vector<Task> TaskCont;
    TaskCont todo;

    _root = 0;
    Interior* parent = 0;
    bool isEqualOrLessChild = false;
    while (true) {
      Node* node = 0;
      const size_t insertCount = std::distance(insertBegin, insertEnd);
      const bool isLeaf = (insertCount <= _conf.getLeafSize());
      if (isLeaf)
        node = new (_arena.allocObjectNoCon<Leaf>())
          Leaf(insertBegin, insertEnd, _arena, calc, _conf);
      else {
        const size_t var =
          (parent == 0 ? static_cast<size_t>(-1) : parent->getVar());
        std::pair<Interior*, Iter> p =
          Node::preSplit(var, insertBegin, insertEnd, _arena, _conf);
        MATHIC_ASSERT(p.second != insertBegin && p.second != insertEnd);
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

      if (parent == 0) {
        ASSERT(_root == 0);
        _root = node;
      } else if (isEqualOrLessChild)
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
        // only strictlyGreater goes on todo
        isEqualOrLessChild = false;
      } else {
        isEqualOrLessChild = true;
        parent = &node->asInterior();
        // continue with equal-or-less as next item      
      }
    }

    if (C::UseTreeDivMask) {
      // Set div tree masks bottom up.
      typedef std::vector<Interior*> NodeCont;
      NodeCont nodes;
      if (_root->isInterior())
        nodes.push_back(&_root->asInterior());
      for (size_t i = 0; i < nodes.size(); ++i) {
        Interior* node = nodes[i];
        if (node->getEqualOrLess().isInterior())
          nodes.push_back(&node->getEqualOrLess().asInterior());
        if (node->getStrictlyGreater().isInterior())
          nodes.push_back(&node->getStrictlyGreater().asInterior());
      }
      typename NodeCont::reverse_iterator it = nodes.rbegin();
      typename NodeCont::reverse_iterator end = nodes.rend();
      for (; it != end; ++it) {
        Interior* node = *it;
        node->updateToLowerBound(node->getEqualOrLess());
        node->updateToLowerBound(node->getStrictlyGreater());
      }
    }
    MATHIC_ASSERT(debugIsValid());
  }

  template<class C>
  typename BinaryKDTree<C>::Entry* BinaryKDTree<C>::findDivisor
    (const ExtMonoRef& extMonomial) {

    MATHIC_ASSERT(debugIsValid());
    MATHIC_ASSERT(_tmp.empty());
    Node* node = _root;
    while (true) {
      while (node->isInterior()) {
        Interior& interior = node->asInterior();
        if (C::UseTreeDivMask &&
            !interior.getDivMask().canDivide(extMonomial.getDivMask()))
          goto next;

        if (interior.getExponent() <
            _conf.getExponent(extMonomial.get(), interior.getVar()))
          _tmp.push_back(&interior.getStrictlyGreater());
        node = &interior.getEqualOrLess();
      }

      {
        MATHIC_ASSERT(node->isLeaf());
        Leaf& leaf = node->asLeaf();
        LeafIt leafIt = leaf.entries().findDivisor(extMonomial, _conf);
        if (leafIt != leaf.entries().end()) {
          MATHIC_ASSERT(_conf.divides(leafIt->get(), extMonomial.get()));
          _tmp.clear();
          return &leafIt->get();
        }
      }
    next:
      if (_tmp.empty())
        break;
      node = _tmp.back();
      _tmp.pop_back();
    }
    MATHIC_ASSERT(_tmp.empty());
    return 0;
  }

  template<class C>
  template<class DO>
  void BinaryKDTree<C>::findAllDivisors(const ExtMonoRef& extMonomial, DO& output) {
    MATHIC_ASSERT(_tmp.empty());
    Node* node = _root;
    while (true) {
      while (node->isInterior()) {
        Interior& interior = node->asInterior();
        if (interior.getExponent() <
            _conf.getExponent(extMonomial.get(), interior.getVar()))
          _tmp.push_back(&interior.getStrictlyGreater());
        node = &interior.getEqualOrLess();
      }
      MATHIC_ASSERT(node->isLeaf());
      Leaf& leaf = node->asLeaf();
      if (!leaf.entries().findAllDivisors(extMonomial, output, _conf)) {
        _tmp.clear();
        break;
      }
      if (_tmp.empty())
        break;
      node = _tmp.back();
      _tmp.pop_back();
    }
    MATHIC_ASSERT(_tmp.empty());
  }

  template<class C>
  template<class EO>
  void BinaryKDTree<C>::forAll(EO& output) {
    MATHIC_ASSERT(_tmp.empty());
    Node* node = _root;
    while (true) {
      while (node->isInterior()) {
        Interior& interior = node->asInterior();
        _tmp.push_back(&interior.getStrictlyGreater());
        node = &interior.getEqualOrLess();
      }
      MATHIC_ASSERT(node->isLeaf());
      Leaf& leaf = node->asLeaf();
      if (!leaf.entries().forAll(output)) {
        _tmp.clear();
        break;
      }
      if (_tmp.empty())
        break;
      node = _tmp.back();
      _tmp.pop_back();
    }
    MATHIC_ASSERT(_tmp.empty());
  }

  template<class C>
  void BinaryKDTree<C>::clear() {
    MATHIC_ASSERT(_tmp.empty());
    if (_root != 0)
      _tmp.push_back(_root);
    while (!_tmp.empty()) {
      Node* node = _tmp.back();
      _tmp.pop_back();
      while (node->isInterior()) {
        _tmp.push_back(&node->asInterior().getStrictlyGreater());
        node = &node->asInterior().getEqualOrLess();
      }
      MATHIC_ASSERT(node->isLeaf());
      node->asLeaf().entries().clear();
    }
    _arena.freeAllAllocs();
    _root = new (_arena.allocObjectNoCon<Leaf>()) Leaf(_arena, _conf);
  }

  template<class C>
  size_t BinaryKDTree<C>::getMemoryUsage() const {
	size_t sum = _arena.getMemoryUsage();
	sum += _tmp.capacity() * sizeof(_tmp.front());
	return sum;
  }

#ifdef MATHIC_DEBUG
  template<class C>
  bool BinaryKDTree<C>::debugIsValid() const {
    MATHIC_ASSERT(_tmp.empty());
    MATHIC_ASSERT(!_conf.getDoAutomaticRebuilds() || _conf.getRebuildRatio() > 0);

    MATHIC_ASSERT(_root != 0);
    if (_root->isLeaf() && _root->asLeaf().entries().empty())
      return true;

    std::vector<Node*> nodes;
    nodes.push_back(_root);
    size_t sizeSum = 0;
    for (size_t i = 0; i < nodes.size(); ++i) {
      Node* node = nodes[i];
      if (node->isInterior()) {
        MATHIC_ASSERT(node->asInterior().getVar() < _conf.getVarCount());
        nodes.push_back(&node->asInterior().getStrictlyGreater());
        nodes.push_back(&node->asInterior().getEqualOrLess());
      } else
        sizeSum += node->asLeaf().entries().size();
    }

    MATHIC_ASSERT(_tmp.empty());
    for (size_t i = 0; i < nodes.size(); ++i) {
      Node* nodei = nodes[i];
      if (nodei->isLeaf()) {
        Leaf& leaf = nodei->asLeaf();
        typedef typename Leaf::const_iterator LeafCIter;
        if (C::UseTreeDivMask) {
          for (LeafCIter it = leaf.entries().begin(); it != leaf.entries().end(); ++it) {
            MATHIC_ASSERT(leaf.entries().getDivMask().canDivide(it->getDivMask()));
          }
        }
        continue;
      }
      Interior& interior = nodei->asInterior();

      ASSERT(_tmp.empty());
      _tmp.push_back(&interior.getEqualOrLess());
      while (!_tmp.empty()) {
        Node* node = _tmp.back();
        _tmp.pop_back();
        if (C::UseTreeDivMask) {
          if (node->isInterior())
            MATHIC_ASSERT(interior.getDivMask().canDivide(node->asInterior().getDivMask()));
          else
            MATHIC_ASSERT(interior.getDivMask().canDivide(node->asLeaf().entries().getDivMask()));
        }
        if (node->isInterior()) {
          _tmp.push_back(&node->asInterior().getStrictlyGreater());
          _tmp.push_back(&node->asInterior().getEqualOrLess());
        } else {
          Leaf& leaf = node->asLeaf();
          typename Leaf::const_iterator stop = leaf.entries().end();
          for (typename Leaf::const_iterator it = leaf.entries().begin(); it != stop; ++it) {
            MATHIC_ASSERT(!(interior.getExponent() <
              _conf.getExponent(it->get(), interior.getVar())));
          }
        }
      }
    }
    return true;
  }
#endif

  template<class C, class EE>
  KDTreeLeaf<C, EE>::KDTreeLeaf(memt::Arena& arena, const C& conf):
  Node(true), _entries(arena, conf) {}

  template<class C, class EE>
  template<class Iter>
  KDTreeLeaf<C, EE>::KDTreeLeaf(
    Iter begin,
    Iter end,
    memt::Arena& arena,
    const DivMaskCalculator& calc, const C& conf):
    Node(true), _entries(begin, end, arena, calc, conf) {}

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
    (size_t var,
     Iter begin,
     Iter end,
     memt::Arena& arena,
     const C& conf) {
    MATHIC_ASSERT(begin != end);
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
        Interior(var, exp);
      SplitEqualOrLess cmp(var, exp, conf);
      Iter middle = std::partition(begin, end, cmp);
      return std::make_pair(interior, middle);
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
  KDTreeInterior<C, EE>&
  KDTreeLeaf<C, EE>::split(Interior* parent, memt::Arena& arena, const C& conf) {
    MATHIC_ASSERT(conf.getVarCount() > 0);
    MATHIC_ASSERT(entries().size() >= 2);
    // MATHIC_ASSERT not all equal
    Leaf& other = *new (arena.allocObjectNoCon<Leaf>()) Leaf(arena, conf);
    size_t var = (parent == 0 ? static_cast<size_t>(-1) : parent->getVar());
    typename C::Exponent exp;
    while (true) {
      var = (var + 1) % conf.getVarCount();

      if (1) {
        typename C::Exponent min = conf.getExponent(entries().front().get(), var);
        typename C::Exponent max = conf.getExponent(entries().front().get(), var);
        for (iterator it = entries().begin(); it != entries().end(); ++it) {
          min = std::min(min, conf.getExponent(it->get(), var));
          max = std::max(max, conf.getExponent(it->get(), var));
        }
        if (min == max && entries().size() > 1) {
          // todo: avoid infinite loop if all equal
          continue;
        }
        exp = min + (max - min) / 2; // this formula for avg avoids overflow

        iterator newEnd = entries().begin();
        for (iterator it = entries().begin(); it != entries().end(); ++it) {
          if (exp < conf.getExponent(it->get(), var))
            other.entries().push_back(*it);
          else {
            if (it != newEnd) {
              *newEnd = *it;
            }
            ++newEnd;
          }
        }
        while (newEnd != entries().end())
          entries().pop_back();
      } else {
        iterator middle = entries().begin() + entries().size() / 2;
        ExpOrder<C, EE> order(var, conf);

        std::nth_element(entries().begin(), middle, entries().end(), order);
        if (middle != entries().end()) {
          exp = conf.getExponent(middle->get(), var);
          while (middle != entries().end() && conf.getExponent(middle->get(), var) == exp)
            ++middle;
        }
        if (middle == entries().end() && entries().size() > 1) {
          // todo: avoid infinite loop if all equal.
          continue; // bad split, use another variable
        }
        MATHIC_ASSERT(middle != entries().end());
        MATHIC_ASSERT(exp != conf.getExponent(middle->get(), var));

#ifdef MATHIC_DEBUG
        for (iterator it = entries().begin(); it != middle; ++it) {
          MATHIC_ASSERT(!(exp < conf.getExponent(it->get(), var)));
        }
        for (iterator it = middle; it != entries().end(); ++it) {
          MATHIC_ASSERT(!(conf.getExponent(it->get(), var) < exp));
        }
#endif
        // nth_element does not guarantee where equal elements go,
        // so we cannot just copy [middle, end()).
        iterator newEnd = entries().begin();
        for (iterator it = entries().begin(); it != entries().end(); ++it) {
          if (exp < conf.getExponent(it->get(), var))
            other.entries().push_back(*it);
          else {
            if (it != newEnd) {
              *newEnd = *it;
            }
            ++newEnd;
          }
        }
        while (newEnd != entries().end())
          entries().pop_back();

      }
      MATHIC_ASSERT(other.entries().size() < conf.getLeafSize());
      MATHIC_ASSERT(entries().size() < conf.getLeafSize());
      break;
    }

    Interior& interior = *new (arena.allocObjectNoCon<Interior>())
      Interior(*this, other, var, exp);
    if (C::UseTreeDivMask) {
      entries().resetDivMask();
      for (const_iterator it = entries().begin(); it != entries().end(); ++it)
        entries().updateToLowerBound(*it);
      interior.updateToLowerBound(entries());
      interior.updateToLowerBound(other.entries());
    }
    if (parent != 0) {
      if (&parent->getEqualOrLess() == this)
        parent->setEqualOrLess(&interior);
      else {
        MATHIC_ASSERT(&parent->getStrictlyGreater() == this);
        parent->setStrictlyGreater(&interior);
      }
    }

    return interior;
  }
}

#endif
