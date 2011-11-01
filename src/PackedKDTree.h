#ifndef MATHIC_PACKED_K_D_TREE_GUARD
#define MATHIC_PACKED_K_D_TREE_GUARD

#include "stdinc.h"
#include "DivMask.h"
#include "KDEntryArray.h"
#include "memtailor/memtailor.h"
#include <ostream>

namespace mathic {
  template<class C>
  class PackedKDTree {
  public:
    typedef typename C::Monomial Monomial;
    typedef typename C::Entry Entry;
    typedef typename C::Exponent Exponent;
    typedef typename DivMask::Extender<Entry, C::UseDivMask> ExtEntry;
    typedef typename DivMask::Extender<const Monomial&,C::UseDivMask> ExtMonoRef;
    typedef typename DivMask::Calculator<C> DivMaskCalculator;

    struct ExpOrder {
    ExpOrder(size_t var, const C& conf): _var(var), _conf(conf) {}
      bool operator()(const ExtEntry& a, const ExtEntry& b) const {
        return _conf.getExponent(a.get(), _var) < _conf.getExponent(b.get(), _var);
      }
    private:
      const size_t _var;
      const C& _conf;
    };

    private:
    typedef C Configuration;
    static const bool UseDivMask = C::UseDivMask;

    class Node {
    public:
      static Node* makeNode(memt::Arena& arena, const C& conf) {
        return new (arena.alloc(sizeOf(0)))
          Node(arena, conf);
      }

      template<class Iter>
      static Node* makeNode(Iter begin, Iter end, memt::Arena& arena,
        const C& conf, size_t childCount) {
        return new (arena.alloc(sizeOf(childCount)))
          Node(begin, end, arena, conf, childCount);
      }

      template<class Iter>
      static Node* makeNode(Iter begin, Iter end, memt::Arena& arena,
        const DivMaskCalculator& calc, const C& conf, size_t childCount) {
        return new (arena.alloc(sizeOf(childCount)))
          Node(begin, end, arena, calc, conf, childCount);
      }

      static size_t sizeOf(size_t childCount) {
        if (childCount > 0)
          --childCount; // array has size 1, so one element already there
        return sizeof(Node) + childCount * sizeof(Child);
      }

      struct Child {
        size_t var;
        Exponent exponent;
        Node* node;
      };
      typedef Child* iterator;
      typedef Child const* const_iterator;
      iterator childBegin() {return _childrenMemoryBegin;}
      const_iterator childBegin() const {return _childrenMemoryBegin;}
      iterator childEnd() {return _childrenEnd;}
      const_iterator childEnd() const {return _childrenEnd;}

      bool hasChildren() const {return childBegin() != childEnd();}
      template<class Ext>
      bool inChild(const_iterator child, const Ext& ext, const C& conf) const {
        return child->exponent < conf.getExponent(ext.get(), child->var);
      }

      mathic::KDEntryArray<C, ExtEntry>& entries() {return _entries;}
      const KDEntryArray<C, ExtEntry>& entries() const {return _entries;}

      Node* split(Child* childFromParent, memt::Arena& arena, const C& conf);

      /** Partitions [begin, end) into two parts. The
          returned node has the information about the split, while the returned
          iterator it is such that the equal-or-less part of the partition
          is [begin, it) and the strictly-greater part is [it, end). These two
          ranges should be used to construct nodes that are then set as the
          children of the returned node. */
/*      template<class Iter>
        static std::pair<Interior*, Iter> preSplit
        (size_t var, Iter begin, Iter end, memt::Arena& arena, const C& conf);*/

    private:
      Node(const Node&); // unavailable
      void operator=(const Node&); // unavailable

      Node(memt::Arena& arena, const C& conf);
      
      template<class Iter>
      Node(Iter begin, Iter end, memt::Arena& arena,
        const C& conf, size_t childCount);

      template<class Iter>
      Node(Iter begin, Iter end, memt::Arena& arena,
        const DivMaskCalculator& calc, const C& conf, size_t childCount);

      class SplitEqualOrLess;

      KDEntryArray<C, ExtEntry> _entries;
      // Array has size 1 to appease compiler since size 0 produces warnings
      // or errors. Actual size can be greater if more memory has been
      // allocated for the node than sizeof(Node).
      Child* _childrenEnd; // points into _childrenMemoryBegin
      Child _childrenMemoryBegin[1];
    };

  public:
    PackedKDTree(const C& configuration);
    ~PackedKDTree();

    template<class MultipleOutput>
    size_t removeMultiples(const ExtMonoRef& monomial, MultipleOutput& out);

    void insert(const ExtEntry& entry);

    template<class Iter>
    void reset(Iter begin, Iter end, const DivMaskCalculator& calc);

    Entry* findDivisor(const ExtMonoRef& monomial);

    template<class DivisorOutput>
    void findAllDivisors(const ExtMonoRef& monomial, DivisorOutput& out);

    template<class EntryOutput>
    void forAll(EntryOutput& out);

    void clear();

    size_t getMemoryUsage() const;

    void print(std::ostream& out) const;

    C& getConfiguration() {return _conf;}

#ifdef MATHIC_DEBUG
    bool debugIsValid() const;
#endif

  private:
    PackedKDTree(const PackedKDTree<C>&); // unavailable
    void operator=(const PackedKDTree<C>&); // unavailable
  
    template<class Iter>
    struct InsertTodo {
      Iter begin;
      Iter end;
      Exponent exp;
      size_t var;
      typename Node::Child* fromParent;
    };

    memt::Arena _arena; // Everything permanent allocated from here.
    C _conf; // User supplied configuration.
    mutable std::vector<Node*> _tmp; // For navigating the tree.
    Node* _root; // Root of the tree. Cannot be null.
  };

  template<class C>
  PackedKDTree<C>::PackedKDTree(const C& configuration):
  _conf(configuration) {
    MATHIC_ASSERT(_conf.getLeafSize() >= 2);
    _root = Node::makeNode(_arena, _conf);
    MATHIC_ASSERT(debugIsValid());
  }

  template<class C>
  PackedKDTree<C>::~PackedKDTree() {
    clear();
  }

  template<class C>
  PackedKDTree<C>::Node::Node(memt::Arena& arena, const C& conf):
  _entries(arena, conf) {
    _childrenEnd = childBegin();
  }

  template<class C>
  template<class Iter>
  PackedKDTree<C>::Node::Node(
    Iter begin,
    Iter end,
    memt::Arena& arena,
    const C& conf,
    size_t childCount):
  _entries(begin, end, arena, conf) {
    _childrenEnd = childBegin() + childCount;
  }

  template<class C>
  template<class Iter>
  PackedKDTree<C>::Node::Node(
    Iter begin,
    Iter end,
    memt::Arena& arena,
    const DivMaskCalculator& calc,
    const C& conf,
    size_t childCount):
    _entries(begin, end, arena, calc, conf) {
    _childrenEnd = childBegin() + childCount;
  }

  template<class C>
  template<class MO>
  size_t PackedKDTree<C>::removeMultiples(
    const ExtMonoRef& extMonomial,
    MO& out
  ) {
    MATHIC_ASSERT(_tmp.empty());
    size_t removedCount = 0;
    Node* node = _root;
    while (true) {
      for (typename Node::const_iterator it = node->childBegin();
        it != node->childEnd(); ++it) {
        _tmp.push_back(it->node);
        if (node->inChild(it, extMonomial, _conf))
          goto stopped;
      }
      removedCount += node->entries().removeMultiples(extMonomial, out, _conf);
stopped:;
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
  void PackedKDTree<C>::insert(const ExtEntry& extEntry) {
    // find node in which to insert extEntry
    typename Node::Child* parentChild = 0;
    Node* node = _root;
    typename Node::iterator child = node->childBegin();
    while (true) {
      if (child == node->childEnd()) {
        MATHIC_ASSERT(node->entries().size() <= _conf.getLeafSize());
        if (node->entries().size() < _conf.getLeafSize())
          break;
        // split node as it is full
        ASSERT(node == _root || parentChild != 0);
        const size_t childOffset = child - node->childBegin();
        node = node->split(parentChild, _arena, _conf);
        child = node->childBegin() + childOffset;
        if (parentChild == 0)
          _root = node;
      } else if (node->inChild(child, extEntry, _conf)) {
        parentChild = &*child;
        node = child->node;
        child = node->childBegin();
      } else
        ++child;
    }

    // insert into node
    MATHIC_ASSERT(child == node->childEnd());
    MATHIC_ASSERT(node->entries().size() < _conf.getLeafSize());
    node->entries().insert(extEntry, _conf);
    MATHIC_ASSERT(debugIsValid());

    //print(std::cerr); std::cerr << std::flush; // todo: debug code, remove
  }

  template<class C>
  template<class Iter>
  void PackedKDTree<C>::reset(Iter insertBegin, Iter insertEnd, const DivMaskCalculator& calc) {
    clear();
    _arena.freeAllAllocs();
    _root = 0; // temporary value

    typedef InsertTodo<Iter> Task;
    typedef std::vector<Task> TaskCont;
    TaskCont todo;
    TaskCont children;

    {
      Task initialTask;
      initialTask.begin = insertBegin;
      initialTask.end = insertEnd;
      initialTask.var = static_cast<size_t>(-1);
      initialTask.fromParent = 0;
      todo.push_back(initialTask);
    }
    while (!todo.empty()) {
      Iter begin = todo.back().begin;
      Iter end = todo.back().end;
      size_t var = todo.back().var;
      typename Node::Child* fromParent = todo.back().fromParent;
      if (fromParent != 0) {
        fromParent->var = var;
        fromParent->exponent = todo.back().exp;
      }
      todo.pop_back();

      // split off children until reaching few enough entries
      while (_conf.getLeafSize() <
        static_cast<size_t>(std::distance(begin, end))) {
        Task child;
        Iter middle = KDEntryArray<C, ExtEntry>::
          split(begin, end, var, child.exp, _conf);
        MATHIC_ASSERT(begin < middle && middle < end);
        MATHIC_ASSERT(var < _conf.getVarCount());
        child.begin = middle;
        child.end = end;
        child.var = var;
        children.push_back(child);
        // set up equal-or-less
        end = middle;
      }
      Node* node = Node::makeNode
        (begin, end, _arena, calc, _conf, children.size());
      if (_root == 0)
        _root = node;
      if (fromParent != 0)
        fromParent->node = node;
      for (size_t child = 0; child < children.size(); ++child) {
        children[child].fromParent = &*(node->childBegin() + child);
        todo.push_back(children[child]);
      }
      children.clear();
    }
    MATHIC_ASSERT(_root != 0);
    MATHIC_ASSERT(debugIsValid());
  }

  template<class C>
  typename PackedKDTree<C>::Entry* PackedKDTree<C>::findDivisor
    (const ExtMonoRef& extMonomial) {
    MATHIC_ASSERT(_tmp.empty());
    Node* node = _root;
    while (true) {
      // look for divisor in entries of node
      typename KDEntryArray<C, ExtEntry>::iterator it =
        node->entries().findDivisor(extMonomial, _conf);
      if (it != node->entries().end()) {
        MATHIC_ASSERT(_conf.divides(it->get(), extMonomial.get()));
        _tmp.clear();
        return &it->get();
      }

      // remember relevant children
      for (typename Node::const_iterator it = node->childBegin();
        it != node->childEnd(); ++it)
        if (node->inChild(it, extMonomial, _conf))
          _tmp.push_back(it->node);

      // grab next node to process
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
  void PackedKDTree<C>::findAllDivisors(const ExtMonoRef& extMonomial, DO& output) {
    MATHIC_ASSERT(_tmp.empty());
    Node* node = _root;
    while (true) {
      if (!node->entries().findAllDivisors(extMonomial, output, _conf)) {
        _tmp.clear();
        break;
      }
      for (typename Node::const_iterator it = node->childBegin();
        it != node->childEnd(); ++it)
        if (node->inChild(it, extMonomial, _conf))
          _tmp.push_back(it->node);
      if (_tmp.empty())
        break;
      node = _tmp.back();
      _tmp.pop_back();
    }
    MATHIC_ASSERT(_tmp.empty());
  }

  template<class C>
  template<class EO>
  void PackedKDTree<C>::forAll(EO& output) {
    MATHIC_ASSERT(_tmp.empty());
    Node* node = _root;
    while (true) {
      if (!node->entries().forAll(output)) {
        _tmp.clear();
        break;
      }
      for (typename Node::iterator it = node->childBegin();
        it != node->childEnd(); ++it)
        _tmp.push_back(it->node);
      if (_tmp.empty())
        break;
      node = _tmp.back();
      _tmp.pop_back();
    }
    MATHIC_ASSERT(_tmp.empty());
  }

  template<class C>
  void PackedKDTree<C>::clear() {
    MATHIC_ASSERT(_tmp.empty());
    // Call Entry destructors
    _tmp.push_back(_root);
    while (!_tmp.empty()) {
      Node* node = _tmp.back();
      _tmp.pop_back();
      node->entries().clear(); // calls destructors
      for (typename Node::iterator it = node->childBegin();
        it != node->childEnd(); ++it)
        _tmp.push_back(it->node);
    }
    _arena.freeAllAllocs();
    _root = Node::makeNode(_arena, _conf);
  }

  template<class C>
  size_t PackedKDTree<C>::getMemoryUsage() const {
    // todo: not accurate
	size_t sum = _arena.getMemoryUsage();
	sum += _tmp.capacity() * sizeof(_tmp.front());
	return sum;
  }

  template<class C>
  void PackedKDTree<C>::print(std::ostream& out) const {
    out << "<<<<<<<< PackedKDTree >>>>>>>>\n";
    MATHIC_ASSERT(_tmp.empty());
    Node* node = _root;
    while (true) {
      out << "**** Node " << node << "\nchildren:\n";
      for (typename Node::iterator it = node->childBegin();
        it != node->childEnd(); ++it) {
        _tmp.push_back(it->node);
        out << "Child " << ((it - node->childBegin()) + 1) << ": "
          << '>' << (it->var + 1) << '^' << it->exponent
          << ' ' << it->node << '\n';
      }
      for (size_t i = 0; i < node->entries().size(); ++i) {
        out << "Entry " << (i + 1) << ": "
          << node->entries().begin()[i].get() << '\n';
      }
      out << '\n';
      if (_tmp.empty())
        break;
      node = _tmp.back();
      _tmp.pop_back();
    }
    MATHIC_ASSERT(_tmp.empty());
  }

#ifdef MATHIC_DEBUG
  template<class C>
  bool PackedKDTree<C>::debugIsValid() const {
    //print(std::cerr); std::cerr << std::flush; // todo: debug remove
    MATHIC_ASSERT(_tmp.empty());
    MATHIC_ASSERT(!_conf.getDoAutomaticRebuilds() || _conf.getRebuildRatio() > 0);

    MATHIC_ASSERT(_root != 0);

    std::vector<Node*> nodes;
    nodes.push_back(_root);
    for (size_t i = 0; i < nodes.size(); ++i) {
      Node* node = nodes[i];
      for (typename Node::iterator it = node->childBegin();
        it != node->childEnd(); ++it) {
        MATHIC_ASSERT(it->var < _conf.getVarCount());
        nodes.push_back(it->node);
      }
    }

    MATHIC_ASSERT(_tmp.empty());
    for (size_t i = 0; i < nodes.size(); ++i) {
      Node* ancestor = nodes[i];
      for (typename Node::iterator ancestorIt = ancestor->childBegin();
        ancestorIt != ancestor->childEnd(); ++ancestorIt) {
        MATHIC_ASSERT(_tmp.empty());
        size_t var = ancestorIt->var;
        Exponent exp = ancestorIt->exponent;
        // check strictly greater than subtree
        _tmp.push_back(ancestorIt->node);
        while (!_tmp.empty()) {
          Node* node = _tmp.back();
          _tmp.pop_back();
          MATHIC_ASSERT(node->entries().
            allStrictlyGreaterThan(var, exp, _conf));
          for (typename Node::iterator it = node->childBegin();
            it != node->childEnd(); ++it)
            _tmp.push_back(it->node);
        }
        // check less than or equal to sub tree.
        MATHIC_ASSERT(ancestor->entries().allLessThanOrEqualTo(var, exp, _conf));
        typename Node::iterator restIt = ancestorIt;
        for (++restIt; restIt != ancestor->childEnd(); ++restIt)
          _tmp.push_back(restIt->node);
        while (!_tmp.empty()) {
          Node* node = _tmp.back();
          _tmp.pop_back();
          MATHIC_ASSERT(node->entries().
            allLessThanOrEqualTo(var, exp, _conf));
          for (typename Node::iterator it = node->childBegin();
            it != node->childEnd(); ++it)
            _tmp.push_back(it->node);
        }
      }
    }
    return true;
  }
#endif

  template<class C>
  class PackedKDTree<C>::Node::SplitEqualOrLess {
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
  typename PackedKDTree<C>::Node* PackedKDTree<C>::Node::split(
    Child* childFromParent,
    memt::Arena& arena,
    const C& conf
  ) {
    MATHIC_ASSERT(conf.getVarCount() > 0);
    MATHIC_ASSERT(entries().size() >= 2);
    size_t var;
    if (hasChildren())
      var = (childEnd() - 1)->var;
    else if (childFromParent != 0)
      var = childFromParent->var;
    else
      var = static_cast<size_t>(-1);
    typename C::Exponent exp;

    // there is not room to add another child, so make a
    // new Node with more space and put the Node we are splitting
    // off into the vacated space. It will fit as it has no
    // children at all.

    typename KDEntryArray<C, ExtEntry>::iterator middle =
      KDEntryArray<C, ExtEntry>::split
      (entries().begin(), entries().end(), var, exp, conf);
    ASSERT(middle != entries().begin());

    Node* copied = makeNode(entries().begin(), middle, arena, conf,
      std::distance(childBegin(), childEnd()) + 1);
    std::copy(childBegin(), childEnd(), copied->childBegin());
    Child newChild;
    newChild.var = var;
    newChild.exponent = exp;
    newChild.node = this;
    *(copied->childEnd() - 1) = newChild;
    if (childFromParent != 0)
      childFromParent->node = copied;
    // OK to call std::copy as begin is not in the range [middle, end).
    std::copy(middle, entries().end(), entries().begin());
    const size_t targetSize = std::distance(middle, entries().end());
    while (entries().size() > targetSize)
      entries().pop_back();
    _childrenEnd = childBegin();

    if (conf.getSortOnInsert())
      std::sort(entries().begin(), entries().end(), Comparer<C>(conf));
    return copied;
  }
}

#endif
