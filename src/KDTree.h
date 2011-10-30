#ifndef MATHIC_K_D_TREE_GUARD
#define MATHIC_K_D_TREE_GUARD

#include "stdinc.h"
//#include "KDTreeLeaf.h"
#include "DivMask.h"
#include "memtailor/memtailor.h"
#include <list>
#include <string>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <vector>

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
  class KDTreeNode : public mathic::DivMask::HasDivMask<C::UseTreeDivMask> {
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
        is [begin, it) and the striclty-greater part is [it, end). */
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
    class KDTreeInterior : public KDTreeNode<C, EE> {
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

    KDTreeLeaf(memt::Arena& arena, const C& conf);
    KDTreeLeaf(memt::Arena& arena, size_t capacity);

    void clear();

    typedef DivMask::Calculator<C> DivMaskCalculator;

    /** Copies [begin, end) into the new leaf. */
    template<class Iter>
      static Leaf* makeLeafCopy
      (Iter begin, Iter end, memt::Arena& arena,
       const DivMaskCalculator& calc, const C& conf);

    iterator begin() {return _begin;}
    const_iterator begin() const {return _begin;}
    iterator end() {return _end;}
    const_iterator end() const {return _end;}

    bool empty() const {return _begin == _end;}
    size_t size() const {return _end - _begin;}
    EE& front() {MATHIC_ASSERT(!empty()); return *_begin;}
    const EE& front() const {MATHIC_ASSERT(!empty()); return *_begin;}
    EE& back() {MATHIC_ASSERT(!empty()); return *(_end - 1);}
    const EE& back() const {MATHIC_ASSERT(!empty()); return *(_end - 1);}

    void push_back(const EE& entry);
    void pop_back();
    void insert(iterator it, const EE& entry);

    void insert(const EE& entry, const C& conf);

    /** Returns how many were removed. */
    template<class EM, class MO>
      size_t removeMultiples(const EM& monomial, MO& out, const C& conf);

    template<class EM>
    iterator findDivisor(const EM& extMonomial, const C& conf);

    /** Calls out.proceed(entry) for each entry that divides extMonomial.
        Stops and returns false if out.proceed(entry) returns false. Returns
        true if all calls out.proceed(entry) returned true. */
    template<class EM, class DO>
    bool findAllDivisors(const EM& extMonomial, DO& out, const C& conf);

    /** Calls out.proceed(entry) for each entry.
        Stops and returns false if out.proceed(entry) returns false. Returns
        true if all calls out.proceed(entry) returned true. */
    template<class EO>
    bool forAll(EO& eo);

    Interior& split(Interior* parent, memt::Arena& arena, const C& conf);

  private:
    KDTreeLeaf(const KDTreeLeaf& t); // unavailable
    void operator=(const KDTreeLeaf&); // unavailable

    iterator _begin;
    iterator _end;
#ifdef MATHIC_DEBUG
    const size_t _capacityDebug;
    const bool _sortOnInsertDebug;
#endif
  };

  template<class C>
  class KDTreeX {
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
    KDTreeX(const C& configuration);
    ~KDTreeX();

    template<class MultipleOutput>
    size_t removeMultiples(const Monomial& monomial, MultipleOutput& out);
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
    KDTreeX(const KDTreeX<C>&); // unavailable
    void operator=(const KDTreeX<C>&); // unavailable
  
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
  KDTreeX<C>::KDTreeX(const C& configuration):
  _conf(configuration) {
    MATHIC_ASSERT(_conf.getLeafSize() >= 2);
    _root = new (_arena.allocObjectNoCon<Leaf>()) Leaf(_arena, _conf);
    MATHIC_ASSERT(debugIsValid());
  }

  template<class C>
  KDTreeX<C>::~KDTreeX() {
    clear();
  }

  template<class C>
  template<class MO>
  size_t KDTreeX<C>::removeMultiples(const Monomial& monomial, MO& out) {
    ExtMonoRef extMonomial(monomial, _divMaskCalculator, _conf);
    return removeMultiples(extMonomial, out);
  }

  template<class C>
  template<class MO>
  size_t KDTreeX<C>::removeMultiples(const ExtMonoRef& extMonomial, MO& out) {

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
      removedCount += node->asLeaf().removeMultiples(extMonomial, out, _conf);
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
  void KDTreeX<C>::insert(const ExtEntry& extEntry) {
    Interior* parent = 0;
    Node* node = _root;
    while (node->isInterior()) {
      node->updateToLowerBound(extEntry);
      parent = &node->asInterior();
      node = &parent->getChildFor(extEntry, _conf);
    }
    Leaf* leaf = &node->asLeaf();

    MATHIC_ASSERT(leaf->size() <= _conf.getLeafSize());
    if (leaf->size() == _conf.getLeafSize()) {
      Interior& interior = leaf->split(parent, _arena, _conf);
      interior.updateToLowerBound(extEntry);
      if (parent == 0) {
        ASSERT(leaf == _root);
        _root = &interior;
      }
      leaf = &interior.getChildFor(extEntry, _conf).asLeaf();
    }
    MATHIC_ASSERT(leaf->size() < _conf.getLeafSize());
    leaf->insert(extEntry, _conf);
    MATHIC_ASSERT(debugIsValid());
  }

  /// @todo: this function is too big and it knows too much about the details
  /// inside nodes. Also, it allocates a std::vector every time.
  template<class C>
  template<class Iter>
  void KDTreeX<C>::insert(Iter insertBegin, Iter insertEnd, const DivMaskCalculator& calc) {
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
        node = Leaf::makeLeafCopy
          (insertBegin, insertEnd, _arena, calc, _conf);
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
  typename KDTreeX<C>::Entry* KDTreeX<C>::findDivisor
    (const ExtMonoRef& extMonomial) {

    MATHIC_ASSERT(debugIsValid());
    MATHIC_ASSERT(_tmp.empty());
    Node* node = _root;
    while (true) {
      while (node->isInterior()) {
        if (C::UseTreeDivMask &&
            !node->getDivMask().canDivide(extMonomial.getDivMask()))
          goto next;

        Interior& interior = node->asInterior();
        if (interior.getExponent() <
            _conf.getExponent(extMonomial.get(), interior.getVar()))
          _tmp.push_back(&interior.getStrictlyGreater());
        node = &interior.getEqualOrLess();
      }

      {
        MATHIC_ASSERT(node->isLeaf());
        Leaf& leaf = node->asLeaf();
        LeafIt leafIt = leaf.findDivisor(extMonomial, _conf);
        if (leafIt != leaf.end()) {
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
  void KDTreeX<C>::findAllDivisors(const ExtMonoRef& extMonomial, DO& output) {
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
      if (!leaf.findAllDivisors(extMonomial, output, _conf)) {
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
  void KDTreeX<C>::forAll(EO& output) {
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
      if (!leaf.forAll(output)) {
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
  void KDTreeX<C>::clear() {
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
      node->asLeaf().clear();
    }
    _arena.freeAllAllocs();
    _root = new (_arena.allocObjectNoCon<Leaf>()) Leaf(_arena, _conf);
  }

  template<class C>
  size_t KDTreeX<C>::getMemoryUsage() const {
	size_t sum = _arena.getMemoryUsage();
	sum += _tmp.capacity() * sizeof(_tmp.front());
	return sum;
  }

#ifdef MATHIC_DEBUG
  template<class C>
  bool KDTreeX<C>::debugIsValid() const {
    MATHIC_ASSERT(_tmp.empty());
    MATHIC_ASSERT(!_conf.getDoAutomaticRebuilds() || _conf.getRebuildRatio() > 0);

    MATHIC_ASSERT(_root != 0);
    if (_root->isLeaf() && _root->asLeaf().empty())
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
        sizeSum += node->asLeaf().size();
    }

    MATHIC_ASSERT(_tmp.empty());
    for (size_t i = 0; i < nodes.size(); ++i) {
      Node* nodei = nodes[i];
      if (nodei->isLeaf()) {
        Leaf& leaf = nodei->asLeaf();
        typedef typename Leaf::const_iterator LeafCIter;
        if (C::UseTreeDivMask) {
          for (LeafCIter it = leaf.begin(); it != leaf.end(); ++it) {
            MATHIC_ASSERT(leaf.getDivMask().canDivide(it->getDivMask()));
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
          MATHIC_ASSERT(interior.getDivMask().canDivide(node->getDivMask()));
        }
        if (node->isInterior()) {
          _tmp.push_back(&node->asInterior().getStrictlyGreater());
          _tmp.push_back(&node->asInterior().getEqualOrLess());
        } else {
          Leaf& leaf = node->asLeaf();
          typename Leaf::const_iterator stop = leaf.end();
          for (typename Leaf::const_iterator it = leaf.begin(); it != stop; ++it) {
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
  Node(true)
#ifdef MATHIC_DEBUG
    ,_capacityDebug(conf.getLeafSize())
    ,_sortOnInsertDebug(conf.getSortOnInsert())
#endif
    {
      _begin = arena.allocArrayNoCon<EE>(conf.getLeafSize()).first;
      _end = _begin;
    }

  template<class C, class EE>
    KDTreeLeaf<C,EE>::KDTreeLeaf(memt::Arena& arena, size_t capacity):
  Node(true, 0)
#ifdef MATHIC_DEBUG
    ,_capacityDebug(capacity)
    ,_sortOnInsertDebug(false)
#endif
    {
      _begin = arena.allocArrayNoCon<EE>(capacity).first;
      _end = _begin;
    }

  template<class C, class EE>
    void KDTreeLeaf<C, EE>::push_back(const EE& entry) {
    MATHIC_ASSERT(size() < _capacityDebug);
    new (_end) EE(entry);
    updateToLowerBound(entry);
    ++_end;
  }

  template<class C, class EE>
    void KDTreeLeaf<C, EE>::pop_back() {
    MATHIC_ASSERT(!empty());
    --_end;
    _end->~EE();
  }

  template<class C, class EE>
    void KDTreeLeaf<C, EE>::insert(iterator it, const EE& entry) {
    MATHIC_ASSERT(size() < _capacityDebug);
    if (it == end()) {
      push_back(entry);
      return;
    }
    push_back(back());
    iterator moveTo = end();
    for (--moveTo; moveTo != it; --moveTo)
      *moveTo = *(moveTo - 1);
    updateToLowerBound(entry);
    *it = entry;
  }

  template<class C, class EE>
    void KDTreeLeaf<C, EE>::insert(const EE& entry, const C& conf) {
    MATHIC_ASSERT(size() < _capacityDebug);
    if (!conf.getSortOnInsert())
      push_back(entry);
    else {
      iterator it = std::upper_bound(begin(), end(), entry, Comparer<C>(conf));
      insert(it, entry);
    }
  }

  template<class C, class EE>
    void KDTreeLeaf<C, EE>::clear() {
      while (!empty())
        pop_back();
  }

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
    template<class Iter>
    KDTreeLeaf<C, EE>* KDTreeLeaf<C, EE>::makeLeafCopy
    (Iter begin, Iter end,
     memt::Arena& arena, const DivMaskCalculator& calc, const C& conf) {
    MATHIC_ASSERT(static_cast<size_t>(std::distance(begin, end)) <=
                  conf.getLeafSize());
    Leaf* leaf = new (arena.allocObjectNoCon<Leaf>()) Leaf(arena, conf);
    // cannot directly copy as memory is not constructed.
    for (; begin != end; ++begin)
      leaf->push_back(EE(*begin, calc, conf));
    if (conf.getSortOnInsert())
      std::sort(leaf->begin(), leaf->end(), Comparer<C>(conf));
    return leaf;
  }

  template<class C, class EE>
  template<class EM, class MO>
  size_t KDTreeLeaf<C, EE>::removeMultiples
    (const EM& monomial, MO& out, const C& conf) {
    iterator it = begin();
    iterator oldEnd = end();
    for (; it != oldEnd; ++it) {
      if (monomial.divides(*it, conf)) {
        out.push_back(it->get());
        break;
      }
    }
    if (it == oldEnd)
      return 0;
    iterator newEnd = it;
    for (++it; it != oldEnd; ++it) {
      if (!monomial.divides(*it, conf)) {
        *newEnd = *it;
        ++newEnd;
      } else
        out.push_back(it->get());
    }
    // cannot just adjust _end as the superfluous
    // entries at the end need to be destructed.
    const size_t newSize = std::distance(begin(), newEnd);
    const size_t removedCount = size() - newSize;
    MATHIC_ASSERT(newSize < size());
    do {
      pop_back();
    } while (newSize < size());
    MATHIC_ASSERT(size() == newSize);
    return removedCount;
  }

  template<class C, class EE>
    template<class EM>
    typename KDTreeLeaf<C, EE>::iterator
    KDTreeLeaf<C, EE>::findDivisor(const EM& extMonomial, const C& conf) {
    if (!conf.getSortOnInsert()) {
      const iterator stop = end();
      for (iterator it = begin(); it != stop; ++it)
        if (it->divides(extMonomial, conf))
          return it;
      return stop;
    } else {
      iterator rangeEnd =
        std::upper_bound(begin(), end(), extMonomial, Comparer<C>(conf));
      iterator it = begin();
      for (; it != rangeEnd; ++it)
        if (it->divides(extMonomial, conf))
          return it;
      return end();
    }
  }

  template<class C, class EE>
  template<class EM, class DO>
  bool KDTreeLeaf<C, EE>::
    findAllDivisors(const EM& extMonomial, DO& out, const C& conf) {
    if (!conf.getSortOnInsert()) {
      const iterator stop = end();
      for (iterator it = begin(); it != stop; ++it)
        if (it->divides(extMonomial, conf))
          if (!out.proceed(it->get()))
            return false;
    } else {
      iterator rangeEnd =
        std::upper_bound(begin(), end(), extMonomial, Comparer<C>(conf));
      iterator it = begin();
      for (; it != rangeEnd; ++it)
        if (it->divides(extMonomial, conf))
          if (!out.proceed(it->get()))
            return false;
    }
    return true;
  }

  template<class C, class EE>
  template<class EO>
  bool KDTreeLeaf<C, EE>::forAll(EO& output) {
    const iterator stop = end();
    for (iterator it = begin(); it != stop; ++it)
      if (!output.proceed(it->get()))
        return false;
    return true;
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
    MATHIC_ASSERT(size() >= 2);
    // MATHIC_ASSERT not all equal
    Leaf& other = *new (arena.allocObjectNoCon<Leaf>()) Leaf(arena, conf);
    size_t var = (parent == 0 ? static_cast<size_t>(-1) : parent->getVar());
    typename C::Exponent exp;
    while (true) {
      var = (var + 1) % conf.getVarCount();

      if (1) {
        typename C::Exponent min = conf.getExponent(front().get(), var);
        typename C::Exponent max = conf.getExponent(front().get(), var);
        for (iterator it = begin(); it != end(); ++it) {
          min = std::min(min, conf.getExponent(it->get(), var));
          max = std::max(max, conf.getExponent(it->get(), var));
        }
        if (min == max && size() > 1) {
          // todo: avoid infinite loop if all equal
          continue;
        }
        exp = min + (max - min) / 2; // this formula for avg avoids overflow

        iterator newEnd = begin();
        for (iterator it = begin(); it != end(); ++it) {
          if (exp < conf.getExponent(it->get(), var))
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
        ExpOrder<C, EE> order(var, conf);

        std::nth_element(begin(), middle, end(), order);
        if (middle != end()) {
          exp = conf.getExponent(middle->get(), var);
          while (middle != end() && conf.getExponent(middle->get(), var) == exp)
            ++middle;
        }
        if (middle == end() && size() > 1) {
          // todoL avoid infinite loop if all equal.
          continue; // bad split, use another variable
        }
        MATHIC_ASSERT(middle != end());
        MATHIC_ASSERT(exp != conf.getExponent(middle->get(), var));

#ifdef MATHIC_DEBUG
        for (iterator it = begin(); it != middle; ++it) {
          MATHIC_ASSERT(!(exp < conf.getExponent(it->get(), var)));
        }
        for (iterator it = middle; it != end(); ++it) {
          MATHIC_ASSERT(!(conf.getExponent(it->get(), var) < exp));
        }
#endif
        // nth_element does not guarantee where equal elements go,
        // so we cannot just copy [middle, end()).
        iterator newEnd = begin();
        for (iterator it = begin(); it != end(); ++it) {
          if (exp < conf.getExponent(it->get(), var))
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
      MATHIC_ASSERT(other.size() < conf.getLeafSize());
      MATHIC_ASSERT(size() < conf.getLeafSize());
      break;
    }

    Interior& interior = *new (arena.allocObjectNoCon<Interior>())
      Interior(*this, other, var, exp);
    if (C::UseTreeDivMask) {
      KDTreeNode<C,EE>::resetDivMask();
      for (const_iterator it = begin(); it != end(); ++it)
        updateToLowerBound(*it);
      interior.updateToLowerBound(*this);
      interior.updateToLowerBound(other);
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

  /** An object that supports queries for divisors of a monomial using
      a KD Tree (K Dimensional Tree). See DivFinder.h for more documentation.

      Extra fields for Configuration:

      * bool getSortOnInsert() const
      Return true to keep the monomials in leaves sorted to speed up queries.

      * size_t getLeafSize() const
      Return the fixed maximal size of a leaf.
  */
  template<class Configuration>
  class KDTree;

  template<class C>
  class KDTree {
  private:
    typedef KDTreeX<C> Tree;
    typedef typename Tree::ExtEntry ExtEntry;
    typedef typename Tree::ExtMonoRef ExtMonoRef;
  public:
    /** Constructs an object with the given configuration. The configuration
        is copied into the object, so a reference to the passed-in object is
        not kept. The configuration is not copied other than the initial copy. */
    KDTree(const C& configuration):
      _divMaskCalculator(configuration),
      _tree(configuration),
      _size(0) {
      resetNumberOfChangesTillRebuild();
      if (getConfiguration().getUseDivisorCache())
        _divisorCache = 0;
    }

    static const bool UseDivMask = C::UseDivMask;
    typedef typename C::Monomial Monomial;
    typedef typename C::Entry Entry;
    typedef typename C::Exponent Exponent;

    /** Returns whether there are any entries. */
    bool empty() const {return size() == 0;}

    /** Returns the number of entries. */
    size_t size() const {return _size;}

    /** Returns a string that describes the data structure. */
    std::string getName() const;

    /** Returns a reference to this object's configuration object. */
    C& getConfiguration() {
      return _tree.getConfiguration();
    }

    /** Returns a reference to this object's configuration object. */
    const C& getConfiguration() const {
      return const_cast<KDTreeX<C>&>(_tree).getConfiguration();
    }

    /** Removes all multiples of monomial. A duplicate counts
        as a multiple. Returns true if any multiples were removed. */
    bool removeMultiples(const Monomial& monomial) {
      DummyMultipleOutput out;
      return removeMultiples(monomial, out);
    }

    /** Removes all multiples of monomial. A duplicate counts
        as a multiple. Returns true if any multiples were removed.
        Calls out.push_back(entry) for each entry that is removed. */
    template<class MultipleOutput>
    bool removeMultiples(const Monomial& monomial, MultipleOutput& out) {
      ExtMonoRef extMonomial(monomial, _divMaskCalculator, getConfiguration());
      size_t removedCount = _tree.removeMultiples(extMonomial, out);
      reportChanges(0, removedCount);
      return removedCount > 0;
    }

    /** Inserts entry into the data structure. Does NOT remove multiples
        of entry and entry is inserted even if it is a multiple of another
        entry. */
    void insert(const Entry& entry) {
      ExtEntry extEntry(entry, _divMaskCalculator, getConfiguration());
      _tree.insert(extEntry);
      reportChanges(1, 0);
    }

    /** Inserts the entries in the range [begin, end) into the data
        structure. Does NOT remove multiples of entry and entry is inserted
        even if it is a multiple of another entry.

        The elements in the range [begin, end) may be rearranged by this
        function, so the range must be mutable. If that is not acceptable,
        call the one element insert method for each element. */
    template<class Iter>
    void insert(Iter begin, Iter end) {
      if (begin == end)
        return;
      const size_t inserted = std::distance(begin, end); 
      if (!empty()) {
        for (; begin != end; ++begin)
          _tree.insert(*begin);
      } else {
        _tree.insert(begin, end);
        // insert into empty container is equivalent to rebuild
        resetNumberOfChangesTillRebuild();
      }
      reportChanges(inserted, 0);
    }

    /** Returns a pointer to an entry that divides monomial. Returns null if no
        entries divide monomial. */
    Entry* findDivisor(const Monomial& monomial) {
      // todo: do this on extended monomials. requires cache to be extended.
      const C& conf = getConfiguration();
      if (conf.getUseDivisorCache() &&
        _divisorCache != 0 &&
        conf.divides(*_divisorCache, monomial))
        return _divisorCache;

      ExtMonoRef extMonomial(monomial, _divMaskCalculator, getConfiguration());
      Entry* divisor = _tree.findDivisor(extMonomial);
      if (conf.getUseDivisorCache() && divisor != 0)
        _divisorCache = divisor;
      return divisor;
    }

    /** Returns the position of a divisor of monomial. Returns end() if no
        entries divide monomial. */
    const Entry* findDivisor(const Monomial& monomial) const {
      return const_cast<KDTree<C>&>(*this).findDivisor(monomial);
    }

    /** Calls out.proceed(entry) for each entry that divides monomial.
        The method returns if proceed returns false, otherwise the
        search for divisors proceeds. */
    template<class DivisorOutput>
    void findAllDivisors(const Monomial& monomial, DivisorOutput& out) {
      ExtMonoRef extMonomial(monomial, _divMaskCalculator, getConfiguration());
      _tree.findAllDivisors(extMonomial, out);
    }

    /** Calls output.proceed(entry) for each entry that divides monomial.
        The method returns if proceed returns false, otherwise the
        search for divisors proceeds. */
    template<class DivisorOutput>
    void findAllDivisors(const Monomial& monomial, DivisorOutput& output) const {
      ConstEntryOutput<DivisorOutput> constOutput(output);
      const_cast<KDTree<C>&>(*this).findAllDivisors(monomial, constOutput);
    }

    /** Calls output.proceed(entry) for each entry.
        The method returns if proceed returns false, otherwise the
        search for divisors proceeds. */
    template<class EntryOutput>
    void forAll(EntryOutput& out) {
      _tree.forAll(out);
    }

    /** Calls out.proceed(entry) for each entry.
        The method returns if proceed returns false, otherwise the
        search for divisors proceeds. */
    template<class EntryOutput>
    void forAll(EntryOutput& output) const {
      ConstEntryOutput<EntryOutput> constOutput(output);
      const_cast<KDTree<C>&>(*this).forAll(constOutput);
    }

    /** Removes all entries. Does not reset the configuration object. */
    void clear() {
      _tree.clear();
      resetNumberOfChangesTillRebuild();
      _tree.calc().rebuildDefault(getConfiguration());
    }

    /** Rebuilds the data structure. */
    void rebuild() {
      EntryRecorder recorder(memt::Arena::getArena(), size());
      _tree.forAll(recorder);
      _tree.clear();
      _divMaskCalculator.rebuild
        (recorder.begin(), recorder.end(), getConfiguration());
      _tree.insert(recorder.begin(), recorder.end(), _divMaskCalculator);
      resetNumberOfChangesTillRebuild();
    }

	/** Returns the number of bytes allocated by this object. Does not
		include sizeof(*this), does not include any additional memory
		that the configuration may have allocated and does not include
		any memory that an Entry may point to. Does include
		sizeof(Entry) as well as unused memory that is being kept to
		avoid frequent allocations. */
    size_t getMemoryUsage() const {return _tree.getMemoryUsage();}

  private:
    KDTree(const KDTree<C>&); // unavailable
    void operator=(const KDTree<C>&); // unavailable

    // For recording all entries in the tree using forAll.
    class EntryRecorder {
    public:
      EntryRecorder(memt::Arena& arena, size_t capacity):
        _entries(arena, capacity) {}
      bool proceed(const Entry& entry) {
        _entries.push_back(entry);
       return true;
      }
      Entry* begin() {return _entries.begin();}
      Entry* end() {return _entries.end();}

    private:
      memt::ArenaVector<Entry, true> _entries;
    };

    /// makes the parameter given to proceed be const.
    template<class EntryOutput>
    class ConstEntryOutput {
    public:
    ConstEntryOutput(EntryOutput& out): _out(out) {}
      bool proceed(const Entry& entry) {return _out.proceed(entry);}
    private:
      EntryOutput& _out;
    };

    /// Ignores everything passed to it. */
    class DummyMultipleOutput {
    public:
      void push_back(Entry& e) {}
    };

    void reportChanges(size_t additions, size_t removals);
    void resetNumberOfChangesTillRebuild();
    bool reportChangesRebuild(size_t additions, size_t removals);
    size_t _changesTillRebuild; /// Update using reportChanges().

    Entry* _divisorCache; /// The divisor in the previous query. Can be null.

    // All DivMasks calculated using this.
    typename Tree::DivMaskCalculator _divMaskCalculator;
    KDTreeX<C> _tree;
    size_t _size;
  };

  template<class C>
  std::string KDTree<C>::getName() const {
    std::stringstream out;
    const C& conf = getConfiguration();
    out << "KDTree leaf:" << conf.getLeafSize();
    if (UseDivMask && conf.getDoAutomaticRebuilds()) {
      out << " autob:" << conf.getRebuildRatio()
          << '/' << conf.getRebuildMin();
    }
    out << (C::UseDivMask && !C::UseTreeDivMask ? " dmask" : "")
        << (C::UseTreeDivMask ? " tree-dmask" : "")
        << (conf.getSortOnInsert() ? " sort" : "")
        << (conf.getUseDivisorCache() ? " cache" : "");
    return out.str();
  }

  template<class C>
  void KDTree<C>::resetNumberOfChangesTillRebuild() {
    const C& conf = getConfiguration();
    if (conf.getUseDivisorCache())
      _divisorCache = 0;
    if (!conf.getDoAutomaticRebuilds())
      return;
    MATHIC_ASSERT(conf.getRebuildRatio() > 0);
    _changesTillRebuild = std::max
      (static_cast<size_t>(size() * conf.getRebuildRatio()),
       conf.getRebuildMin());
  }

  template<class C>
  void KDTree<C>::reportChanges(size_t additions, size_t removals) {
    if (getConfiguration().getUseDivisorCache() && (additions | removals) != 0)
      _divisorCache = 0;
    if (reportChangesRebuild(additions, removals))
      rebuild();
  }

  template<class C>
  bool KDTree<C>::reportChangesRebuild(size_t additions, size_t removals) {
    // note how negative value/overflow of _changesTillRebuild cannot
    // happen this way.
    MATHIC_ASSERT(removals <= size() + additions);
    _size = (size() + additions) - removals;
    if (!getConfiguration().getDoAutomaticRebuilds())
      return false;
    const size_t changesMadeCount = additions + removals;
    if (_changesTillRebuild > changesMadeCount) {
      _changesTillRebuild -= changesMadeCount;
      return false;
    } else
      return true;
  }
}

#endif
