#ifndef MATHIC_K_D_TREE_GUARD
#define MATHIC_K_D_TREE_GUARD

#include "stdinc.h"
#include "KDTreeLeaf.h"
#include "DivMask.h"
#include "memtailor/memtailor.h"
#include <list>
#include <string>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <vector>

namespace mathic {
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
  public:
    typedef C Configuration;
    static const bool UseDivMask = C::UseDivMask;
    typedef typename C::Monomial Monomial;
    typedef typename C::Entry Entry;
    typedef typename C::Exponent Exponent;

  private:
    typedef typename DivMask::Extender<Entry, C::UseDivMask> ExtEntry;
    typedef typename DivMask::Extender<const Monomial&,C::UseDivMask> ExtMonoRef;
    typedef typename DivMask::Calculator<C> DivMaskCalculator;

    typedef KDTreeNode<C, ExtEntry> Node;
    typedef KDTreeInterior<C, ExtEntry> Interior;
    typedef KDTreeLeaf<C, ExtEntry> Leaf;

    typedef typename Leaf::iterator LeafIt;

  public:
    /** Constructs an object with the given configuration. The configuration
        is copied into the object, so a reference to the passed-in object is
        not kept. The configuration is not copied other than the initial copy. */
    KDTree(const C& configuration);
    ~KDTree();

    /** Removes all multiples of monomial. A duplicate counts
        as a multiple. Returns true if any multiples were removed. */
    bool removeMultiples(const Monomial& monomial);

    /** Removes all multiples of monomial. A duplicate counts
        as a multiple. Returns true if any multiples were removed.
        Calls out.push_back(entry) for each entry that is removed. */
    template<class MultipleOutput>
    bool removeMultiples(const Monomial& monomial, MultipleOutput& out);

    /** Inserts entry into the data structure. Does NOT remove multiples
        of entry and entry is inserted even if it is a multiple of another
        entry. */
    void insert(const Entry& entry);

    /** Inserts the entries in the range [begin, end) into the data
        structure. Does NOT remove multiples of entry and entry is inserted
        even if it is a multiple of another entry.

        The elements in the range [begin, end) may be rearranged by this
        function, so the range must be mutable. If that is not acceptable,
        call the one element insert method for each element. */
    template<class Iter>
    void insert(Iter begin, Iter end);

    /** Returns a pointer to an entry that divides monomial. Returns null if no
        entries divide monomial. */
    Entry* findDivisor(const Monomial& monomial);

    /** Returns the position of a divisor of monomial. Returns end() if no
        entries divide monomial. */
    const Entry* findDivisor(const Monomial& monomial) const {
      return const_cast<KDTree<C>&>(*this).findDivisor(monomial);
    }

    /** Calls out.proceed(entry) for each entry that divides monomial.
        The method returns if proceed returns false, otherwise the
        search for divisors proceeds. */
    template<class DivisorOutput>
    void findAllDivisors(const Monomial& monomial, DivisorOutput& out);

    /** Calls out.proceed(entry) for each entry that divides monomial.
        The method returns if proceed returns false, otherwise the
        search for divisors proceeds. */
    template<class DivisorOutput>
    void findAllDivisors(const Monomial& monomial, DivisorOutput& out) const;

    /** Calls out.proceed(entry) for each entry.
        The method returns if proceed returns false, otherwise the
        search for divisors proceeds. */
    template<class EntryOutput>
    void forAll(EntryOutput& out);

    /** Calls out.proceed(entry) for each entry.
        The method returns if proceed returns false, otherwise the
        search for divisors proceeds. */
    template<class EntryOutput>
    void forAll(EntryOutput& out) const;

    /** Removes all entries. Does not reset the configuration object. */
    void clear();

    /** Rebuilds the data structure. */
    void rebuild();

    /** Returns whether there are entries. */
    bool empty() const {return size() == 0;}

    /** Returns the number of entries. */
    size_t size() const {return _size;}

    /** Returns a string that describes the data structure. */
    std::string getName() const;

    /** Returns a reference to this object's configuration object. */
    C& getConfiguration() {return _conf;}

    /** Returns a reference to this object's configuration object. */
    const C& getConfiguration() const {return _conf;}

	/** Returns the number of bytes allocated by this object. Does not
		include sizeof(*this), does not include any additional memory
		that the configuration may have allocated and does not include
		any memory that an Entry may point to. Does include
		sizeof(Entry) as well as unused memory that is being kept to
		avoid frequent allocations. */
    size_t getMemoryUsage() const;

  private:
    KDTree(const KDTree<C>&); // unavailable
    void operator=(const KDTree<C>&); // unavailable

    void resetNumberOfChangesTillRebuild();
    void reportChanges(size_t additions, size_t removals);

    /** Transfers proceed from entry to const entry. */
    template<class DO>
      class ConstEntryOutput {
    public:
    ConstEntryOutput(DO& out): _out(out) {}
      bool proceed(const Entry& entry) {return _out.proceed(entry);}
    private:
      DO& _out;
    };

    /** Ignores everything passed to it. */
    class DummyMultipleOutput {
    public:
      void push_back(Entry& e) {}
    };

    template<class Iter>
      struct InsertTodo {
        Iter begin;
        Iter end;
        Interior* parent;
      };

#ifdef MATHIC_DEBUG
    bool debugIsValid() const;
#endif

    mutable std::vector<Node*> _tmp; // For navigating the tree.
    memt::Arena _arena; // Everything permanent allocated from here.
    C _conf; // User supplied configuration.
    Node* _root; // Root of the tree.
    size_t _size; // Number of entries.
    Entry* _divisorCache; /// The divisor in the previous query. Can be null.
    size_t _changesTillRebuild; /// Update using reportChanges().
    DivMaskCalculator _divMaskCalculator; // All DivMasks calculated using this.
  };

  template<class C>
    KDTree<C>::KDTree(const C& configuration):
  _conf(configuration),
    _size(0),
    _divMaskCalculator(configuration) {
      MATHIC_ASSERT(_conf.getLeafSize() >= 2);
      _root = new (_arena.allocObjectNoCon<Leaf>()) Leaf(_arena, _conf);
      if (_conf.getUseDivisorCache())
        _divisorCache = 0;
      resetNumberOfChangesTillRebuild();
      MATHIC_ASSERT(debugIsValid());
    }

  template<class C>
  KDTree<C>::~KDTree() {
    clear();
  }

  template<class C>
  std::string KDTree<C>::getName() const {
    std::stringstream out;
    out << "KDTree leaf:" << _conf.getLeafSize();
    if (UseDivMask && _conf.getDoAutomaticRebuilds()) {
      out << " autob:" << _conf.getRebuildRatio()
          << '/' << _conf.getRebuildMin();
    }
    out << (C::UseDivMask && !C::UseTreeDivMask ? " dmask" : "")
        << (C::UseTreeDivMask ? " tree-dmask" : "")
        << (_conf.getSortOnInsert() ? " sort" : "")
        << (_conf.getUseDivisorCache() ? " cache" : "");
    return out.str();
  }

  template<class C>
    template<class MO>
    bool KDTree<C>::removeMultiples(const Monomial& monomial, MO& out) {
    ExtMonoRef extMonomial(monomial, _divMaskCalculator, _conf);

    MATHIC_ASSERT(_tmp.empty());
    size_t removedCount = 0;
    Node* node = _root;
    while (true) {
      while (node->isInterior()) {
        Interior& interior = node->asInterior();
        if (!(interior.getExponent() <
              _conf.getExponent(monomial, interior.getVar())))
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
    if (_conf.getUseDivisorCache() && removedCount > 0)
      _divisorCache = 0;
    reportChanges(0, removedCount);
    MATHIC_ASSERT(debugIsValid());
    return removedCount > 0;
  }

  template<class C>
    bool KDTree<C>::removeMultiples(const Monomial& monomial) {
    DummyMultipleOutput out;
    return removeMultiples(monomial, out);
  }

  template<class C>
  void KDTree<C>::insert(const Entry& entry) {
    ExtEntry extEntry(entry, _divMaskCalculator, _conf);

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

    if (_conf.getUseDivisorCache())
      _divisorCache = 0;
    reportChanges(1, 0);
    MATHIC_ASSERT(debugIsValid());
  }

  /// @todo: this function is too big and it knows too much about the details
  /// inside nodes. Also, it allocates a std::vector every time.
  template<class C>
    template<class Iter>
    void KDTree<C>::insert(Iter insertBegin, Iter insertEnd) {
    if (!empty()) {
      for (; insertBegin != insertEnd; ++insertBegin)
        insert(*insertBegin);
      return;
    } else if (insertBegin == insertEnd)
      return;

    _arena.freeAllAllocs();
    _root = 0;
    _size = std::distance(insertBegin, insertEnd);
    _divMaskCalculator.rebuild(insertBegin, insertEnd, _conf);
    if (_conf.getUseDivisorCache())
      _divisorCache = 0;

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
          (insertBegin, insertEnd, _arena, _divMaskCalculator, _conf);
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
    if (_conf.getUseDivisorCache())
      _divisorCache = 0;

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
    // range insert into empty container is equivalent to a rebuild.
    resetNumberOfChangesTillRebuild();
  }

  template<class C>
  typename KDTree<C>::Entry* KDTree<C>::findDivisor
    (const Monomial& monomial) {
    if (_conf.getUseDivisorCache() &&
        _divisorCache != 0 && _conf.divides(*_divisorCache, monomial))
      return _divisorCache;

    ExtMonoRef extMonomial(monomial, _divMaskCalculator, _conf);

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
            _conf.getExponent(monomial, interior.getVar()))
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
          if (_conf.getUseDivisorCache())
            _divisorCache = &leafIt->get();
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
    void KDTree<C>::findAllDivisors(const Monomial& monomial, DO& output) {
    ExtMonoRef extMonomial(monomial, _divMaskCalculator, _conf);

    MATHIC_ASSERT(_tmp.empty());
    Node* node = _root;
    while (true) {
      while (node->isInterior()) {
        Interior& interior = node->asInterior();
        if (interior.getExponent() <
            _conf.getExponent(monomial, interior.getVar()))
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
  template<class DO>
  void KDTree<C>::findAllDivisors(const Monomial& monomial, DO& output) const {
    ConstEntryOutput<DO> constOutput(output);
    const_cast<KDTree<C>&>(*this).findAllDivisors(monomial, constOutput);
  }

  template<class C>
  template<class EO>
  void KDTree<C>::forAll(EO& output) {
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
  template<class DO>
  void KDTree<C>::forAll(DO& output) const {
    ConstEntryOutput<DO> constOutput(output);
    const_cast<KDTree<C>&>(*this).forAll(constOutput);
  }

  template<class C>
  void KDTree<C>::clear() {
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
    _size = 0;
    _root = new (_arena.allocObjectNoCon<Leaf>()) Leaf(_arena, _conf);
    resetNumberOfChangesTillRebuild();
    _divMaskCalculator.rebuildDefault(_conf);
    if (_conf.getUseDivisorCache())
      _divisorCache = 0;
  }

  template<class C>
    void KDTree<C>::rebuild() {
    const size_t totalSize = size();
    typedef memt::ArenaVector<Entry, true> TmpContainer;
    TmpContainer tmpCopy(memt::Arena::getArena(), totalSize);

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
      Leaf& leaf = node->asLeaf();
      typename Leaf::const_iterator stop = leaf.end();
      for (typename Leaf::const_iterator it = leaf.begin(); it != stop; ++it)
        tmpCopy.push_back(it->get());
      leaf.clear();
    }
    _arena.freeAllAllocs();
    _size = 0;
    _root = new (_arena.allocObjectNoCon<Leaf>()) Leaf(_arena, _conf);
    resetNumberOfChangesTillRebuild();
    // range insert into empty container IS a rebuild.
    insert(tmpCopy.begin(), tmpCopy.end());
  }

  template<class C>
    void KDTree<C>::resetNumberOfChangesTillRebuild() {
    if (!_conf.getDoAutomaticRebuilds())
      return;
    MATHIC_ASSERT(_conf.getRebuildRatio() > 0);
    _changesTillRebuild = std::max
      (static_cast<size_t>(size() * _conf.getRebuildRatio()),
       _conf.getRebuildMin());
  }

  template<class C>
    void KDTree<C>::reportChanges(size_t additions, size_t removals) {
    // note how negative value/overflow of _changesTillRebuild cannot
    // happen this way.
    MATHIC_ASSERT(removals <= _size + additions);
    _size = (_size + additions) - removals;
    if (!_conf.getDoAutomaticRebuilds())
      return;
    const size_t changesMadeCount = additions + removals;
    if (_changesTillRebuild > changesMadeCount)
      _changesTillRebuild -= changesMadeCount;
    else
      rebuild();
  }

  template<class C>
  size_t KDTree<C>::getMemoryUsage() const {
	size_t sum = _arena.getMemoryUsage();
	sum += _tmp.capacity() * sizeof(_tmp.front());
	return sum;
  }

#ifdef MATHIC_DEBUG
  template<class C>
  bool KDTree<C>::debugIsValid() const {
    MATHIC_ASSERT(_tmp.empty());
    MATHIC_ASSERT(!_conf.getDoAutomaticRebuilds() || _conf.getRebuildRatio() > 0);

    MATHIC_ASSERT(_root != 0);
    if (empty())
      return true;

    std::vector<Node*> nodes;
    nodes.reserve(size());
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
    MATHIC_ASSERT(sizeSum == size());

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
}

#endif
