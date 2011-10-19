#ifndef MATHIC_BIT_MASK_GUARD
#define MATHIC_BIT_MASK_GUARD

#include <vector>
#include <utility>

namespace mathic {
  /** Class representing a div mask. This is a set of bits that can
      be used to determine that one monomial cannot divide another
      monomial. */
  class DivMask {
  public:
    /** Calculates div masks. Don't change NullCalculator
        from its default value. The actual code are in partial specializations
        selecting the right version based on NullCalculator. */
    template<class Configuration,
      bool NullCalculator = Configuration::UseDivMask>
      class Calculator;

  DivMask(): _mask(0) {}
    static DivMask getMaxMask() {return ~static_cast<MaskType>(0);}

    template<class T, class Configuration>
      DivMask(const T& t,
              const Calculator<Configuration>& calc, const Configuration& conf):
    _mask(calc.compute(t, conf)) {}

    template<class T, class Configuration>
      void recalculate(const T& t,
                       const Calculator<Configuration>& calc,
                       const Configuration& conf) {
      _mask = calc.compute(t, conf);
    }

    bool canDivide(const DivMask& mask) const {
      return (_mask & ~mask._mask) == 0;
    }

    void combineAnd(const DivMask& mask) {_mask &= mask._mask;}

    bool operator==(DivMask& mask) const {return _mask == mask._mask;}
    bool operator!=(DivMask& mask) const {return _mask != mask._mask;}

    /** Extender extends T with a div mask if UseDivMask is true. It is
        allowed for T to be a reference type or const. */
    template<class T, bool UseDivMask>
      class Extender;

    /** Base class to include a DivMask into a class at compile time
        based on the template parameter UseDivMask. The class offers
        the same methods either way, but they are replaced by do-nothing
        or asserting versions if UseDivMask is false. */
    template<bool UseDivMask>
      class HasDivMask;

  protected:
    typedef unsigned long MaskType;
  DivMask(MaskType mask): _mask(mask) {}

  private:
    /** To eliminate warnings about T& if T is already a reference type. */
    template<class T> struct Ref {typedef T& RefType;};
    template<class T> struct Ref<T&> {typedef T& RefType;};
    MaskType _mask;
  };

  template<class C>
    class DivMask::Calculator<C, true> {
  public:
    Calculator(const C& conf) {rebuildDefault(conf);}

    /** Change the meaning of the bits in the div masks produced by this
        object to something that will likely work well for entries
        in [begin, end). All div masks will have to be recomputed
        after this. Mixing div masks computed before a call to
        rebuild() with ones after has unpredictable results. */
    template<class Iter>
      void rebuild(Iter begin, Iter end, const C& conf);

    /** Rebuilds without the benefit of knowing a range of entries
        that the div masks are supposed to work well for. */
    void rebuildDefault(const C& conf);

    /** Computes a div mask for t. */
    template<class T>
      DivMask::MaskType compute(const T& t, const C& conf) const;

  private:
    typedef typename C::Exponent Exponent;
    /** If entry at index i is the pair (var,exp) then the bit at
        index var in a bit mask is 1 if the exponent of var in the
        monomial is strictly greater than exp. */
    typedef std::vector<std::pair<size_t, Exponent> > BitContainer;
    BitContainer _bits;
  };

  template<class C>
    template<class Iter>
    void DivMask::Calculator<C, true>::
    rebuild(Iter begin, Iter end, const C& conf) {
    if (begin == end) {
      rebuildDefault(conf);
      return;
    }

    _bits.clear();
    const size_t varCount = conf.getVarCount();
    const size_t TotalBits = sizeof(MaskType) * BitsPerByte;
    for (size_t var = 0; var < varCount; ++var) {
      const size_t bitsForVar =
        TotalBits / varCount + (var < TotalBits % varCount);
      if (bitsForVar == 0)
        continue;

      // determine minimum and maximum
      Exponent min = conf.getExponent(*begin, 0);
      Exponent max = min;
      for (Iter it = begin; it != end; ++it) {
        if (max < conf.getExponent(*begin, var))
          max = conf.getExponent(*begin, var);
        if (conf.getExponent(*begin, var) < min)
          min = conf.getExponent(*begin, var);
      }

      // divide the range [a,b] into bitsForVar equal pieces
      // and use the left end points of those ranges
      // as the points for the bits.
      Exponent increment = (max - min) / bitsForVar;
      if (increment == 0)
        increment = 1;
      for (size_t i = 0; i < bitsForVar; ++i)
        _bits.push_back(std::make_pair(var, min + increment * i));
    }
  }

  template<class C>
    void DivMask::Calculator<C, true>::
    rebuildDefault(const C& conf) {
    _bits.clear();
    const size_t varCount = conf.getVarCount();
    const size_t TotalBits = sizeof(MaskType) * BitsPerByte;
    for (size_t var = 0; var < varCount; ++var) {
      const size_t bitsForVar =
        TotalBits / varCount + (var < TotalBits % varCount);
      Exponent exp = 0;
      for (size_t i = 0; i < bitsForVar; ++i) {
        _bits.push_back(std::make_pair(var, exp));
        exp = (i == 0 ? 1 : exp * 2);
      }
    }
  }

  template<class C>
    template<class T>
    DivMask::MaskType DivMask::Calculator<C, true>::
    compute(const T& t, const C& conf) const {
    typedef typename BitContainer::const_iterator const_iterator;
    const const_iterator end = _bits.end();
    DivMask::MaskType mask = 0;
    for (const_iterator it = _bits.begin(); it != end; ++it)
      mask = (mask << 1) | (conf.getExponent(t, it->first) > it->second);
    return mask;
  }

  template<class C>
    class DivMask::Calculator<C, false> {
  public:
    Calculator(const C& conf) {}

    template<class Iter>
      void rebuild(Iter begin, Iter end, const C& conf) {}
    void rebuildDefault(const C& conf) {}
  };

  template<>
    class DivMask::HasDivMask<true> {
  public:
    template<class T, class C>
      HasDivMask(const T& t, const Calculator<C>& calc, const C& conf):
    _mask(t, calc, conf) {}
    HasDivMask() {resetDivMask();}

    DivMask& getDivMask() {return _mask;}
    const DivMask& getDivMask() const {return _mask;}
    void resetDivMask() {_mask = DivMask::getMaxMask();}
    bool canDivide(const HasDivMask<true>& t) const {
      return getDivMask().canDivide(t.getDivMask());
    }

    void updateToLowerBound(const HasDivMask<true>& t) {
      _mask.combineAnd(t.getDivMask());
    }

  protected:
    template<class T, class C>
      void recalculateDivMask
      (const T& t, const Calculator<C>& calc, const C& conf) {
      _mask.recalculate(t, calc, conf);
    }

  private:
    DivMask _mask;
  };

  template<>
    class DivMask::HasDivMask<false> {
  public:
    void resetDivMask() {MATHIC_ASSERT(false);}
    DivMask getDivMask() const {MATHIC_ASSERT(false); return DivMask();}
    bool canDivide(const HasDivMask<false>& t) const {return true;}
    template<bool B>
      void updateToLowerBound(const HasDivMask<B>& entry) {}
  };

  template<class T>
    class DivMask::Extender<T, true> : public HasDivMask<true> {
  private:
    typedef typename Ref<T>::RefType Reference;
    typedef typename Ref<const T>::RefType ConstReference;
  public:
  Extender(): HasDivMask<true>(), _t() {}
    template<class C>
      Extender(ConstReference t, const Calculator<C>& calc, const C& conf):
    HasDivMask<true>(t, calc, conf), _t(t) {}

    template<class S, class C>
      bool divides(const Extender<S, true>& t, const C& conf) const {
      return this->canDivide(t) && conf.divides(get(), t.get());
    }

    template<class C>
      void recalculateDivMask(const Calculator<C>& calc, const C& conf) {
      this->HasDivMask<true>::recalculateDivMask(get(), calc, conf);
    }

    Reference get() {return _t;}
    ConstReference get() const {return _t;}

  private:
    T _t;
  };

  template<class T>
    class DivMask::Extender<T, false> : public HasDivMask<false> {
  private:
    typedef typename Ref<T>::RefType Reference;
    typedef typename Ref<const T>::RefType ConstReference;
  public:
  Extender(): HasDivMask<false>(), _t() {}
    template<class C>
      Extender(ConstReference t, const Calculator<C>& calc, const C& conf):
    HasDivMask<false>(), _t(t) {}

    template<class S, class C>
      bool divides(const Extender<S, false>& t, const C& conf) const {
      return conf.divides(get(), t.get());
    }

    template<class C>
      void recalculateDivMask(const Calculator<C>& calc, const C& conf) {}

    Reference get() {return _t;}
    ConstReference get() const {return _t;}

  private:
    T _t;
  };
}

#endif
