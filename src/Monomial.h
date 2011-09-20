#ifndef MONOMIAL_GUARD
#define MONOMIAL_GUARD

#include <vector>

/*
typedef std::vector<int> Monomial;
//*/
class Monomial {
public:
  typedef int Exponent;

  Monomial(): _exponents(0) {IF_DEBUG(_size = 0);}
  Monomial(std::vector<Exponent>& v): _exponents(&v[0]) {
    IF_DEBUG(_size = v.size());
  }

  inline Exponent& operator[](size_t index) {
    ASSERT(index < _size);
    return _exponents[index];
  }
  inline const Exponent& operator[](size_t index) const {
    ASSERT(index < _size);
    return _exponents[index];
  }

  const Exponent* getPointer() const {return _exponents;}

#ifdef DEBUG
  size_t size() const {return _size;}
  bool operator==(const Monomial& m) const {return _exponents == m._exponents;}
  bool operator<(const Monomial& m) const {return _exponents < m._exponents;}
#endif

private:
#ifdef DEBUG
  size_t _size;
#endif
  Exponent* _exponents;
};
//*/

/*
class Monomial {
public:
  typedef int Exponent;

  Monomial(std::vector<Exponent>& v): _exponents(v) {}
  operator std::vector<Exponent>&() {return _exponents;}
  operator const std::vector<Exponent>&() const {return _exponents;}

  Exponent& operator[](size_t index) {
    ASSERT(index < _exponents.size());
    return _exponents[index];
  }
  const Exponent& operator[](size_t index) const {
    ASSERT(index < _exponents.size());
    return _exponents[index];
  }
  void push_back(const Exponent& exponent) {
    _exponents.push_back(exponent);
  }
  bool operator==(const Monomial& monomial) const {
    return _exponents == monomial._exponents;
  }

private:
  std::vector<Exponent> _exponents;
};
//*/
#endif
