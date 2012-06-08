// This is an example of how to make a divisor query data structure using
// Mathic. We make both a KDTree and a DivList. We use a std::vector<int>
// as a monomial to keep things simple.

#include <mathic/KDTree.h>
#include <mathic/DivList.h>

#include <iostream>
#include <vector>

// Divisor query data structures a parametrized on a configuration type
// that tells the data structure how to represent monomials and how to do
// a few basic operations like extracting an exponent from a monomial.
//
// BasicConfiguration contains the methods necessary to teach Mathic
// how to use a std::vector as a representation for a monomial.
class BasicConfiguration {
public: // the methods that Mathic calls must be public

  // The configuration needs to know how many variables there are
  // so that it can tell Mathic. Variables are indexed from 0 to
  // numberOfVariables - 1.
  BasicConfiguration(size_t numberOfVariables):
    varCount(numberOfVariables) {}

  // The type of an exponent in a monomial.
  typedef int Exponent;

  // The type of monomials.
  typedef std::vector<Exponent> Monomial;

  // The type of monomials that are entries in the data structure.
  // For some applications we might want this to be a completely
  // different type from Monomial, but for this simple example we will
  // just let this be the same type as Monomial.
  typedef Monomial Entry;

  // This method is how Mathic determines how many variables there are.
  // Mathic does not store this number and instead calls this method every time
  // it needs to know this information, so this method should run quickly.
  size_t getVarCount() const {
    return varCount;
  }

  // Mathic uses this function to extract an Exponent from a Monomial.
  Exponent getExponent(Monomial m, size_t var) const {
    MATHIC_ASSERT(var < m.size());
    return m[var];
  }

  // Mathic uses this function to determine if a divides b.
  // This method can be implemented using the other methods on a
  // configuration as done here. The point of this method is that
  // some representations of monomials allow faster divisiblity
  // checks than is achieved by extracting every exponent.
  //
  // We only need one method called divides since Monomial and Entry
  // are the same type. Otherwise we would need 3 versions
  // taking parameters (Monomial, Entry), (Entry, Monomial)
  // and (Monomial,Monomial). Thankfully here we only need
  // (Monomial,Monomial).
  bool divides(const Monomial& a, const Monomial& b) const {
    size_t const size = getVarCount();
    for (size_t var = 0; var < size; ++var)
      if (getExponent(a, var) > getExponent(b, var))
        return false; // a does not divide b
    return true; // a divides b
  }

private:
  size_t varCount; // the number of variables
};

// We can now combine the BasicConfiguration with options for a KDTree
// to get a configuration type for a KDTree.
class KDTreeConfiguration :
  public BasicConfiguration, // sets std::vector as the monomial representation
  public mathic::KDTreeSuggestedOptions // sets default options for the KDTree
{
public:
  KDTreeConfiguration(size_t varCount): BasicConfiguration(varCount) {}
};

// With a configuration type, we can now create a KDTree type that uses
// the std::vector representation and using default options.
class MyKDTree : public mathic::KDTree<KDTreeConfiguration> {
public:
  MyKDTree(size_t varCount): 
    KDTree<KDTreeConfiguration>(KDTreeConfiguration(varCount)) {}
};

// We can do the same thing for DivList
class DivListConfiguration :
  public BasicConfiguration, // sets std::vector as the monomial representation
  public mathic::DivListSuggestedOptions // sets default options for the DivList
{
public:
  DivListConfiguration(size_t varCount): BasicConfiguration(varCount) {}
};

// note how mic:: is an alias for the mathic:: namespace. You can disable
// the mic:: namespace definition by defining the preprocessor
// symbol MATHIC_NO_MIC_NAMESPACE in case you have some other library
// that uses the mic:: namespace.
typedef mic::DivList<DivListConfiguration> MyDivList;

int main() {
  MyKDTree tree(2);
  std::vector<int> a(2);
  a[0] = 3;
  a[1] = 2;
  tree.insert(a);

  std::vector<int> b(2);
  b[0] = 1;
  b[1] = 4;
  tree.insert(b);

  std::vector<int> c(2);
  c[0] = 5;
  c[1] = 0;
  // returns null as there is no divisor
  MATHIC_ASSERT(tree.findDivisor(c) == 0);

  // returns a pointer to a copy of a
  c[1] = 2;
  MATHIC_ASSERT(*tree.findDivisor(c) == a);

  return 0;
}
