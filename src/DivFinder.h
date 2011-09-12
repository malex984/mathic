#ifndef DIV_FINDER_GUARD
#define DIV_FINDER_GUARD

/** An object that supports queries for divisors of a monomial.

 This is an object for documentation purposes. Use the concrete
 implementations of this interface to get an actual DivFinder-like
 object.

 DivFinder stores configurable Entry objects that represent a monomial
 together with any required additional information. 

 The class is parameterized on a Configuration. The Configuration class
 must have the following members.

 * A type Monomial
  Represents a monomial.

 * A type Entry
  These are the things added to the finder. Entry represents a monomial
  along with any additional required information that client code needs
  to associate with monomials in the finder.

 * A type Exponent
  These are exponents of the monomials. Exponent must have a copy constructor,
  an operator= and an operator< for comparision with other exponents.

 * A function Exponent getExponent(Monomial m, size_t var) const
  Returns the exponent of the variable var in m. m can be a const reference.

 * A function bool divides(Entry a, Monomial b) const
 * A function bool divides(Monomial a, Entry b) const
 * A function bool divides(Entry a, Entry b) const
  Returns whether a divides b. a and b can be const references. If Entry
  and Monomial are the same type then only one function is needed.

 * A function size_t getVarCount() const
  Returns the number of variables. Variables are indexed from 0 to
  getVarCount() - 1. Must be a positive number.

 DivFinder's that make use of sorting additionally require these fields:

 * A type Comparer with a function bool operator()(Entry a, Monomial b) const
 * A type Comparer with a function bool operator()(Monomial a, Entry b) const
 * A type Comparer with a function bool operator()(Entry a, Entry b) const
  Comparer must define a total order < on monomials. a and b can be const
  references. If Entry and Monomial are the same type then only
  one function is needed.

 * A function Comparer getComparer() const
  Returns a Comparer object. All Comparer objects must define the same
  ordering.

 It is possible to obtain non-const Entries from a DivFinder. It is allowed
 to change these, but the monomial they represent must not change. When the
 DivFinder calls some method X on the Configuration, the Configuration must
 not in return cause a method to be called on the DivFinder until the method
 X has returned. This is because the DivFinder may not be in a valid state
 at the time it is calling X.
*/
template<class Configuration>
class DivFinder; // no implementation

#endif
