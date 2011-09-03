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
  These are the things added to the finder. Entry represents a Monomial
  along with any additional required information that client code needs
  to associate with monomials in the finder.

 * A type Exponent
  These are exponents of the monomials. Exponent must have a copy constructor
  and operator <.

 * A function Exponent getExponent(size_t var, Monomial m) const
  Returns the exponent of the variable var in m. m can be a const reference.

 * A function bool divides(Monomial a, Monomial b) const
  Returns whether a divides b. a and b can be const references.

 * A function size_t getVarCount() const
  Returns the number of variables. Variables are indexed from 0 to
  getVarCount() - 1.
*/
template<class Configuration>
class DivFinder; // no implementation

#endif
