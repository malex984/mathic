#ifndef MATHIC_BIT_TRIANGLE_GUARD
#define MATHIC_BIT_TRIANGLE_GUARD

#include "stdinc.h"
#include <vector>

namespace mathic {
  // Object that stores a triangular 2-dimensional array of bits. For example:
  //
  // row
  //  3|      
  //  2|      1
  //  1|    0 1
  //  0|  0 0 0
  //    -------
  //    0 1 2 3 column
  //
  // A bit is addressed by a pair (column, row) where the column goes first.
  // All valid address pairs have 0 <= row < column < columnCount().
  // Columns can be added dynamically.
  class BitTriangle {
  public:
	// Returns how many columns the triangle has
	std::size_t columnCount() const {return mColumns.size();}

	// Returns true if there are no columns in the triangle
	bool empty() const {return mColumns.empty();}

	// Adds a new column of the triangle. This increases columnCount() by
	// one, and the index of the new column is the previous value of
	// columnCount(). The new bits are all set to false initially.
	void addColumn() {
	  std::size_t const oldSize = mColumns.size();
	  mColumns.resize(oldSize + 1);
	  mColumns[oldSize].resize(oldSize);
	}

	// Returns the bit in the given column and row. As this is a triangle it
	// must be true that row < column.
	bool bit(std::size_t column, std::size_t row) const {
	  MATHIC_ASSERT(column < columnCount());
	  MATHIC_ASSERT(row < column);
	  return mColumns[column][row];
	}

	// As bit(), but uses max(x,y) as the column and min(x,y) as the
	// row.
	bool bitUnordered(std::size_t x, std::size_t y) const {
	  MATHIC_ASSERT(x < columnCount());
	  MATHIC_ASSERT(y < columnCount());
	  MATHIC_ASSERT(x != y);
	  if (x < y)
		std::swap(x, y);
	  return bit(x, y);
	}

	// Sets the bit in the given column and row. As this is a triangle
	// it must be true that column >= row.
	void setBit(std::size_t column, std::size_t row, bool value) {
	  MATHIC_ASSERT(column < columnCount());
	  MATHIC_ASSERT(row < column);
	  mColumns[column][row] = value;
	}

	// As setBit, but uses max(x,y) as the column and min(x,y) as the
	// row.
	void setBitUnordered(std::size_t x, std::size_t y, bool value) {
	  MATHIC_ASSERT(x < columnCount());
	  MATHIC_ASSERT(y < columnCount());
	  MATHIC_ASSERT(x != y);
	  if (x < y)
		std::swap(x, y);
	  setBit(x, y, value);
	}

	std::size_t getMemoryUse() const;

  private:
	// @todo: use one big array instead of an array of arrays.
	std::vector<std::vector<bool> > mColumns;
  };
}

#endif
