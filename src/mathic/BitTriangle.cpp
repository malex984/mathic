#include "stdinc.h"
#include "BitTriangle.h"

namespace mathic {
  std::size_t BitTriangle::getMemoryUse() const
  {
	std::size_t sum = mColumns.capacity() * sizeof(mColumns.front());
	const std::size_t stop = mColumns.size();
	for (std::size_t i = 0; i != stop; ++i) {
	  std::size_t const capacity = mColumns[i].capacity();
	  sum += (capacity + 7) / 8; // 8 bits per byte rounded up
	}
	return sum;
  }
}
