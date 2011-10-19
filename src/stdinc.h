#ifndef MATHIC_STDINC_GUARD
#define MATHIC_STDINC_GUARD

#include <cstddef>
#include <memory>

#if defined DEBUG || defined _DEBUG
#define MATHIC_DEBUG
#include <cassert>
#define MATHIC_DEBUG(X) ::assert(X);
#endif

#ifndef MATHIC_ASSERT
#define MATHIC_ASSERT(X)
#endif

namespace mathic {
  static const size_t BitsPerByte = 8;
  static const size_t MemoryAlignment = sizeof(void*);
}

#endif
