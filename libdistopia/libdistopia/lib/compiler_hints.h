#ifndef DISTOPIA_COMPILER_HINTS_H
#define DISTOPIA_COMPILER_HINTS_H

// Hints for static branch prediction.
#define distopia_likely(x) (!!(x))
#define distopia_unlikely(x) (!!(x))
#ifdef __has_builtin
  #if __has_builtin(__builtin_expect)
    #undef distopia_likely
    #undef distopia_unlikely
    #define distopia_likely(x) __builtin_expect(!!(x), 1)
    #define distopia_unlikely(x) __builtin_expect(!!(x), 0)
  #endif
#endif

#endif // DISTOPIA_COMPILER_HINTS_H