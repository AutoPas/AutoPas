/**
 * @file AlignedAllocator.h
 * @author tchipev
 * @date 07.02.2018
 */

#pragma once

#ifdef __SSE3__
#include <xmmintrin.h>
#endif

#include <cstdlib>
#include <limits>
#include <new>
#include <utility>

namespace autopas {

/**
 * Default size for a cache line.
 * @todo C++17: replace value by std::hardware_destructive_interference_size
 * not yet possible: as of 27.06.2019 this(P0154R1) is in neither libstc++ (gnu) nor libc++ (clang)
 */
constexpr unsigned int DEFAULT_CACHE_LINE_SIZE{64};

/**
 * AlignedAllocator class
 * @tparam T
 * @tparam Alignment
 */
template <class T, size_t Alignment = DEFAULT_CACHE_LINE_SIZE>
class AlignedAllocator {
 public:
  // needed for compatibility with stl::allocator
  /// value type
  using value_type = T;
  /// pointer type
  using pointer = T *;
  /// const pointer type
  using const_pointer = const T *;
  /// reference type
  using reference = T &;
  /// const reference type
  using const_reference = const T &;
  /// size type
  using size_type = size_t;

  /**
   * Equivalent allocator for other types
   * Class whose member other is an alias of allocator for type U.
   * (from cplusplus.com)
   * @tparam U
   */
  template <class U>
  struct rebind {
    /// other
    using other = AlignedAllocator<U, Alignment>;
  };

  /**
   * \brief Default empty constructor
   */
  AlignedAllocator() = default;

  /**
   * \brief Copy constructor
   */
  template <class U>
  AlignedAllocator(const AlignedAllocator<U, Alignment> &) {}

  /**
   * \brief Default destructor
   */
  ~AlignedAllocator() = default;

  /**
   * \brief Returns maximum possible value of n, with which we can call
   * allocate(n)
   * \return maximum size possible to allocate
   */
  size_t max_size() const noexcept { return (std::numeric_limits<size_t>::max() - size_t(Alignment)) / sizeof(T); }

  /**
   * \brief Allocate aligned memory for n objects of type T
   * \param n size to allocate
   * \return Pointer to the allocated memory
   */
  T *allocate(std::size_t n) {
    if (n <= max_size()) {
#if defined(_SX)
      T *ptr = static_cast<T *>(malloc(sizeof(T) * n));
#elif defined(__SSE3__) && !defined(__PGI)
      T *ptr = static_cast<T *>(_mm_malloc(sizeof(T) * n, Alignment));
#else
      // non-standard (obsolete):
      // T *ptr = static_cast<T *>(memalign(Alignment, sizeof(T) * n));

      // standard, but only in cstdlib 9 (gcc9) and libc++7 (clang7)
      /// @todo c++20: enable this!
      // T *ptr = static_cast<T *>(aligned_alloc(Alignment, sizeof(T) * n));

      // non-standard v2 (preferred over v1).
      // From gnu.org: The memalign function is obsolete and aligned_alloc or posix_memalign should be used instead.
      void *raw_ptr;
      posix_memalign(&raw_ptr, Alignment, sizeof(T) * n);
      T *ptr = static_cast<T *>(raw_ptr);
#endif
      if (ptr == nullptr) {
        throw std::bad_alloc();
      }
      return ptr;
    }
    throw std::bad_alloc();
  }

  /**
   * \brief Deallocate memory pointed to by ptr
   * \param ptr pointer to deallocate
   */
  void deallocate(T *ptr, std::size_t /*n*/) {
#if defined(__SSE3__) && !defined(__PGI)
    _mm_free(ptr);
#else
    free(ptr);
#endif
  }

  /**
   * \brief Construct object of type U at already allocated memory, pointed to
   * by p
   * \param p pointer to the object
   * \param args arguments for the construction
   */
  template <class U, class... Args>
  void construct(U *p, Args &&... args) {
    ::new ((void *)p) U(std::forward<Args>(args)...);
  }

  /**
   * \brief Destroy object pointed to by p, but does not deallocate the memory
   * \param p pointer to the object that should be destroyed
   */
  template <class U>
  void destroy(U *p) {
    p->~U();
  }
};

/**
 * equals operator
 * @tparam T Type of the first allocator
 * @tparam TAlignment alignment of the first allocator
 * @tparam U type of the second allocator
 * @tparam UAlignment alignment of the second allocator
 * @return true if alignment is equal
 */
template <typename T, size_t TAlignment, typename U, size_t UAlignment>
inline bool operator==(const AlignedAllocator<T, TAlignment> &, const AlignedAllocator<U, UAlignment> &) {
  return TAlignment == UAlignment;
}

/**
 * unequals operator
 * @tparam T Type of the first allocator
 * @tparam TAlignment alignment of the first allocator
 * @tparam U type of the second allocator
 * @tparam UAlignment alignment of the second allocator
 * @param a first allocator
 * @param b second allocator
 * @return true if unequal
 */
template <typename T, size_t TAlignment, typename U, size_t UAlignment>
inline bool operator!=(const AlignedAllocator<T, TAlignment> &a, const AlignedAllocator<U, UAlignment> &b) {
  return !(a == b);
}

}  // namespace autopas
