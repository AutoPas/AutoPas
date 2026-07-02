/* This is only the stb_arr part of stb.h and the required other functions
 
   stb.h - v2.37 - Sean's Tool Box -- public domain -- http://nothings.org/stb.h
          no warranty is offered or implied; use this code at your own risk

   This is a single header file with a bunch of useful utilities
   for getting stuff done in C/C++.

   Documentation: http://nothings.org/stb/stb_h.html
   Unit tests:    http://nothings.org/stb/stb.c

 ============================================================================
   You MUST

      #define STB_DEFINE

   in EXACTLY _one_ C or C++ file that includes this header, BEFORE the
   include, like this:

      #define STB_DEFINE
      #include "stb.h"

   All other files should just #include "stb.h" without the #define.
 ============================================================================

Version History

   2.36   various fixes
   2.35   fix clang-cl issues with swprintf
   2.34   fix warnings
   2.33   more fixes to random numbers
   2.32   stb_intcmprev, stb_uidict, fix random numbers on Linux
   2.31   stb_ucharcmp
   2.30   MinGW fix
   2.29   attempt to fix use of swprintf()
   2.28   various new functionality
   2.27   test _WIN32 not WIN32 in STB_THREADS
   2.26   various warning & bugfixes
   2.25   various warning & bugfixes
   2.24   various warning & bugfixes
   2.23   fix 2.22
   2.22   64-bit fixes from '!='; fix stb_sdict_copy() to have preferred name
   2.21   utf-8 decoder rejects "overlong" encodings; attempted 64-bit improvements
   2.20   fix to hash "copy" function--reported by someone with handle "!="
   2.19   ???
   2.18   stb_readdir_subdirs_mask
   2.17   stb_cfg_dir
   2.16   fix stb_bgio_, add stb_bgio_stat(); begin a streaming wrapper
   2.15   upgraded hash table template to allow:
            - aggregate keys (explicit comparison func for EMPTY and DEL keys)
            - "static" implementations (so they can be culled if unused)
   2.14   stb_mprintf
   2.13   reduce identifiable strings in STB_NO_STB_STRINGS
   2.12   fix STB_ONLY -- lots of uint32s, TRUE/FALSE things had crept in
   2.11   fix bug in stb_dirtree_get() which caused "c://path" sorts of stuff
   2.10   STB_F(), STB_I() inline constants (also KI,KU,KF,KD)
   2.09   stb_box_face_vertex_axis_side
   2.08   bugfix stb_trimwhite()
   2.07   colored printing in windows (why are we in 1985?)
   2.06   comparison functions are now functions-that-return-functions and
          accept a struct-offset as a parameter (not thread-safe)
   2.05   compile and pass tests under Linux (but no threads); thread cleanup
   2.04   stb_cubic_bezier_1d, smoothstep, avoid dependency on registry
   2.03   ?
   2.02   remove integrated documentation
   2.01   integrate various fixes; stb_force_uniprocessor
   2.00   revised stb_dupe to use multiple hashes
   1.99   stb_charcmp
   1.98   stb_arr_deleten, stb_arr_insertn
   1.97   fix stb_newell_normal()
   1.96   stb_hash_number()
   1.95   hack stb__rec_max; clean up recursion code to use new functions
   1.94   stb_dirtree; rename stb_extra to stb_ptrmap
   1.93   stb_sem_new() API cleanup (no blockflag-starts blocked; use 'extra')
   1.92   stb_threadqueue--multi reader/writer queue, fixed size or resizeable
   1.91   stb_bgio_* for reading disk asynchronously
   1.90   stb_mutex uses CRITICAL_REGION; new stb_sync primitive for thread
          joining; workqueue supports stb_sync instead of stb_semaphore
   1.89   support ';' in constant-string wildcards; stb_mutex wrapper (can
          implement with EnterCriticalRegion eventually)
   1.88   portable threading API (only for win32 so far); worker thread queue
   1.87   fix wildcard handling in stb_readdir_recursive
   1.86   support ';' in wildcards
   1.85   make stb_regex work with non-constant strings;
               beginnings of stb_introspect()
   1.84   (forgot to make notes)
   1.83   whoops, stb_keep_if_different wasn't deleting the temp file
   1.82   bring back stb_compress from stb_file.h for cmirror
   1.81   various bugfixes, STB_FASTMALLOC_INIT inits FASTMALLOC in release
   1.80   stb_readdir returns utf8; write own utf8-utf16 because lib was wrong
   1.79   stb_write
   1.78   calloc() support for malloc wrapper, STB_FASTMALLOC
   1.77   STB_FASTMALLOC
   1.76   STB_STUA - Lua-like language; (stb_image, stb_csample, stb_bilinear)
   1.75   alloc/free array of blocks; stb_hheap bug; a few stb_ps_ funcs;
          hash*getkey, hash*copy; stb_bitset; stb_strnicmp; bugfix stb_bst
   1.74   stb_replaceinplace; use stdlib C function to convert utf8 to UTF-16
   1.73   fix performance bug & leak in stb_ischar (C++ port lost a 'static')
   1.72   remove stb_block, stb_block_manager, stb_decompress (to stb_file.h)
   1.71   stb_trimwhite, stb_tokens_nested, etc.
   1.70   back out 1.69 because it might problemize mixed builds; stb_filec()
   1.69   (stb_file returns 'char *' in C++)
   1.68   add a special 'tree root' data type for stb_bst; stb_arr_end
   1.67   full C++ port. (stb_block_manager)
   1.66   stb_newell_normal
   1.65   stb_lex_item_wild -- allow wildcard items which MUST match entirely
   1.64   stb_data
   1.63   stb_log_name
   1.62   stb_define_sort; C++ cleanup
   1.61   stb_hash_fast -- Paul Hsieh's hash function (beats Bob Jenkins'?)
   1.60   stb_delete_directory_recursive
   1.59   stb_readdir_recursive
   1.58   stb_bst variant with parent pointer for O(1) iteration, not O(log N)
   1.57   replace LCG random with Mersenne Twister (found a public domain one)
   1.56   stb_perfect_hash, stb_ischar, stb_regex
   1.55   new stb_bst API allows multiple BSTs per node (e.g. secondary keys)
   1.54   bugfix: stb_define_hash, stb_wildmatch, regexp
   1.53   stb_define_hash; recoded stb_extra, stb_sdict use it
   1.52   stb_rand_define, stb_bst, stb_reverse
   1.51   fix 'stb_arr_setlen(NULL, 0)'
   1.50   stb_wordwrap
   1.49   minor improvements to enable the scripting language
   1.48   better approach for stb_arr using stb_malloc; more invasive, clearer
   1.47   stb_lex (lexes stb.h at 1.5ML/s on 3Ghz P4; 60/70% of optimal/flex)
   1.46   stb_wrapper_*, STB_MALLOC_WRAPPER
   1.45   lightly tested DFA acceleration of regexp searching
   1.44   wildcard matching & searching; regexp matching & searching
   1.43   stb_temp
   1.42   allow stb_arr to use stb_malloc/realloc; note this is global
   1.41   make it compile in C++; (disable stb_arr in C++)
   1.40   stb_dupe tweak; stb_swap; stb_substr
   1.39   stb_dupe; improve stb_file_max to be less stupid
   1.38   stb_sha1_file: generate sha1 for file, even > 4GB
   1.37   stb_file_max; partial support for utf8 filenames in Windows
   1.36   remove STB__NO_PREFIX - poor interaction with IDE, not worth it
          streamline stb_arr to make it separately publishable
   1.35   bugfixes for stb_sdict, stb_malloc(0), stristr
   1.34   (streaming interfaces for stb_compress)
   1.33   stb_alloc; bug in stb_getopt; remove stb_overflow
   1.32   (stb_compress returns, smaller&faster; encode window & 64-bit len)
   1.31   stb_prefix_count
   1.30   (STB__NO_PREFIX - remove stb_ prefixes for personal projects)
   1.29   stb_fput_varlen64, etc.
   1.28   stb_sha1
   1.27   ?
   1.26   stb_extra
   1.25   ?
   1.24   stb_copyfile
   1.23   stb_readdir
   1.22   ?
   1.21   ?
   1.20   ?
   1.19   ?
   1.18   ?
   1.17   ?
   1.16   ?
   1.15   stb_fixpath, stb_splitpath, stb_strchr2
   1.14   stb_arr
   1.13   ?stb, stb_log, stb_fatal
   1.12   ?stb_hash2
   1.11   miniML
   1.10   stb_crc32, stb_adler32
   1.09   stb_sdict
   1.08   stb_bitreverse, stb_ispow2, stb_big32
          stb_fopen, stb_fput_varlen, stb_fput_ranged
          stb_fcmp, stb_feq
   1.07   (stb_encompress)
   1.06   stb_compress
   1.05   stb_tokens, (stb_hheap)
   1.04   stb_rand
   1.03   ?(s-strings)
   1.02   ?stb_filelen, stb_tokens
   1.01   stb_tolower
   1.00   stb_hash, stb_intcmp
          stb_file, stb_stringfile, stb_fgets
          stb_prefix, stb_strlower, stb_strtok
          stb_image
          (stb_array), (stb_arena)

Parenthesized items have since been removed.

LICENSE

 See end of file for license information.

CREDITS

 Written by Sean Barrett.

 Fixes:
  Philipp Wiesemann
  Robert Nix
  r-lyeh
  blackpawn
  github:Mojofreem
  Ryan Whitworth
  Vincent Isambart
  Mike Sartain
  Eugene Opalev
  Tim Sjostrand
  github:infatum
  Dave Butler (Croepha)
  Ethan Lee (flibitijibibo)
  Brian Collins
  Kyle Langley
*/
#include <stdarg.h>

#ifndef STB__INCLUDE_STB_H
#define STB__INCLUDE_STB_H

#include <stdlib.h>     // stdlib could have min/max

#ifdef STB_DEFINE
   #include <assert.h>
   #include <stdarg.h>
   #include <stddef.h>
   #include <ctype.h>
   #include <math.h>
   #include <string.h>
#endif

#ifdef __cplusplus
   #define STB_EXTERN   extern "C"
#else
   #define STB_EXTERN   extern
#endif

// if we're STB_ONLY, can't rely on uint32 or even uint, so all the
// variables we'll use herein need typenames prefixed with 'stb':
typedef unsigned char stb_uchar;
typedef unsigned char stb_uint8;
typedef unsigned int  stb_uint;
typedef unsigned short stb_uint16;
typedef          short stb_int16;
typedef   signed char  stb_int8;
#if defined(STB_USE_LONG_FOR_32_BIT_INT) || defined(STB_LONG32)
  typedef unsigned long  stb_uint32;
  typedef          long  stb_int32;
#else
  typedef unsigned int   stb_uint32;
  typedef          int   stb_int32;
#endif
typedef char stb__testsize2_16[sizeof(stb_uint16)==2 ? 1 : -1];
typedef char stb__testsize2_32[sizeof(stb_uint32)==4 ? 1 : -1];

#ifdef _MSC_VER
  typedef unsigned __int64 stb_uint64;
  typedef          __int64 stb_int64;
  #define STB_IMM_UINT64(literalui64) (literalui64##ui64)
  #define STB_IMM_INT64(literali64) (literali64##i64)
#else
  // ??
  typedef unsigned long long stb_uint64;
  typedef          long long stb_int64;
  #define STB_IMM_UINT64(literalui64) (literalui64##ULL)
  #define STB_IMM_INT64(literali64) (literali64##LL)
#endif
typedef char stb__testsize2_64[sizeof(stb_uint64)==8 ? 1 : -1];


// add platform-specific ways of checking for sizeof(char*) == 8,
// and make those define STB_PTR64
#if defined(_WIN64) || defined(__x86_64__) || defined(__ia64__) || defined(__LP64__)
  #define STB_PTR64
#endif

#ifdef STB_PTR64
typedef char stb__testsize2_ptr[sizeof(char *) == 8];
typedef stb_uint64 stb_uinta;
typedef stb_int64  stb_inta;
#else
typedef char stb__testsize2_ptr[sizeof(char *) == 4];
typedef stb_uint32 stb_uinta;
typedef stb_int32  stb_inta;
#endif
typedef char stb__testsize2_uinta[sizeof(stb_uinta)==sizeof(char*) ? 1 : -1];

// if so, we should define an int type that is the pointer size. until then,
// we'll have to make do with this (which is not the same at all!)

typedef union
{
   unsigned int i;
   void * p;
} stb_uintptr;

#define stb_min(a,b)   ((a) < (b) ? (a) : (b))
#define stb_max(a,b)   ((a) > (b) ? (a) : (b))

//////////////////////////////////////////////////////////////////////////////
//
//                         bit operations
//

#define stb_big32(c)    (((c)[0]<<24) + (c)[1]*65536 + (c)[2]*256 + (c)[3])
#define stb_little32(c) (((c)[3]<<24) + (c)[2]*65536 + (c)[1]*256 + (c)[0])
#define stb_big16(c)    ((c)[0]*256 + (c)[1])
#define stb_little16(c) ((c)[1]*256 + (c)[0])

STB_EXTERN          int stb_bitcount(unsigned int a);
STB_EXTERN unsigned int stb_bitreverse8(unsigned char n);
STB_EXTERN unsigned int stb_bitreverse(unsigned int n);

STB_EXTERN          int stb_is_pow2(size_t);
STB_EXTERN          int stb_log2_ceil(size_t);
STB_EXTERN          int stb_log2_floor(size_t);

STB_EXTERN          int stb_lowbit8(unsigned int n);
STB_EXTERN          int stb_highbit8(unsigned int n);

#ifdef STB_DEFINE
int stb_bitcount(unsigned int a)
{
   a = (a & 0x55555555) + ((a >>  1) & 0x55555555); // max 2
   a = (a & 0x33333333) + ((a >>  2) & 0x33333333); // max 4
   a = (a + (a >> 4)) & 0x0f0f0f0f; // max 8 per 4, now 8 bits
   a = (a + (a >> 8)); // max 16 per 8 bits
   a = (a + (a >> 16)); // max 32 per 8 bits
   return a & 0xff;
}

unsigned int stb_bitreverse8(unsigned char n)
{
   n = ((n & 0xAA) >> 1) + ((n & 0x55) << 1);
   n = ((n & 0xCC) >> 2) + ((n & 0x33) << 2);
   return (unsigned char) ((n >> 4) + (n << 4));
}

unsigned int stb_bitreverse(unsigned int n)
{
  n = ((n & 0xAAAAAAAA) >>  1) | ((n & 0x55555555) << 1);
  n = ((n & 0xCCCCCCCC) >>  2) | ((n & 0x33333333) << 2);
  n = ((n & 0xF0F0F0F0) >>  4) | ((n & 0x0F0F0F0F) << 4);
  n = ((n & 0xFF00FF00) >>  8) | ((n & 0x00FF00FF) << 8);
  return (n >> 16) | (n << 16);
}

int stb_is_pow2(size_t n)
{
   return (n & (n-1)) == 0;
}

// tricky use of 4-bit table to identify 5 bit positions (note the '-1')
// 3-bit table would require another tree level; 5-bit table wouldn't save one
#if defined(_WIN32) && !defined(__MINGW32__)
#pragma warning(push)
#pragma warning(disable: 4035)  // disable warning about no return value
int stb_log2_floor(size_t n)
{
   #if _MSC_VER > 1700
   unsigned long i;
   #ifdef STB_PTR64
   _BitScanReverse64(&i, n);
   #else
   _BitScanReverse(&i, n);
   #endif
   return i != 0 ? i : -1;
   #else
   __asm {
      bsr eax,n
      jnz done
      mov eax,-1
   }
   done:;
   #endif
}
#pragma warning(pop)
#else
int stb_log2_floor(size_t n)
{
   static signed char log2_4[16] = { -1,0,1,1,2,2,2,2,3,3,3,3,3,3,3,3 };

#ifdef STB_PTR64
   if (n >= ((size_t) 1u << 32))
        return stb_log2_floor(n >> 32);
#endif

   // 2 compares if n < 16, 3 compares otherwise
   if (n < (1U << 14))
        if (n < (1U <<  4))        return     0 + log2_4[n      ];
        else if (n < (1U <<  9))      return  5 + log2_4[n >>  5];
             else                     return 10 + log2_4[n >> 10];
   else if (n < (1U << 24))
             if (n < (1U << 19))      return 15 + log2_4[n >> 15];
             else                     return 20 + log2_4[n >> 20];
        else if (n < (1U << 29))      return 25 + log2_4[n >> 25];
             else                     return 30 + log2_4[n >> 30];
}
#endif

// define ceil from floor
int stb_log2_ceil(size_t n)
{
   if (stb_is_pow2(n))  return     stb_log2_floor(n);
   else                 return 1 + stb_log2_floor(n);
}

int stb_highbit8(unsigned int n)
{
   return stb_log2_ceil(n&255);
}

int stb_lowbit8(unsigned int n)
{
   static signed char lowbit4[16] = { -1,0,1,0, 2,0,1,0, 3,0,1,0, 2,0,1,0 };
   int k = lowbit4[n & 15];
   if (k >= 0) return k;
   k = lowbit4[(n >> 4) & 15];
   if (k >= 0) return k+4;
   return k;
}
#endif


//////////////////////////////////////////////////////////////////////////////
//
//                   stb_alloc - hierarchical allocator
//
//                                     inspired by http://swapped.cc/halloc
//
//
// When you alloc a given block through stb_alloc, you have these choices:
//
//       1. does it have a parent?
//       2. can it have children?
//       3. can it be freed directly?
//       4. is it transferrable?
//       5. what is its alignment?
//
// Here are interesting combinations of those:
//
//                              children   free    transfer     alignment
//  arena                          Y         Y         N           n/a
//  no-overhead, chunked           N         N         N         normal
//  string pool alloc              N         N         N            1
//  parent-ptr, chunked            Y         N         N         normal
//  low-overhead, unchunked        N         Y         Y         normal
//  general purpose alloc          Y         Y         Y         normal
//
// Unchunked allocations will probably return 16-aligned pointers. If
// we 16-align the results, we have room for 4 pointers. For smaller
// allocations that allow finer alignment, we can reduce the pointers.
//
// The strategy is that given a pointer, assuming it has a header (only
// the no-overhead allocations have no header), we can determine the
// type of the header fields, and the number of them, by stepping backwards
// through memory and looking at the tags in the bottom bits.
//
// Implementation strategy:
//     chunked allocations come from the middle of chunks, and can't
//     be freed. thefore they do not need to be on a sibling chain.
//     they may need child pointers if they have children.
//
// chunked, with-children
//     void *parent;
//
// unchunked, no-children -- reduced storage
//     void *next_sibling;
//     void *prev_sibling_nextp;
//
// unchunked, general
//     void *first_child;
//     void *next_sibling;
//     void *prev_sibling_nextp;
//     void *chunks;
//
// so, if we code each of these fields with different bit patterns
// (actually same one for next/prev/child), then we can identify which
// each one is from the last field.

STB_EXTERN void  stb_free(void *p);
STB_EXTERN void *stb_malloc_global(size_t size);
STB_EXTERN void *stb_malloc(void *context, size_t size);
STB_EXTERN void *stb_malloc_nofree(void *context, size_t size);
STB_EXTERN void *stb_malloc_leaf(void *context, size_t size);
STB_EXTERN void *stb_malloc_raw(void *context, size_t size);
STB_EXTERN void *stb_realloc(void *ptr, size_t newsize);

STB_EXTERN void stb_reassign(void *new_context, void *ptr);
STB_EXTERN void stb_malloc_validate(void *p, void *parent);

extern int stb_alloc_chunk_size ;
extern int stb_alloc_count_free ;
extern int stb_alloc_count_alloc;
extern int stb_alloc_alignment  ;

#ifdef STB_DEFINE

int stb_alloc_chunk_size  = 65536;
int stb_alloc_count_free  = 0;
int stb_alloc_count_alloc = 0;
int stb_alloc_alignment   = -16;

typedef struct stb__chunk
{
   struct stb__chunk *next;
   int                data_left;
   int                alloc;
} stb__chunk;

typedef struct
{
   void *  next;
   void ** prevn;
} stb__nochildren;

typedef struct
{
   void ** prevn;
   void *  child;
   void *  next;
   stb__chunk *chunks;
} stb__alloc;

typedef struct
{
   stb__alloc *parent;
} stb__chunked;

#define STB__PARENT          1
#define STB__CHUNKS          2

typedef enum
{
   STB__nochildren = 0,
   STB__chunked    = STB__PARENT,
   STB__alloc      = STB__CHUNKS,

   STB__chunk_raw  = 4,
} stb__alloc_type;

// these functions set the bottom bits of a pointer efficiently
#define STB__DECODE(x,v)  ((void *) ((char *) (x) - (v)))
#define STB__ENCODE(x,v)  ((void *) ((char *) (x) + (v)))

#define stb__parent(z)       (stb__alloc *) STB__DECODE((z)->parent, STB__PARENT)
#define stb__chunks(z)       (stb__chunk *) STB__DECODE((z)->chunks, STB__CHUNKS)

#define stb__setparent(z,p)  (z)->parent = (stb__alloc *) STB__ENCODE((p), STB__PARENT)
#define stb__setchunks(z,c)  (z)->chunks = (stb__chunk *) STB__ENCODE((c), STB__CHUNKS)

static stb__alloc stb__alloc_global =
{
   NULL,
   NULL,
   NULL,
   (stb__chunk *) STB__ENCODE(NULL, STB__CHUNKS)
};

static stb__alloc_type stb__identify(void *p)
{
   void **q = (void **) p;
   return (stb__alloc_type) ((stb_uinta) q[-1] & 3);
}

static void *** stb__prevn(void *p)
{
   if (stb__identify(p) == STB__alloc) {
      stb__alloc      *s = (stb__alloc *) p - 1;
      return &s->prevn;
   } else {
      stb__nochildren *s = (stb__nochildren *) p - 1;
      return &s->prevn;
   }
}

void stb_free(void *p)
{
   if (p == NULL) return;

   // count frees so that unit tests can see what's happening
   ++stb_alloc_count_free;

   switch(stb__identify(p)) {
      case STB__chunked:
         // freeing a chunked-block with children does nothing;
         // they only get freed when the parent does
         // surely this is wrong, and it should free them immediately?
         // otherwise how are they getting put on the right chain?
         return;
      case STB__nochildren: {
         stb__nochildren *s = (stb__nochildren *) p - 1;
         // unlink from sibling chain
         *(s->prevn) = s->next;
         if (s->next)
            *stb__prevn(s->next) = s->prevn;
         free(s);
         return;
      }
      case STB__alloc: {
         stb__alloc *s = (stb__alloc *) p - 1;
         stb__chunk *c, *n;
         void *q;

         // unlink from sibling chain, if any
         *(s->prevn) = s->next;
         if (s->next)
            *stb__prevn(s->next) = s->prevn;

         // first free chunks
         c = (stb__chunk *) stb__chunks(s);
         while (c != NULL) {
            n = c->next;
            stb_alloc_count_free += c->alloc;
            free(c);
            c = n;
         }

         // validating
         stb__setchunks(s,NULL);
         s->prevn = NULL;
         s->next = NULL;

         // now free children
         while ((q = s->child) != NULL) {
            stb_free(q);
         }

         // now free self
         free(s);
         return;
      }
      default:
         assert(0); /* NOTREACHED */
   }
}

void stb_malloc_validate(void *p, void *parent)
{
   if (p == NULL) return;

   switch(stb__identify(p)) {
      case STB__chunked:
         return;
      case STB__nochildren: {
         stb__nochildren *n = (stb__nochildren *) p - 1;
         if (n->prevn)
            assert(*n->prevn == p);
         if (n->next) {
            assert(*stb__prevn(n->next) == &n->next);
            stb_malloc_validate(n, parent);
         }
         return;
      }
      case STB__alloc: {
         stb__alloc *s = (stb__alloc *) p - 1;

         if (s->prevn)
            assert(*s->prevn == p);

         if (s->child) {
            assert(*stb__prevn(s->child) == &s->child);
            stb_malloc_validate(s->child, p);
         }

         if (s->next) {
            assert(*stb__prevn(s->next) == &s->next);
            stb_malloc_validate(s->next, parent);
         }
         return;
      }
      default:
         assert(0); /* NOTREACHED */
   }
}

static void * stb__try_chunk(stb__chunk *c, int size, int align, int pre_align)
{
   char *memblock = (char *) (c+1), *q;
   stb_inta iq;
   int start_offset;

   // we going to allocate at the end of the chunk, not the start. confusing,
   // but it means we don't need both a 'limit' and a 'cur', just a 'cur'.
   // the block ends at: p + c->data_left
   //   then we move back by size
   start_offset = c->data_left - size;

   // now we need to check the alignment of that
   q = memblock + start_offset;
   iq = (stb_inta) q;
   assert(sizeof(q) == sizeof(iq));

   // suppose align = 2
   // then we need to retreat iq far enough that (iq & (2-1)) == 0
   // to get (iq & (align-1)) = 0 requires subtracting (iq & (align-1))

   start_offset -= iq & (align-1);
   assert(((stb_uinta) (memblock+start_offset) & (align-1)) == 0);

   // now, if that + pre_align works, go for it!
   start_offset -= pre_align;

   if (start_offset >= 0) {
      c->data_left = start_offset;
      return memblock + start_offset;
   }

   return NULL;
}

static void stb__sort_chunks(stb__alloc *src)
{
   // of the first two chunks, put the chunk with more data left in it first
   stb__chunk *c = stb__chunks(src), *d;
   if (c == NULL) return;
   d = c->next;
   if (d == NULL) return;
   if (c->data_left > d->data_left) return;

   c->next = d->next;
   d->next = c;
   stb__setchunks(src, d);
}

static void * stb__alloc_chunk(stb__alloc *src, int size, int align, int pre_align)
{
   void *p;
   stb__chunk *c = stb__chunks(src);

   if (c && size <= stb_alloc_chunk_size) {

      p = stb__try_chunk(c, size, align, pre_align);
      if (p) { ++c->alloc; return p; }

      // try a second chunk to reduce wastage
      if (c->next) {
         p = stb__try_chunk(c->next, size, align, pre_align);
         if (p) { ++c->alloc; return p; }

         // put the bigger chunk first, since the second will get buried
         // the upshot of this is that, until it gets allocated from, chunk #2
         // is always the largest remaining chunk. (could formalize
         // this with a heap!)
         stb__sort_chunks(src);
         c = stb__chunks(src);
      }
   }

   // allocate a new chunk
   {
      stb__chunk *n;

      int chunk_size = stb_alloc_chunk_size;
      // we're going to allocate a new chunk to put this in
      if (size > chunk_size)
         chunk_size = size;

      assert(sizeof(*n) + pre_align <= 16);

      // loop trying to allocate a large enough chunk
      // the loop is because the alignment may cause problems if it's big...
      // and we don't know what our chunk alignment is going to be
      while (1) {
         n = (stb__chunk *) malloc(16 + chunk_size);
         if (n == NULL) return NULL;

         n->data_left = chunk_size - sizeof(*n);

         p = stb__try_chunk(n, size, align, pre_align);
         if (p != NULL) {
            n->next = c;
            stb__setchunks(src, n);

            // if we just used up the whole block immediately,
            // move the following chunk up
            n->alloc = 1;
            if (size == chunk_size)
               stb__sort_chunks(src);

            return p;
         }

         free(n);
         chunk_size += 16+align;
      }
   }
}

static stb__alloc * stb__get_context(void *context)
{
   if (context == NULL) {
      return &stb__alloc_global;
   } else {
      int u = stb__identify(context);
      // if context is chunked, grab parent
      if (u == STB__chunked) {
         stb__chunked *s = (stb__chunked *) context - 1;
         return stb__parent(s);
      } else {
         return (stb__alloc *) context - 1;
      }
   }
}

static void stb__insert_alloc(stb__alloc *src, stb__alloc *s)
{
   s->prevn = &src->child;
   s->next  = src->child;
   src->child = s+1;
   if (s->next)
      *stb__prevn(s->next) = &s->next;
}

static void stb__insert_nochild(stb__alloc *src, stb__nochildren *s)
{
   s->prevn = &src->child;
   s->next  = src->child;
   src->child = s+1;
   if (s->next)
      *stb__prevn(s->next) = &s->next;
}

static void * malloc_base(void *context, size_t size, stb__alloc_type t, int align)
{
   void *p;

   stb__alloc *src = stb__get_context(context);

   if (align <= 0) {
      // compute worst-case C packed alignment
      // e.g. a 24-byte struct is 8-aligned
      int align_proposed = 1 << stb_lowbit8((unsigned int) size);

      if (align_proposed < 0)
         align_proposed = 4;

      if (align_proposed == 0) {
         if (size == 0)
            align_proposed = 1;
         else
            align_proposed = 256;
      }

      // a negative alignment means 'don't align any larger
      // than this'; so -16 means we align 1,2,4,8, or 16

      if (align < 0) {
         if (align_proposed > -align)
            align_proposed = -align;
      }

      align = align_proposed;
   }

   assert(stb_is_pow2(align));

   // don't cause misalignment when allocating nochildren
   if (t == STB__nochildren && align > 8)
      t = STB__alloc;

   switch (t) {
      case STB__alloc: {
         stb__alloc *s = (stb__alloc *) malloc(size + sizeof(*s));
         if (s == NULL) return NULL;
         p = s+1;
         s->child = NULL;
         stb__insert_alloc(src, s);

         stb__setchunks(s,NULL);
         break;
      }

      case STB__nochildren: {
         stb__nochildren *s = (stb__nochildren *) malloc(size + sizeof(*s));
         if (s == NULL) return NULL;
         p = s+1;
         stb__insert_nochild(src, s);
         break;
      }

      case STB__chunk_raw: {
         p = stb__alloc_chunk(src, (int) size, align, 0);
         if (p == NULL) return NULL;
         break;
      }

      case STB__chunked: {
         stb__chunked *s;
         if (align < sizeof(stb_uintptr)) align = sizeof(stb_uintptr);
         s = (stb__chunked *) stb__alloc_chunk(src, (int) size, align, sizeof(*s));
         if (s == NULL) return NULL;
         stb__setparent(s, src);
         p = s+1;
         break;
      }

      default: p = NULL; assert(0); /* NOTREACHED */
   }

   ++stb_alloc_count_alloc;
   return p;
}

void *stb_malloc_global(size_t size)
{
   return malloc_base(NULL, size, STB__alloc, stb_alloc_alignment);
}

void *stb_malloc(void *context, size_t size)
{
   return malloc_base(context, size, STB__alloc, stb_alloc_alignment);
}

void *stb_malloc_nofree(void *context, size_t size)
{
   return malloc_base(context, size, STB__chunked, stb_alloc_alignment);
}

void *stb_malloc_leaf(void *context, size_t size)
{
   return malloc_base(context, size, STB__nochildren, stb_alloc_alignment);
}

void *stb_malloc_raw(void *context, size_t size)
{
   return malloc_base(context, size, STB__chunk_raw, stb_alloc_alignment);
}

char *stb_malloc_string(void *context, size_t size)
{
   return (char *) malloc_base(context, size, STB__chunk_raw, 1);
}

void *stb_realloc(void *ptr, size_t newsize)
{
   stb__alloc_type t;

   if (ptr == NULL) return stb_malloc(NULL, newsize);
   if (newsize == 0) { stb_free(ptr); return NULL; }

   t = stb__identify(ptr);
   assert(t == STB__alloc || t == STB__nochildren);

   if (t == STB__alloc) {
      stb__alloc *s = (stb__alloc *) ptr - 1;

      s = (stb__alloc *) realloc(s, newsize + sizeof(*s));
      if (s == NULL) return NULL;

      ptr = s+1;

      // update pointers
      (*s->prevn) = ptr;
      if (s->next)
         *stb__prevn(s->next) = &s->next;

      if (s->child)
         *stb__prevn(s->child) = &s->child;

      return ptr;
   } else {
      stb__nochildren *s = (stb__nochildren *) ptr - 1;

      s = (stb__nochildren *) realloc(ptr, newsize + sizeof(s));
      if (s == NULL) return NULL;

      // update pointers
      (*s->prevn) = s+1;
      if (s->next)
         *stb__prevn(s->next) = &s->next;

      return s+1;
   }
}

void *stb_realloc_c(void *context, void *ptr, size_t newsize)
{
   if (ptr == NULL) return stb_malloc(context, newsize);
   if (newsize == 0) { stb_free(ptr); return NULL; }
   // @TODO: verify you haven't changed contexts
   return stb_realloc(ptr, newsize);
}

void stb_reassign(void *new_context, void *ptr)
{
   stb__alloc *src = stb__get_context(new_context);

   stb__alloc_type t = stb__identify(ptr);
   assert(t == STB__alloc || t == STB__nochildren);

   if (t == STB__alloc) {
      stb__alloc *s = (stb__alloc *) ptr - 1;

      // unlink from old
      *(s->prevn) = s->next;
      if (s->next)
         *stb__prevn(s->next) = s->prevn;

      stb__insert_alloc(src, s);
   } else {
      stb__nochildren *s = (stb__nochildren *) ptr - 1;

      // unlink from old
      *(s->prevn) = s->next;
      if (s->next)
         *stb__prevn(s->next) = s->prevn;

      stb__insert_nochild(src, s);
   }
}

#endif

//////////////////////////////////////////////////////////////////////////////
//
//                                stb_arr
//
//  An stb_arr is directly useable as a pointer (use the actual type in your
//  definition), but when it resizes, it returns a new pointer and you can't
//  use the old one, so you have to be careful to copy-in-out as necessary.
//
//  Use a NULL pointer as a 0-length array.
//
//     float *my_array = NULL, *temp;
//
//     // add elements on the end one at a time
//     stb_arr_push(my_array, 0.0f);
//     stb_arr_push(my_array, 1.0f);
//     stb_arr_push(my_array, 2.0f);
//
//     assert(my_array[1] == 2.0f);
//
//     // add an uninitialized element at the end, then assign it
//     *stb_arr_add(my_array) = 3.0f;
//
//     // add three uninitialized elements at the end
//     temp = stb_arr_addn(my_array,3);
//     temp[0] = 4.0f;
//     temp[1] = 5.0f;
//     temp[2] = 6.0f;
//
//     assert(my_array[5] == 5.0f);
//
//     // remove the last one
//     stb_arr_pop(my_array);
//
//     assert(stb_arr_len(my_array) == 6);


#ifdef STB_MALLOC_WRAPPER
  #define STB__PARAMS    , char *file, int line
  #define STB__ARGS      ,       file,     line
#else
  #define STB__PARAMS
  #define STB__ARGS
#endif

// calling this function allocates an empty stb_arr attached to p
// (whereas NULL isn't attached to anything)
STB_EXTERN void stb_arr_malloc(void **target, void *context);

// call this function with a non-NULL value to have all successive
// stbs that are created be attached to the associated parent. Note
// that once a given stb_arr is non-empty, it stays attached to its
// current parent, even if you call this function again.
// it turns the previous value, so you can restore it
STB_EXTERN void* stb_arr_malloc_parent(void *p);

// simple functions written on top of other functions
#define stb_arr_empty(a)       (  stb_arr_len(a) == 0 )
#define stb_arr_add(a)         (  stb_arr_addn((a),1) )
#define stb_arr_push(a,v)      ( *stb_arr_add(a)=(v)  )

typedef struct
{
   int len, limit;
   int stb_malloc;
   unsigned int signature;
} stb__arr;

#define stb_arr_signature      0x51bada7b  // ends with 0123 in decimal

// access the header block stored before the data
#define stb_arrhead(a)         /*lint --e(826)*/ (((stb__arr *) (a)) - 1)
#define stb_arrhead2(a)        /*lint --e(826)*/ (((stb__arr *) (a)) - 1)

#ifdef STB_DEBUG
#define stb_arr_check(a)       assert(!a || stb_arrhead(a)->signature == stb_arr_signature)
#define stb_arr_check2(a)      assert(!a || stb_arrhead2(a)->signature == stb_arr_signature)
#else
#define stb_arr_check(a)       ((void) 0)
#define stb_arr_check2(a)      ((void) 0)
#endif

// ARRAY LENGTH

// get the array length; special case if pointer is NULL
#define stb_arr_len(a)         (a ? stb_arrhead(a)->len : 0)
#define stb_arr_len2(a)        ((stb__arr *) (a) ? stb_arrhead2(a)->len : 0)
#define stb_arr_lastn(a)       (stb_arr_len(a)-1)

// check whether a given index is valid -- tests 0 <= i < stb_arr_len(a)
#define stb_arr_valid(a,i)     (a ? (int) (i) < stb_arrhead(a)->len : 0)

// change the array length so is is exactly N entries long, creating
// uninitialized entries as needed
#define stb_arr_setlen(a,n)  \
            (stb__arr_setlen((void **) &(a), sizeof(a[0]), (n)))

// change the array length so that N is a valid index (that is, so
// it is at least N entries long), creating uninitialized entries as needed
#define stb_arr_makevalid(a,n)  \
            (stb_arr_len(a) < (n)+1 ? stb_arr_setlen(a,(n)+1),(a) : (a))

// remove the last element of the array, returning it
#define stb_arr_pop(a)         ((stb_arr_check(a), (a))[--stb_arrhead(a)->len])

// access the last element in the array
#define stb_arr_last(a)        ((stb_arr_check(a), (a))[stb_arr_len(a)-1])

// is iterator at end of list?
#define stb_arr_end(a,i)       ((i) >= &(a)[stb_arr_len(a)])

// (internal) change the allocated length of the array
#define stb_arr__grow(a,n)     (stb_arr_check(a), stb_arrhead(a)->len += (n))

// add N new uninitialized elements to the end of the array
#define stb_arr__addn(a,n)     /*lint --e(826)*/ \
                               ((stb_arr_len(a)+(n) > stb_arrcurmax(a))      \
                                 ? (stb__arr_addlen((void **) &(a),sizeof(*a),(n)),0) \
                                 : ((stb_arr__grow(a,n), 0)))

// add N new uninitialized elements to the end of the array, and return
// a pointer to the first new one
#define stb_arr_addn(a,n)      (stb_arr__addn((a),n),(a)+stb_arr_len(a)-(n))

// add N new uninitialized elements starting at index 'i'
#define stb_arr_insertn(a,i,n) (stb__arr_insertn((void **) &(a), sizeof(*a), (i), (n)))

// insert an element at i
#define stb_arr_insert(a,i,v)  (stb__arr_insertn((void **) &(a), sizeof(*a), (i), (1)), ((a)[i] = v))

// delete N elements from the middle starting at index 'i'
#define stb_arr_deleten(a,i,n) (stb__arr_deleten((void **) &(a), sizeof(*a), (i), (n)))

// delete the i'th element
#define stb_arr_delete(a,i)   stb_arr_deleten(a,i,1)

// delete the i'th element, swapping down from the end
#define stb_arr_fastdelete(a,i)  \
   (stb_swap(&a[i], &a[stb_arrhead(a)->len-1], sizeof(*a)), stb_arr_pop(a))


// ARRAY STORAGE

// get the array maximum storage; special case if NULL
#define stb_arrcurmax(a)       (a ? stb_arrhead(a)->limit : 0)
#define stb_arrcurmax2(a)      (a ? stb_arrhead2(a)->limit : 0)

// set the maxlength of the array to n in anticipation of further growth
#define stb_arr_setsize(a,n)   (stb_arr_check(a), stb__arr_setsize((void **) &(a),sizeof((a)[0]),n))

// make sure maxlength is large enough for at least N new allocations
#define stb_arr_atleast(a,n)   (stb_arr_len(a)+(n) > stb_arrcurmax(a)      \
                                 ? stb_arr_setsize((a), (n)) : 0)

// make a copy of a given array (copies contents via 'memcpy'!)
#define stb_arr_copy(a)        stb__arr_copy(a, sizeof((a)[0]))

// compute the storage needed to store all the elements of the array
#define stb_arr_storage(a)     (stb_arr_len(a) * sizeof((a)[0]))

#define stb_arr_for(v,arr)     for((v)=(arr); (v) < (arr)+stb_arr_len(arr); ++(v))

// IMPLEMENTATION

STB_EXTERN void stb_arr_free_(void **p);
STB_EXTERN void *stb__arr_copy_(void *p, int elem_size);
STB_EXTERN void stb__arr_setsize_(void **p, int size, int limit  STB__PARAMS);
STB_EXTERN void stb__arr_setlen_(void **p, int size, int newlen  STB__PARAMS);
STB_EXTERN void stb__arr_addlen_(void **p, int size, int addlen  STB__PARAMS);
STB_EXTERN void stb__arr_deleten_(void **p, int size, int loc, int n  STB__PARAMS);
STB_EXTERN void stb__arr_insertn_(void **p, int size, int loc, int n  STB__PARAMS);

#define stb_arr_free(p)            stb_arr_free_((void **) &(p))
#define stb__arr_copy              stb__arr_copy_

#ifndef STB_MALLOC_WRAPPER
  #define stb__arr_setsize         stb__arr_setsize_
  #define stb__arr_setlen          stb__arr_setlen_
  #define stb__arr_addlen          stb__arr_addlen_
  #define stb__arr_deleten         stb__arr_deleten_
  #define stb__arr_insertn         stb__arr_insertn_
#else
  #define stb__arr_addlen(p,s,n)    stb__arr_addlen_(p,s,n,__FILE__,__LINE__)
  #define stb__arr_setlen(p,s,n)    stb__arr_setlen_(p,s,n,__FILE__,__LINE__)
  #define stb__arr_setsize(p,s,n)   stb__arr_setsize_(p,s,n,__FILE__,__LINE__)
  #define stb__arr_deleten(p,s,i,n) stb__arr_deleten_(p,s,i,n,__FILE__,__LINE__)
  #define stb__arr_insertn(p,s,i,n) stb__arr_insertn_(p,s,i,n,__FILE__,__LINE__)
#endif

#ifdef STB_DEFINE
static void *stb__arr_context;

void *stb_arr_malloc_parent(void *p)
{
   void *q = stb__arr_context;
   stb__arr_context = p;
   return q;
}

void stb_arr_malloc(void **target, void *context)
{
   stb__arr *q = (stb__arr *) stb_malloc(context, sizeof(*q));
   q->len = q->limit = 0;
   q->stb_malloc = 1;
   q->signature = stb_arr_signature;
   *target = (void *) (q+1);
}

static void * stb__arr_malloc(int size)
{
   if (stb__arr_context)
      return stb_malloc(stb__arr_context, size);
   return malloc(size);
}

void * stb__arr_copy_(void *p, int elem_size)
{
   stb__arr *q;
   if (p == NULL) return p;
   q = (stb__arr *) stb__arr_malloc(sizeof(*q) + elem_size * stb_arrhead2(p)->limit);
   stb_arr_check2(p);
   memcpy(q, stb_arrhead2(p), sizeof(*q) + elem_size * stb_arrhead2(p)->len);
   q->stb_malloc = !!stb__arr_context;
   return q+1;
}

void stb_arr_free_(void **pp)
{
   void *p = *pp;
   stb_arr_check2(p);
   if (p) {
      stb__arr *q = stb_arrhead2(p);
      if (q->stb_malloc)
         stb_free(q);
      else
         free(q);
   }
   *pp = NULL;
}

static void stb__arrsize_(void **pp, int size, int limit, int len  STB__PARAMS)
{
   void *p = *pp;
   stb__arr *a;
   stb_arr_check2(p);
   if (p == NULL) {
      if (len == 0 && size == 0) return;
      a = (stb__arr *) stb__arr_malloc(sizeof(*a) + size*limit);
      a->limit = limit;
      a->len   = len;
      a->stb_malloc = !!stb__arr_context;
      a->signature = stb_arr_signature;
   } else {
      a = stb_arrhead2(p);
      a->len = len;
      if (a->limit < limit) {
         void *p;
         if (a->limit >= 4 && limit < a->limit * 2)
            limit = a->limit * 2;
         if (a->stb_malloc)
            p = stb_realloc(a, sizeof(*a) + limit*size);
         else
            #ifdef STB_MALLOC_WRAPPER
            p = stb__realloc(a, sizeof(*a) + limit*size, file, line);
            #else
            p = realloc(a, sizeof(*a) + limit*size);
            #endif
         if (p) {
            a = (stb__arr *) p;
            a->limit = limit;
         } else {
            // throw an error!
         }
      }
   }
   a->len   = stb_min(a->len, a->limit);
   *pp = a+1;
}

void stb__arr_setsize_(void **pp, int size, int limit  STB__PARAMS)
{
   void *p = *pp;
   stb_arr_check2(p);
   stb__arrsize_(pp, size, limit, stb_arr_len2(p)  STB__ARGS);
}

void stb__arr_setlen_(void **pp, int size, int newlen  STB__PARAMS)
{
   void *p = *pp;
   stb_arr_check2(p);
   if (stb_arrcurmax2(p) < newlen || p == NULL) {
      stb__arrsize_(pp, size, newlen, newlen  STB__ARGS);
   } else {
      stb_arrhead2(p)->len = newlen;
   }
}

void stb__arr_addlen_(void **p, int size, int addlen  STB__PARAMS)
{
   stb__arr_setlen_(p, size, stb_arr_len2(*p) + addlen  STB__ARGS);
}

void stb__arr_insertn_(void **pp, int size, int i, int n  STB__PARAMS)
{
   void *p = *pp;
   if (n) {
      int z;

      if (p == NULL) {
         stb__arr_addlen_(pp, size, n  STB__ARGS);
         return;
      }

      z = stb_arr_len2(p);
      stb__arr_addlen_(&p, size, n  STB__ARGS);
      memmove((char *) p + (i+n)*size, (char *) p + i*size, size * (z-i));
   }
   *pp = p;
}

void stb__arr_deleten_(void **pp, int size, int i, int n  STB__PARAMS)
{
   void *p = *pp;
   if (n) {
      memmove((char *) p + i*size, (char *) p + (i+n)*size, size * (stb_arr_len2(p)-(i+n)));
      stb_arrhead2(p)->len -= n;
   }
   *pp = p;
}

#endif

#endif //STB_INCLUDE_STB_H
