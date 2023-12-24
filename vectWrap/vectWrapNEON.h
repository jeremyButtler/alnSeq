/*#########################################################
# Name: vecWrapNEONFun.h
# Use:
#  - Holds macro wrappers and variables for NEON 128 bit
#    intrinsics
# Note:
#  - I would avoid the min/max and the unsigned ifs. They
#    are not very good on SSE2, SSE4, or AVX2.
#########################################################*/

#ifndef VECTWRAPNEON_H
#define VECTWRAPNEON_H

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' vectWrapNEON.h TOC: Table Of Contents
'  - Holds macro wrapper for NEON 128 bit intrinsics
'  o sec-01:
'    - Definitions and variable declerations
'  o sec-02:
'    - Macros wrapping each intrinsic
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-01:
^  - Definitions and variable declerations
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include <arm_neon.h>
#include <stdint.h>

#define defVectBits 128
#define defVectBytes 16
#define defNum8BitElms 16
#define defNum16BitElms 8
#define defNum32BitElms 4
#define defMvMask8 0xffffffffffffffffull
#define defMvMask16 0xffffffffffffffffull
#define defMvMask32 0xffffffffffffffffull

typedef int8x16_t vectI8;  /*vector of bytes (8bits)*/
typedef int16x8_t vectI16; /*vector of bytes (16bits)*/
typedef int32x4_t vectI32; /*vector of bytes (32bits)*/

typedef uint8x16_t vectU8;  /*vector of bytes (8bits)*/
typedef uint16x8_t vectU16; /*vector of bytes (16bits)*/
typedef uint32x4_t vectU32; /*vector of bytes (32bits)*/

/*For the masks, I am going to use casting behind the
` scenes to get them to a single data type
*/
typedef int8x16_t mask8;
typedef int16x8_t mask16;
typedef int32x4_t mask32;

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-02:
^  - Macros
^  o sec-02 sub-01:
^    - Input/output functions (load, set, store)
^  o sec-02 sub-02:
^    - Output functions
^  o sec-02 sub-03:
^    - Comparison functions (if)
^  o sec-02 sub-04:
^    - Logical functions (andNot, and, or, xor)
^  o sec-02 sub-05:
^    - Bit manipulation functions (shifts)
^  o sec-02 sub-06:
^    - Math functions
^  o sec-02 sub-07:
^    - Casts [for NEON & masks] and Conversions
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*********************************************************\
* Sec-02 Sub-01:
*  - Input functions
*  o sec-02 sub-01 cat-01:
*    - aligned loading
*  o sec-02 sub-01 cat-02:
*    - Unaligned loading
*  o sec-02 sub-01 cat-03:
*    - make zero vectors
*  o sec-02 sub-01 cat-04:
*    - Make vectors of one element
*  o sec-02 sub-01 cat-05:
*    - Insert an element into an vector
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-01:
+  - Aligned loading
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define load_I8Ary_retVectI8(arrayI8) \
   (vld1q_s8((int8_t *) (arrayI8));)

#define load_I16Ary_retVectI16(arrayI16) \
  (vld1q_s16((int16_t *) (arrayI16));)

#define load_I32Ary_retVectI32(arrayI32) \
   (vld1q_s32((int32_t *) (arrayI32));)


#define load_U8Ary_retVectU8(arrayU8) \
  (vld1q_s8((uint8_t *) (arrayU8));)

#define load_U16Ary_retVectU16(arrayU16) \
   (vld1q_u16((uint16_t *) (arrayU16));)

#define load_U32Ary_retVectU32(arrayU32) \
   (vld1q_u32((uint32_t *) (arrayU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-02:
+  - Unaligned loading
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define loadu_I8Ary_retVectI8(arrayI8) \
   (vld1q_s8((int8_t *) (arrayI8));)

#define loadu_I16Ary_retVectI16(arrayI16) \
   (vld1q_s16((int16_t *) (arrayI16));)

#define loadu_I32Ary_retVectI32(arrayI32) \
   (vld1q_s32((int32_t *) (arrayI32));)


#define loadu_U8Ary_retVectU8(arrayU8) \
   (vld1q_u8((uint8_t *) (arrayU8));)

#define loadu_U16Ary_retVectU16(arrayU16) \
   (vld1q_u16((uint16_t *) (arrayU16));)

#define loadu_U32Ary_retVectU32(arrayU32) \
   (vld1q_u32((uint32_t *) (arrayU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-03:
+  - make zero vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define zero_retVectI8() (vdupq_n_s8(0))
#define zero_retVectI16() (vdupq_n_s16(0))
#define zero_retVectI32() (vdupq_n_s32(0))

#define zero_retVectU8() (vdupq_n_u8(0))
#define zero_retVectU16() (vdupq_n_u16(0))
#define zero_retVectU32() (vdupq_n_u32(0))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-04:
+  - Make vectors of one element
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define set1_I8_retVectI8(valC) (vdupq_n_s8((valC)))
#define set1_I16_retVectI16(valC) (vdupq_n_s16((valC)))
#define set1_I32_retVectI32(valC) (vdupq_n_s32((valC)))

#define set1_U8_retVectU8(valC) (vdupq_n_u8((valC)))
#define set1_U16_retVectU16(valC) (vdupq_n_u16((valC)))
#define set1_U32_retVectU32(valC) (vdupq_n_u32((valC)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-05:
+  - Insert an element into an vector
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define insert_I16_retVectI16(insVectI16, insI16, posI)\
   (vld1q_lane_s16((insI16), (insVectI16), posI);)

#define insert_I32_retVectI32(insI32, posI)\
   (vld1q_lane_s32((insI32), (insVectI32), posI);)

#define insert_U16_retVectU16(insVectU16, insU16, posU)\
   (vld1q_lane_u16((insU16), (insVectU16), posU);)

#define insert_U32_retVectU32(insVectU32,insU32,posU)\
   (vld1q_lane_u32((insU32), (insVectU32), posU);)

/*********************************************************\
* Sec-02 Sub-02:
*  - Output functions
*  o sec-02 sub-02 cat-01:
*    - Store vector ouput into an array (aligned)
*  o sec-02 sub-02 cat-02:
*    - Store vector ouput into an array (unaligned)
*  o sec-02 sub-02 cat-03:
*    - Extract one element
*  o sec-02 sub-02 cat-04:
*    - Store masks into longs
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-01:
+  - Store vector ouput into an array (aligned)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define store_vectI8_retAryI8(retAry, inVectI8) \
   (vstlq_s8((int8_t *) (retAry), (inVectI8));)

#define store_vectI16_retAryI16(retAry, inVectI16) \
  (vstlq_s16((int16_t *) (retAry), (inVectI16));)

#define store_vectI32_retAryI32(retAry, inVectI32) \
  (vstlq_s32((int32_t *) (retAry), (inVectI32));)


#define store_vectU8_retAryU8(retAry, inVectU8) \
  (vstlq_u8((uint8_t *) (retAry), (inVectU8));)

#define store_vectU16_retAryU16(retAry, inVectU16) \
  (vstlq_u16((uint16_t *) (retAry), (inVectU16));)

#define store_vectU32_retAryU32(retAry, inVectU32) \
  (vstlq_u32((uint32_t *) (retAry), (inVectU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-02:
+  - Store vector ouput into an array (unaligned)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define storeu_vectI8_retAryI8(retAry, inVectI8) \
   (vstlq_s8((int8_t *) (retAry), (inVectI8));)

#define storeu_vectI16_retAryI16(retAry, inVectI16) \
  (vstlq_s16((int16_t *) (retAry), (inVectI16));)

#define storeu_vectI32_retAryI32(retAry, inVectI32) \
  (vstlq_s32((int32_t *) (retAry), (inVectI32));)


#define storeu_vectU8_retAryU8(retAry, inVectU8) \
  (vstlq_u8((uint8_t *) (retAry), (inVectU8));)

#define storeu_vectU16_retAryU16(retAry, inVectU16) \
  (vstlq_u16((uint16_t *) (retAry), (inVectU16));)

#define storeu_vectU32_retAryU32(retAry, inVectU32) \
  (vstlq_u32((uint32_t *) (retAry), (inVectU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-03:
+  - Store masks into longs
+  - The logic for these store functions came from
+  - https://community.arm.com/arm-community-blogs/b/infrastructure-solutions-blog/posts/porting-x86-vector-bitmask-optimizations-to-arm-neon
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*storing in retUL because this is how AVX512 does it*/

#define store_mask8_retUL(retUL, mask8) ({ \
  uint16x8_t res16x8 = vreinterpretq_u16_s8((mask8)); \
  uint8x8_t res8x8 = vshrn_n_u16(res16x8, 4); \
  uint64x1_t res64x1 = vreinterpret_u64_u8(res8x8); \
  (retUL) = vget_lane_u64(res64x1, 0); \
})

#define store_mask16_retUL(retUL, mask16) ({ \
  uint32x4_t res32x4 = vreinterpretq_u32_s16((mask16));\
  uint16x4_t res16x4 = vshrn_n_u32(res32x4, 4); \
  uint64x1_t res64x1 = vreinterpret_u64_u16(res16x4); \
  (retUL) = vget_lane_u64(res64x1, 0); \
}) /*Mask is 16 bit integers*/

#define store_mask32_retUL(retUL, mask32) ({ \
   uint64x2_t res64x2 =vreinterpretq_u64_s32((mask32));\
   uint32x2_t res32x2 = vshrn_n_u64(tmp64x2, 4); \
   uint64x1_t res64x1 = vreinterpret_u64_u32(tmp32x2); \
   (retUL) = vget_lane_u64(res64x1, 0); \
}) /*Mask is 32 bit integers*/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-04:
+  - Extract an element
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define extract_vectI8_I8(inVectI8, posI)\
   (vgetq_lane_s8((inVectI8), (posI)))

#define extract_vectU8_U8(inVectU8, posI)\
   (vgetq_lane_u8((inVectU8), (posI)))

#define extract_vectI16_I16(inVectI16, posI)\
   (vgetq_lane_s16((inVectI16), (posI)))

#define extract_vectU16_IU6(inVectU16, posI)\
   (vgetq_lane_u16((inVectU16), (posI)))

#define extract_vectI32_I32(inVectI32, posI)\
   (vgetq_lane_s32((inVectI32), (posI)))

#define extract_vectU32_U32(inVectU32, posI)\
   (vgetq_lane_u32((inVectU32), (posI)))

/*********************************************************\
* Sec-02 Sub-03:
*  - Comparison functions (if)
*  o sec-02 sub-03 cat-01:
*    - Comparisons (return masks)
*  o sec-02 sub-03 cat-02:
*    - Unsinged equal comparisions (returns masks)
*  o sec-02 sub-03 cat-03:
*    - Unsinged < comparisions(SSE/AVX2 less efficent;mask)
*  o sec-02 sub-03 cat-04:
*    - Unsinged > comparisions (SSE/AVX2 less efficent)
*  o sec-02 sub-03 cat-05:
*    - Fix differences in population counts (total 1's)
*    - This is here because I am also supporting SSE2 and
*      NEON.
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-01:
+  - Comparisons (return masks)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*8-bit comparisions; returns a mask8*/
#define cmpeq_vectI8_retMask8(vectOneI8, vectTwoI8) \
   ((mask8) vceqq_s8((vectOneI8), (vectTwoI8));)

#define cmpgt_vectI8_retMask8(vectOneI8, vectTwoI8) \
  ((mask8) vcgtq_s8((vectOneI8), (vectTwoI8));)

#define cmplt_vectI8_retMask8(vectOneI8, vectTwoI8) \
   ((mask8) vcltq_s8((vectOneI8), (vectTwoI8));)

/*16-bit comparisions; returns a mask16*/

#define cmpeq_vectI16_retMask16(vectOneI16, vectTwoI16) \
   ((mask16) vceqq_s16((vectOneI16), (vectTwoI16));)

#define cmpgt_vectI16_retMask16(vectOneI16, vectTwoI16) \
   ((mask16) vcgtq_s16((vectOneI16), (vectTwoI16));)

#define cmplt_vectI16_retMask16(vectOneI16, vectTwoI16) \
   ((mask16) vcltq_s16((vectOneI16), (vectTwoI16));)

/*32-bit comparisions; returns a mask32*/

#define cmpeq_vectI32_retMask32(vectOneI32, vectTwoI32) \
   ((mask32) vceqq_s32((vectOneI32), (vectTwoI32));)

#define cmpgt_vectI32_retMask32(vectOneI32, vectTwoI32) \
   ((mask32) vcgtq_s32((vectOneI32), (vectTwoI32));)

#define cmplt_vectI32_retMask32(vectOneI32, vectTwoI32) \
   ((mask32) vcltq_s32((vectOneI32), (vectTwoI32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-02:
+  - Unsinged comparisions (SSE less efficent, masks)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cmpeq_vectU8_retMask8(vectOneU8, vectTwoU8) ({\
   (mask8) vreinterpretq_u8_s8(\
      vceqq_u8((vectOneU8),(vectTwoU8))\
   );\
})

#define cmpeq_vectU16_retMask16(vectOneU16, vectTwoU16) ({\
   (mask16) vreinterpretq_u16_s16(\
      vceqq_u16((vectOneU16), (vectTwoU16)\
   );\
})

#define cmpeq_vectU32_retMask32(vectOneU32, vectTwoU32) ({\
   (mask32) vreinterpretq_u32_s32(\
      vceqq_u32((vectOneU32), (vectTwoU32)\
   );\
})

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-03:
+  - Unsinged > comparisions (SSE less efficent; mask)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*This needs 9 operations & 2 additional vectors in SSE2*/

#define cmplt_vectU8_retMask8(vectOneU8, vectTwoU8) ({\
   (mask8) vreinterpretq_u8_s8(\
      vcltq_u8((vectOneU8),(vectTwoU8))\
   );\
})

#define cmplt_vectU16_retMask16(vectOneU16, vectTwoU16) ({\
   (mask16) vreinterpretq_u16_s16(\
      vcltq_u16((vectOneU16), (vectTwoU16)\
   );\
})

#define cmplt_vectU32_retMask32(vectOneU32, vectTwoU32) ({\
   (mask32) vreinterpretq_u32_s32(\
      vcltq_u32((vectOneU32), (vectTwoU32)\
   );\
})

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-04:
+  - Unsinged greater than comparisions (less efficent)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*This needs 9 operations & 2 additional vectors in SSE2*/

#define cmplt_vectU8_retMask8(vectOneU8, vectTwoU8) ({\
   (mask8) vreinterpretq_u8_s8(\
      vcgtq_u8((vectOneU8),(vectTwoU8))\
   );\
})

#define cmplt_vectU16_retMask16(vectOneU16, vectTwoU16) ({\
   (mask16) vreinterpretq_u16_s16(\
      vcgtq_u16((vectOneU16), (vectTwoU16)\
   );\
})

#define cmplt_vectU32_retMask32(vectOneU32, vectTwoU32) ({\
   (mask32) vreinterpretq_u32_s32(\
      vcgtq_u32((vectOneU32), (vectTwoU32)\
   );\
})

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-05:
+  - Fix differences in population counts (total 1's)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/* This is here to support SSE2, AVX2, and NEON. 
`  The returned mask in AVX128 has only one bit per
`  data type in comparison
`  SSE (all) and AVX2 returns a 128 bit vector that is
`   converted to characters
*/
#define fix_mask8_popcount(inUL)  ((inUL) >> 2)
#define fix_mask16_popcount(inUL) ((inUL) >> 3)
#define fix_mask32_popcount(inUL) ((inUL) >> 4)

/*********************************************************\
* Sec-02 Sub-04:
*  - Logical functions (andNot, and, or, xor)
*  o sec-02 sub-02 cat-01:
*    - and not functions
*  o sec-02 sub-02 cat-02:
*    - and functions
*  o sec-02 sub-02 cat-03:
*    - or functions
*  o sec-02 sub-02 cat-04:
*    - xor functions
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-01:
+  - and not functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define andNot_vectI8_retVectI8(notVectI8, inVectI8) \
   (vandq_s8(vmvnq_s8((notVectI8)), (inVectI8));)

#define andNot_vectU8_retVectU8(notVectU8, inVectU8) \
   (vandq_u8(vmvnq_u8((notVectU8)), (inVectU8));)

#define andNot_vectI16_retVectI16(notVectI16, inVectI16) \
   (vandq_s16(vmvnq_s16((notVectI16)), (inVectI16));)

#define andNot_vectU16_retVectU16(notVectU16, inVectU16) \
   (vandq_u16(vmvnq_u16((notVectU16)), (inVectI16));)

#define andNot_vectI32_retVectI32(notVectI32, inVectI32) \
   (vandq_s32(vmvnq_s32((notVectI32)), (inVectI32));)

#define andNot_vectU32_retVectU32(notVectU32, inVectU32) \
   (vandq_u32(vmvnq_u32((notVectU32)), (inVectU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-02:
+  - and functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define and_vectI8_retVectI8(vectOneI8, vectTwoI8) \
    (vandq_s8((vectOneI8), (vectTwoI8));)

#define and_vectU8_retVectU8(vectOneU8, vectTwoU8) \
    (vandq_u8((vectOneU8), (vectTwoU8));)

#define and_vectI16_retVectI16(vectOneI16, vectTwoI16) \
    (vandq_s16((vectOneI16), (vectTwoI16));)

#define and_vectU16_retVectU16(vectOneU16, vectTwoU16) \
    (vandq_u16((vectOneU16), (vectTwoU16));)

#define and_vectI32_retVectI32(vectOneI32, vectTwoI32) \
    (vandq_s32((vectOneI32), (vectTwoI32));)

#define and_vectU32_retVectU32(vectOneU32, vectTwoU32) \
    (vandq_u32((vectOneU32), (vectTwoU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-03:
+  - or functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define or_vectI8_retVectI8(vectOneI8, vectTwoI8) \
    (vorrq_s8((vectOneI8), (vectTwoI8));)

#define or_vectU8_retVectU8(vectOneU8, vectTwoU8) \
    (vorrq_u8((vectOneU8), (vectTwoU8));)

#define or_vectI16_retVectI16(vectOneI16, vectTwoI16) \
    (vorrq_s16((vectOneI16), (vectTwoI16));)

#define or_vectU16_retVectU16(vectOneU16, vectTwoU16) \
    (vorrq_u16((vectOneU16), (vectTwoU16));)

#define or_vectI32_retVectI32(vectOneI32, vectTwoI32) \
    (vorrq_s32((vectOneI32), (vectTwoI32));)

#define or_vectU32_retVectU32(vectOneU32, vectTwoU32) \
    (vorrq_u32((vectOneU32), (vectTwoU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-04:
+  - xor functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define xor_vectI8_retVectI8(vectOneI8, vectTwoI8) \
    (veorq_s8((vectOneI8), (vectTwoI8));)

#define xor_vectU8_retVectU8(vectOneU8, vectTwoU8) \
    (veorq_u8((vectOneU8), (vectTwoU8));)

#define xor_vectI16_retVectI16(vectOneI16, vectTwoI16) \
    (veorq_s16((vectOneI16), (vectTwoI16));)

#define xor_vectU16_retVectU16(vectOneU16, vectTwoU16) \
    (veorq_u16((vectOneU16), (vectTwoU16));)

#define xor_vectI32_retVectI32(vectOneI32, vectTwoI32) \
    (veorq_s32((vectOneI32), (vectTwoI32));)

#define xor_vectU32_retVectU32(vectOneU32, vectTwoU32) \
    (veorq_u32((vectOneU32), (vectTwoU32));)

/*********************************************************\
* Sec-02 Sub-05:
*  - Bit manipulation functions (shifts)
*  o sec-02 sub-05 cat-01:
*    - Shift numbers in vectors right by x bits
*  o sec-02 sub-05 cat-03:
*    - Shift a vector right by x bytes
*  o sec-02 sub-05 cat-02:
*    - Shift numbers in vectors left by x bits
*  o sec-02 sub-05 cat-04:
*    - Shift a vector left by x bytes
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-01:
+  - Shift each element right by x bytes
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sr_vectI8_retVectI8(inVectI8, numBitsI)\
    (vshrq_s8((inVectI8), (numBitsI));)

#define sr_vectU8_retVectU8(inVectU8, numBitsI)\
    (vshrq_u8((inVectU8), (numBitsI));)

#define sr_vectI16_retVectI16(inVectI16, numBitsI)\
    (vshrq_s16((inVectI16), (numBitsI));)

#define sr_vectU16_retVectU16(inVectU16, numBitsI)\
    (vshrq_u16((inVectU16), (numBitsI));)

#define sr_vectI32_retVectI32(inVectI32, numBitsI)\
    (vshrq_s32((inVectI32), (numBitsI));)

#define sr_vectU32_retVectU32(inVectU32, numBitsI)\
    (vshrq_u32((inVectU32), (numBitsI));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-02:
+  - Shift a vector right by x bytes
+  - https://stackoverflow.com/questions/11854819/attempt-to-convert-sse2-fast-corner-score-code-to-arm-neon
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define srvect_vectI8_retVectI8(inVectI8, numBytesI)({\
    vectI8 zeroVect = zero_retVectI8();\
    vextq_s8(zeroVect, inVectI8, 16 - numBytesI);\
      /*Copies 8 - numBytesI zeros from zeroVect to the
      ` the 8 - numBytesI positions in inVectI8.
      ` the numBytesI position controls how many elements
      ` are kept in inVectI8. So 16 - numBytesI ensures
      ` that the correct number of bytes are kept.
      */\
})

#define srvect_vectU8_retVectU8(inVectU8, numBytesI)({\
    vectU8 zeroVect = zero_retVectU8();\
    vextq_u8(zeroVect, inVectU8, 16 - numBytesU);\
})

#define srvect_vectI16_retVectI16(inVectI16, numShortsI)({\
    vectI16 zeroVect = zero_retVectI16();\
    vextq_s16(zeroVect, inVectI16, 8 - numShortsI);\
})

#define srvect_vectU16_retVectU16(inVectU16, numShortsI)({\
    vectU16 zeroVect = zero_retVectU16();\
    vextq_u16(zeroVect, inVectU16, 8 - numShortsI);\
})

#define srvect_vectI32_retVectI32(inVectI32, numIntsI)({\
    vectI32 zeroVect = zero_retVectI32();\
    vextq_s32(zeroVect, inVectI32, 4 - numIntsI);\
})

#define srvect_vectU32_retVectU32(inVectU32, numIntsI)({\
    vectU32 zeroVect = zero_retVectU32();\
    vextq_u32(zeroVect, inVectU32, 4 - numIntsI);\
})

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-03:
+  - Shift numbers in vectors left by x bits
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sl_vectI8_retVectI8(inVectI8, numBitsI)\
    (vshlq_s8((inVectI8), (numBitsI));)

#define sl_vectU8_retVectU8(inVectU8, numBitsI)\
    (vshlq_u8((inVectU8), (numBitsI));)

#define sl_vectI16_retVectI16(inVectI16, numBitsI)\
    (vshlq_s16((inVectI16), (numBitsI));)

#define sl_vectU16_retVectU16(inVectU16, numBitsI)\
    (vshlq_u16((inVectU16), (numBitsI));)

#define sl_vectI32_retVectI32(inVectI32, numBitsI)\
    (vshlq_s32((inVectI32), (numBitsI));)

#define sl_vectU32_retVectU32(inVectU32, numBitsI)\
    (vshlq_u32((inVectU32), (numBitsI));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-04:
+  - Shift a vector left by x bytes
+  - https://stackoverflow.com/questions/11854819/attempt-to-convert-sse2-fast-corner-score-code-to-arm-neon
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define slvect_vectI8_retVectI8(inVectI8, numBytesI)({\
    vectI8 zeroVect = zero_retVectI8();\
    vextq_s8(inVectI8, zeroVect, numBytesI);\
      /*Copies numBytesI zeros from zeroVect to the
      ` the last two positions in inVectI8.
      ` the numBytesI position controls how many elements
      ` are discarded in inVect.
      */\
})

#define slvect_vectU8_retVectU8(inVectU8, numBytesI)({\
    vectU8 zeroVect = zero_retVectU8();\
    vextq_u8(inVectU8, zeroVect, numBytesU);\
})

#define slvect_vectI16_retVectI16(inVectI16, numShortsI)({\
    vectI16 zeroVect = zero_retVectI16();\
    vextq_s16(inVectI16, zeroVect, umShortsI);\
})

#define slvect_vectU16_retVectU16(inVectU16, numShortsI)({\
    vectU16 zeroVect = zero_retVectU16();\
    vextq_u16(inVectU16, zeroVect, umShortsI);\
})

#define slvect_vectI32_retVectI32(inVectI32, numIntsI)({\
    vectI32 zeroVect = zero_retVectI32();\
    vextq_s32(inVectI32, zeroVect, umIntsI);\
})

#define slvect_vectU32_retVectU32(inVectU32, numIntsI)({\
    vectU32 zeroVect = zero_retVectU32();\
    vextq_u32(inVectU32, zeroVect, umIntsI);\
})

/*********************************************************\
* Sec-02 Sub-06:
*  - Math functions
*  o sec-02 sub-06 cat-01:
*    - addition max 64 bit; [sse2-avx128]; no epu
*  o sec-02 sub-06 cat-02:
*    - staturation addition; max 16 bit [sse2-avx128]; epu
*  o sec-02 sub-06 cat-03:
*    - subtraction; max 64bit [sse2-avx128]; no epu
*  o sec-02 sub-06 cat-04:
*    - staturation subtraction; max 16 bit; [sse2-avx128];u
*  o sec-02 sub-06 cat-05:
*    - Max
*  o sec-02 sub-06 cat-06:
*    - min
*  o sec-02 sub-06 cat-07:
*    - Horizontal max
*  o sec-02 sub-06 cat-08:
*    - Horizontal min
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-01:
+  - addition max 64 bit; [sse2-avx128]; no epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define add_vectI8_retVectI8(vectOneI8, vectTwoI8) \
   (vaddq_s8((vectOneI8), (vectTwoI8));)

#define add_vectI16_retVectI16(vectOneI16, vectTwoI16) \
   (vaddq_s16((vectOneI16), (vectTwoI16));)

#define add_vectI32_retVectI32(vectOneI32, vectTwoI32) \
   (vaddq_s32((vectOneI32), (vectTwoI32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-02:
+  - staturation addition; max 16 bit [sse2 to avx128]; epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define addSat_vectU8_retVectU8(vectOneU8, vectTwoU8) \
   (vqaddq_u8((vectOneU8), (vectTwoU8));)

#define addSat_vectU16_retVectU16(vectOneU16, vectTwoU16)\
   (vqaddq_u16((vectOneU16), (vectTwoU16));)

#define addSat_vectI8_retVectI8(vectOneI8, vectTwoI8) \
   (vqaddq_s8((vectOneI8), (vectTwoI8));)

#define addSat_vectI16_retVectI16(vectOneI16, vectTwoI16) \
   (vqaddq_s16((vectOneI16), (vectTwoI16));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-03:
+  - subtraction; max 64bit [sse2-avx128]; no epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sub_vectI8_retVectI8(vectOneI8, vectTwoI8) \
     (vsubq_s8((vectOneI8),(vectTwoI8));)

#define sub_vectI16_retVectI16(vectOneI16, vectTwoI16) \
     (vsubq_s16((vectOneI16),(vectTwoI16));)

#define sub_vectI32_retVectI32(vectOneI32, vectTwoI32) \
     (vsubq_s32((vectOneI32),(vectTwoI32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-04:
+  - staturation subtraction; max 16 bit; [sse2-avx128]); u
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define subSat_vectU8_retVectU8(vectOneU8, vectTwoU8) \
   (vqsubq_u8((vectOneU8), (vectTwoU8));)

#define subSat_vectU16_retVectU16(vectOneU16, vectTwoU16) \
   (vqsubq_u16((vectOneU16), (vectTwoU16));)

#define subSat_vectI8_retVectI8(vectOneI8, vectTwoI8) \
   (vqsubq_s8((vectOneI8), (vectTwoI8));)

#define subSat_vectI16_retVectI16(vectOneI16, vectTwoI16) \
   (vqsubq_s16((vectOneI16), (vectTwoI16));)

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-05:
+  - Max
\+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define max_vectI8_retVectI8(firstVectI8,secVectI8)\
   (vmaxq_s8(firstVectI8, secVectI8);)

#define max_vectU8_retVectU8(firstVectU8,secVectU8)\
   (vmaxq_u8(firstVectU8, secVectU8);)

#define max_vectI16_retVectI16(firstVectI16,secVectI16)\
   (vmaxq_s16(firstVectI16, secVectI16);)

#define max_vectU16_retVectU16(firstVectU16,secVectU16)\
   (vmaxq_u16(firstVectU16, secVectU16);)

#define max_vectI32_retVectI32(firstVectI32,secVectI32)\
   (vmaxq_s32(firstVectI32, secVectI32);)

#define max_vectU32_retVectU32(firstVectU32,secVectU32)\
   (vmaxq_u32(firstVectU32, secVectU32);)

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-06:
+  - Min
\+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define min_vectI8_retVectI8(firstVectI8,secVectI8)\
   (vminq_s8(firstVectI8, secVectI8);)

#define min_vectU8_retVectU8(firstVectU8,secVectU8)\
   (vminq_u8(firstVectU8, secVectU8);)

#define min_vectI16_retVectI16(firstVectI16,secVectI16)\
   (vminq_s16(firstVectI16, secVectI16);)

#define min_vectU16_retVectU16(firstVectU16,secVectU16)\
   (vminq_u16(firstVectU16, secVectU16);)

#define min_vectI32_retVectI32(firstVectI32,secVectI32)\
   (vminq_s32(firstVectI32, secVectI32);)

#define min_vectU32_retVectU32(firstVectU32,secVectU32)\
   (vminq_u32(firstVectU32, secVectU32);)

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-07:
+  - Horizontal max
\+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define hmax_vectI8_retI8(retI8, inVectI8)\
   (retI8 = vmaxvq_s8((inVectI8)))

#define hmax_vectU8_retU8(retU8, inVectU8)\
   (retU8 = vmaxvq_u8((inVectU8)))

#define hmax_vectI16_retI16(retI16, inVectI16)\
   (retI16 = vmaxvq_s16((inVectI16)))

#define hmax_vectU16_retU16(retU16, inVectU16)\
   (retU16 = vmaxvq_u16((inVectU16)))

#define hmax_vectI32_retI32(retI32, inVectI32)\
   (retI32 = vmaxvq_s32((inVectI32)))

#define hmax_vectU32_retU32(retU32, inVectU32)\
   (retU32 = vmaxvq_u32((inVectU32)))

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-08:
+  - Horizontal min
\+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define hmin_vectI8_retI8(retI8, inVectI8)\
   (retI8 = vminvq_s8((inVectI8)))

#define hmin_vectU8_retU8(retU8, inVectU8)\
   (retU8 = vminvq_u8((inVectU8)))

#define hmin_vectI16_retI16(retI16, inVectI16)\
   (retI16 = vminvq_s16((inVectI16)))

#define hmin_vectU16_retU16(retU16, inVectU16)\
   (retU16 = vminvq_u16((inVectU16)))

#define hmin_vectI32_retI32(retI32, inVectI32)\
   (retI32 = vminvq_s32((inVectI32)))

#define hmin_vectU32_retU32(retU32, inVectU32)\
   (retU32 = vminvq_u32((inVectU32)))

/*********************************************************\
* Sec-02 Sub-07:
*  - Conversions [Here for NEON support]
*  o sec02 sub-07 cat-01:
*    - Convert signed vectors to other signed vectors
*  o sec02 sub-07 cat-02:
*    - Convert unsigned vectors to unsigned vectors
*  o sec02 sub-07 cat-03:
*    - Convert unsigned vectors to signed vectors
*  o sec02 sub-07 cat-04:
*    - Convert signed vectors to unsigned vectors
*  o sec02 sub-07 cat-05:
*    - Convert masks to vectors
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-01:
+  - Convert signed vectors to other signed vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cast_vectI8_to_VectI16(inVectI8) \
  (vreinterpretq_s8_s16((inVectI8));)

#define cast_vectI8_to_VectI32(inVectI8) \
  (vreinterpretq_s8_s32((inVectI8));)

#define cast_vectI16_to_VectI8(inVectI16) \
  (vreinterpretq_s16_s8((inVectI16));)

#define cast_vectI16_to_VectI32(inVectI16) \
  (vreinterpretq_s16_s32((inVectI16));)

#define cast_vectI32_to_VectI8(inVectI32) \
  (vreinterpretq_s32_s8((inVectI32));)

#define cast_vectI32_to_VectI16(inVectI32) \
  (vreinterpretq_s32_s16((inVectI32));)

/*Pointer casting (hope works)*/

#define ptrcast_vectI8_to_VectI16(inVectI8) \
  (vreinterpretq_s8_s16((inVectI8));)

#define ptrcast_vectI8_to_VectI32(inVectI8) \
  (vreinterpretq_s8_s32((inVectI8));)

#define ptrcast_vectI16_to_VectI8(inVectI16) \
  (vreinterpretq_s16_s8((inVectI16));)

#define ptrcast_vectI16_to_VectI32(inVectI16) \
  (vreinterpretq_s16_s32((inVectI16));)

#define ptrcast_vectI32_to_VectI8(inVectI32) \
  (vreinterpretq_s32_s8((inVectI32));)

#define ptrcast_vectI32_to_VectI16(inVectI32) \
  (vreinterpretq_s32_s16((inVectI32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-02:
+  - Convert unsigned vectors to unsigned vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cast_vectU8_to_VectU16(inVectU8) \
  (vreinterpretq_u8_u16((inVectU8));)

#define cast_vectU8_to_VectU32(inVectU8) \
  (vreinterpretq_u8_u32((inVectU8));)

#define cast_vectU16_to_VectU8(inVectU16) \
  (vreinterpretq_u16_u8((inVectU16));)

#define cast_vectU16_to_VectU32(inVectU16) \
  (vreinterpretq_u16_u32((inVectU16));)

#define cast_vectU32_to_VectU8(inVectU32) \
  (vreinterpretq_u32_u8((inVectU32));)

#define cast_vectU32_to_VectU16(inVectU32) \
  (vreinterpretq_u32_u16((inVectU32));)

/*Pointer casting (hope works)*/

#define ptrcast_vectU8_to_VectU16(inVectU8) \
  (vreinterpretq_u8_u16((inVectU8));)

#define ptrcast_vectU8_to_VectU32(inVectU8) \
  (vreinterpretq_u8_u32((inVectU8));)

#define ptrcast_vectU16_to_VectU8(inVectU16) \
  (vreinterpretq_u16_u8((inVectU16));)

#define ptrcast_vectU16_to_VectU32(inVectU16) \
  (vreinterpretq_u16_u32((inVectU16));)

#define ptrcast_vectU32_to_VectU8(inVectU32) \
  (vreinterpretq_u32_u8((inVectU32));)

#define ptrcast_vectU32_to_VectU16(inVectU32) \
  (vreinterpretq_u32_u16((inVectU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-03:
+  - Convert unsigned vectors to signed vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cast_vectU8_to_VectI8(inVectU8) \
  (vreinterpretq_u8_s8((inVectU8));)

#define cast_vectU8_to_VectI16(inVectU8) \
  (vreinterpretq_u8_s16((inVectU8));)

#define cast_vectU8_to_VectI32(inVectU8) \
  (vreinterpretq_u8_s32((inVectU8));)

#define cast_vectU16_to_VectI8(inVectU16) \
  (vreinterpretq_u16_s8((inVectU16));)

#define cast_vectU16_to_VectI16(inVectU16) \
  (vreinterpretq_u16_s16((inVectU16));)

#define cast_vectU16_to_VectI32(inVectU16) \
  (vreinterpretq_u16_s32((inVectU16));)

#define cast_vectU32_to_VectI8(inVectU32) \
  (vreinterpretq_u32_s8((inVectU32));)

#define cast_vectU32_to_VectI16(inVectU32) \
  (vreinterpretq_u32_s16((inVectU32));)

#define cast_vectU32_to_VectI32(inVectU32) \
  (vreinterpretq_u32_s32((inVectU32));)

/*Pointer casting (hope works)*/

#define ptrcast_vectU8_to_VectI8(inVectU8) \
  (vreinterpretq_u8_s8((inVectU8));)

#define ptrcast_vectU8_to_VectI16(inVectU8) \
  (vreinterpretq_u8_s16((inVectU8));)

#define ptrcast_vectU8_to_VectI32(inVectU8) \
  (vreinterpretq_u8_s32((inVectU8));)

#define ptrcast_vectU16_to_VectI8(inVectU16) \
  (vreinterpretq_u16_s8((inVectU16));)

#define ptrcast_vectU16_to_VectI16(inVectU16) \
  (vreinterpretq_u16_s16((inVectU16));)

#define ptrcast_vectU16_to_VectI32(inVectU16) \
  (vreinterpretq_u16_s32((inVectU16));)

#define ptrcast_vectU32_to_VectI8(inVectU32) \
  (vreinterpretq_u32_s8((inVectU32));)

#define ptrcast_vectU32_to_VectI16(inVectU32) \
  (vreinterpretq_u32_s16((inVectU32));)

#define ptrcast_vectU32_to_VectI32(inVectU32) \
  (vreinterpretq_u32_s32((inVectU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-04:
+  - Convert signed vectors to unsigned vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cast_vectI8_to_VectU8(inVectI8) \
  (vreinterpretq_s8_u8((inVectI8));)

#define cast_vectI8_to_VectU16(inVectI8) \
  (vreinterpretq_s8_u16((inVectI8));)

#define cast_vectI8_to_VectU32(inVectI8) \
  (vreinterpretq_s8_u32((inVectI8));)

#define cast_vectI16_to_VectU8(inVectI16) \
  (vreinterpretq_s16_u8((inVectI16));)

#define cast_vectI16_to_VectU16(inVectI16) \
  (vreinterpretq_s16_u16((inVectI16));)

#define cast_vectI16_to_VectU32(inVectI16) \
  (vreinterpretq_s16_u32((inVectI16));)

#define cast_vectI32_to_VectU8(inVectI32) \
  (vreinterpretq_s32_u8((inVectI32));)

#define cast_vectI32_to_VectU16(inVectI32) \
  (vreinterpretq_s32_u16((inVectI32));)

#define cast_vectI32_to_VectU32(inVectI32) \
  (vreinterpretq_s32_u32((inVectI32));)

/*Pointer casting (hope works)*/

#define ptrcast_vectI8_to_VectU8(inVectI8) \
  (vreinterpretq_s8_u8((inVectI8));)

#define ptrcast_vectI8_to_VectU16(inVectI8) \
  (vreinterpretq_s8_u16((inVectI8));)

#define ptrcast_vectI8_to_VectU32(inVectI8) \
  (vreinterpretq_s8_u32((inVectI8));)

#define ptrcast_vectI16_to_VectU8(inVectI16) \
  (vreinterpretq_s16_u8((inVectI16));)

#define ptrcast_vectI16_to_VectU16(inVectI16) \
  (vreinterpretq_s16_u16((inVectI16));)

#define ptrcast_vectI16_to_VectU32(inVectI16) \
  (vreinterpretq_s16_u32((inVectI16));)

#define ptrcast_vectI32_to_VectU8(inVectI32) \
  (vreinterpretq_s32_u8((inVectI32));)

#define ptrcast_vectI32_to_VectU16(inVectI32) \
  (vreinterpretq_s32_u16((inVectI32));)

#define ptrcast_vectI32_to_VectU32(inVectI32) \
  (vreinterpretq_s32_u32((inVectI32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-05:
+  - Convert masks to vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define mask8_retVectI8(inMask8) ((vectI8) (inMask8);)
#define mask16_retVectI16(inMask16) ((vectI16) (inMask16);)
#define mask32_retVectI32(inMask32) ((vectI32) (inMask32);)

#define mask8_retVectU8(inMask8) \
   (cast_vectI8_to_VectU8((vectI8) (inMask8));)

#define mask16_retVectU16(inMask16)\
   (cast_vectI16_to_VectU16((vectI16) (inMask16));)

#define mask32_retVectU32(inMask32)\
   (cast_vectI32_to_VectU32((vectI32) (inMask32));)

#endif
