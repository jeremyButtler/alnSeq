/*########################################################
# Name: cStrToNumberFun
# Use:
#  - These are my numeric functions to convert numbers to
#    specific data types.
# Libraries:
#   - "dataTypeShortHand.h"
# C Standard Libraries:
# Note:
#  - Currently this is set up only for base 10.
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' cStrToNumberFun SOF: Start Of Header
'  - Functions to converts c-strings to numbers
'  o header:
'    - Has includes and the defined size of each number
'  o fun-01 base10StrToUL:
'    - Converts a base 10 c-string to an unsigned long
'  o fun-02 base10StrToUI:
'    - Converts a base 10 c-string to an unsigned integer
'  o fun-03 base10StrToUS:
'    - Converts a base 10 c-string to an unsigned short
'  o fun-04 base10StrToUC:
'    - Converts a base 10 c-string to an unsigned char
'  o fun-05 base10StrToU64:
'    - Converts a base 10 c-string to an uint64_t
'  o fun-06 base10StrToU32:
'    - Converts a base 10 c-string to an uint32_t
'  o fun-07 base10StrToU16:
'    - Converts a base 10 c-string to an uint16_t
'  o fun-08 base10StrToU8:
'    - Converts a base 10 c-string to an uint8_t
'  o fun-09 base10StrToUL:
'    - Converts a base 10 c-string to an unsigned long
'  o fun-10 base10StrToSI:
'    - Converts a base 10 c-string to an signed integer
'  o fun-11 base10StrToSS:
'    - Converts a base 10 c-string to an signed short
'  o fun-12 base10StrToSC:
'    - Converts a base 10 c-string to an signed char
'  o fun-13 base10StrToS64:
'    - Converts a base 10 c-string to an int64_t
'  o fun-14 base10StrToS32:
'    - Converts a base 10 c-string to an int32_t
'  o fun-15 base10StrToS16:
'    - Converts a base 10 c-string to an int16_t
'  o fun-16 base10StrToS8:
'    - Converts a base 10 c-string to an int8_t
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - Has includes and the defined size of each number
|  o header sec-01:
|    - Includes
|  o header sec-02:
|    - Maximum digits and values for unsigned data types
|  o header sec-03:
|    - Maximum digits and values for signed data types
\-------------------------------------------------------*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec-01:
^  - includes 
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#ifndef STRTONUM_H
#define STRTONUM_H

#include "dataTypeShortHand.h" /*unsigned var to uvar*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec-02:
^  - Maximum digits and values for unsigned data types
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*Some tricks I use that will be optimized out at
` compiler time:
*/
#define defMaxUL()(\
     ((ulong) 1 << ((sizeof(ulong) << 3) - 1))\
   | (((ulong) 1 << ((sizeof(ulong) << 3) - 1)) - 1)\
)

#define defNumBase10DigInUL()(\
      (char) (((sizeof(ulong) << 3)) / 3.333333333)\
)

#define defMaxUI()(\
     ((uint) 1 << ((sizeof(uint) << 3) - 1))\
   | (((uint) 1 << ((sizeof(uint) << 3) - 1)) - 1)\
)

#define defNumBase10DigInUI()(\
      (char) (((sizeof(uint) << 3)) / 3.333333333)\
)

#define defMaxUS()(\
     (1 << ((sizeof(ushort) << 3) - 1))\
   | ((1 << ((sizeof(ushort) << 3) - 1)) - 1)\
)

#define defNumBase10DigInUS()(\
      (char) (((sizeof(ushort) << 3)) / 3.333333333)\
)

#define defMaxUC()(\
     (1 << ((sizeof(uchar) << 3) - 1))\
   | ((1 << ((sizeof(uchar) << 3) - 1)) - 1)\
)

#define defNumBase10DigInUC()(\
      (char) (((sizeof(uchar) << 3)) / 3.333333333)\
)

/* (1 << ((sizeof(datatype) << 3) - 1)) |
`  ((1 << ((sizeof(datatype) << 3) - 1))) - 1)
`    - This gets me the maximum value of a datatype
`      without using a C library and avoides overflow
`      issues. This will be optimized out.
`    o sizeof(datatype) << 3
`      - gives me the number of bits
`    o 1 << ((sizeof(datatype) << 3) - 1)
`      - Gives me the highest power of two (2 ^ (bits - 1))
`    o highest power of 2 - 1
`      - Sets every bit to 1 except the highest bit
`        (highest power of two) (2 ^ (bits - 2))
`    o highest power of 2 | (highest power of 2 - 1)
`      - Sets all bits, thus is the maximum number
`
` (char) (((sizeof(datatype) << 3)) / 3.333333333)
`    - This gives me the number of digits that will not
`      overflow (number digits - 1) in a base 10 number.
`    o (sizeof(datatype) << 3
`      - Get number of bits in a datatype
`    o number bits / 3.333333333
`      - Converts powers of 2 to powers of 10
`    o (cast) powers 10
`      - Make sure is a whole number
*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec-03:
^  - Maximum digits and values for signed data types
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#define defMaxSL()(\
     ((long) 1 << ((sizeof(long) << 3) - 2))\
   | (((long) 1 << ((sizeof(long) << 3) - 2)) - 1)\
)

#define defNumBase10DigInSL()(\
      (char) (((sizeof(long) << 3) - 1) / 3.333333333)\
)

#define defMaxSI()(\
     ((int) 1 << ((sizeof(int) << 3) - 2))\
   | (((int) 1 << ((sizeof(int) << 3) - 2)) - 1)\
)

#define defNumBase10DigInSI()(\
      (char) (((sizeof(int) << 3) - 1) / 3.333333333)\
)

#define defMaxSS()(\
     (1 << ((sizeof(short) << 3) - 2))\
   | ((1 << ((sizeof(short) << 3) - 2)) - 1)\
)

#define defNumBase10DigInSS()(\
      (char) (((sizeof(short) << 3) - 1) / 3.333333333)\
)

#define defMaxSC()(\
     (1 << ((sizeof(char) << 3) - 2))\
   | ((1 << ((sizeof(char) << 3) - 2)) - 1)\
)

#define defNumBase10DigInSC()(\
      (char) (((sizeof(char) << 3) - 1) / 3.333333333)\
)

/* (1 << ((sizeof(datatype) << 3) - 2)) |
`  ((1 << ((sizeof(datatype) << 3) - 2))) - 1)
`    - This gets me the maximum value of a datatype
`      without using a C library and does not have
`      overflow issues. This will be optimized out.
`    o sizeof(datatype) << 3
`      - gives me the number of bits
`    o 1 << ((sizeof(datatype) << 3) - 2)
`      - Gives me the highest power of two (2 ^ (bits-1))
`    o 2^(bits - 1) - 1
`      - Sets every bit to 1 except the highest bit
`        (highest power of two)
`    o 2^(bits - 1) | (2^(bits - 1) - 1)
`      - Sets all bits, thus is the maximum number
`
` (char) (((sizeof(datatype) << 3) - 1) / 3.333333333)
`    - This gives me the number of digits that will not
`      overflow (number digits - 1) in a base 10 number.
`    o (sizeof(datatype) << 3
`      - Get number of bits in a datatype
`    o number bits - 1
`      - Gives me the number of bits in an signed datatype
`    o (number bits - 1) / 3.333333333
`      - Converts powers of 2 to powers of 10
`    o (cast) powers 10
`      - Make sure is a whole number
*/

/*-------------------------------------------------------\
| Fun-01: base10StrToUL
|  - Converts a base 10 c-string to an unsigned long.
|    Your system determines how many bits are in an long.
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retUL:
|    o an unisnged long to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUL to hold the unsigned long
| Note:
|  - This function will only convert the first number of
|    digits in an unsigned long or till non-numeric
|    character
\-------------------------------------------------------*/
#define base10StrToUL(inStr, retUL)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   ulong macroTmpUL = 0;\
   (retUL) = 0;\
   \
   for(\
      macroCDig = 0;\
      macroCDig < defNumBase10DigInUL();\
      ++macroCDig\
   ){ /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retUL) = (retUL) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a last digit (could overflow)*/\
      if((retUL) <= defMaxUL() / 10)\
      { /*If: number is < the max unsigned int size/10*/\
         macroTmpUL = 10 * (retUL) + *macroRetStr - 48;\
         \
         if(macroTmpUL > 9)\
         { /*If: I have room for 10th digit*/\
            (retUL) = macroTmpUL;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: number is < the max unsigned int size/10*/\
   } /*If have a last digit (could overflow)*/\
   \
   macroRetStr; /*Let the user see were ended at*/\
}) /*base10StrToUL*/

/*-------------------------------------------------------\
| Fun-02: base10StrToUI
|  - Converts a base 10 c-string to an unsigned integer.
|    Your system determines how many bits are in an int.
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retUI:
|    o an unisnged integer to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUI to hold the unsigned integer
| Note:
|  - This function will only convert the first number of
|    digits in an unsigned integer or till non-numeric
|    character
\-------------------------------------------------------*/
#define base10StrToUI(inStr, retUI)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   uint macroTmpUI = 0;\
   (retUI) = 0;\
   \
   for(\
      macroCDig = 0;\
      macroCDig < defNumBase10DigInUI();\
      ++macroCDig\
   ){ /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retUI) = (retUI) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a last digit (could overflow)*/\
      if((retUI) <= defMaxUI() / 10)\
      { /*If: number is < the max unsigned int size/10*/\
         macroTmpUI = 10 * (retUI) + *macroRetStr - 48;\
         \
         if(macroTmpUI > 9)\
         { /*If: I have room for 10th digit*/\
            (retUI) = macroTmpUI;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: number is < the max unsigned int size/10*/\
   } /*If have a last digit (could overflow)*/\
   \
   macroRetStr; /*Let the user see were ended at*/\
}) /*base10StrToUI*/

/*-------------------------------------------------------\
| Fun-03: base10StrToUS
|  - Converts a base 10 c-string to an unsigned short.
|    Your system determines how many bits are in an short.
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retUS:
|    o an unsigned short to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUS to hold the unsigned short
| Note:
|  - This function will only convert the first number of
|    digits in an unsigned short or till non-numeric
|    character
\-------------------------------------------------------*/
#define base10StrToUS(inStr, retUS)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   ushort macroTmpUS = 0;\
   (retUS) = 0;\
   \
   for(\
      macroCDig = 0;\
      macroCDig < defNumBase10DigInUS();\
      ++macroCDig\
   ){ /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retUS) = (retUS) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a last digit (could overflow)*/\
      if((retUS) <= defMaxUS() / 10)\
      { /*If: number is < the max unsigned int size/10*/\
         macroTmpUS = 10 * (retUS) + *macroRetStr - 48;\
         \
         if(macroTmpUS > 9)\
         { /*If: I have room for 10th digit*/\
            (retUS) = macroTmpUS;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: number is < the max unsigned int size/10*/\
   } /*If have a last digit (could overflow)*/\
   \
   macroRetStr; /*Let the user see were ended at*/\
}) /*base10StrToUS*/

/*-------------------------------------------------------\
| Fun-04: base10StrToUC
|  - Converts a base 10 c-string to an unsigned char.
|    Your system determines how many bits are in an char.
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retUC:
|    o an unsigned char to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUC to hold the unsigned char
| Note:
|  - This function will only convert the first number of
|    digits in an unsigned char or till non-numeric
|    character
\-------------------------------------------------------*/
#define base10StrToUC(inStr, retUC)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   uchar macroTmpUC = 0;\
   (retUC) = 0;\
   \
   for(\
      macroCDig = 0;\
      macroCDig < defNumBase10DigInUC();\
      ++macroCDig\
   ){ /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retUC) = (retUC) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a last digit (could overflow)*/\
      if((retUC) <= defMaxUC() / 10)\
      { /*If: number is < the max unsigned int size/10*/\
         macroTmpUC = 10 * (retUC) + *macroRetStr - 48;\
         \
         if(macroTmpUC > 9)\
         { /*If: I have room for 10th digit*/\
            (retUC) = macroTmpUC;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: number is < the max unsigned int size/10*/\
   } /*If have a last digit (could overflow)*/\
   \
   macroRetStr; /*Let the user see were ended at*/\
}) /*base10StrToUC*/

/*-------------------------------------------------------\
| Fun-05: base10StrToU64
|  - Converts a base 10 c-string to an 64 bit unsigned int
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retU64:
|    o 64 bit unisnged integer to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUI64 to hold the uint64_t integer
| Note:
|  - This function will only convert the first number of
|    digits in an uint64_t unsigned integer or till
|    non-numeric character
\-------------------------------------------------------*/
#define base10StrToU64(inStr, retU64)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   (retU64) = 0;\
   \
   for(macroCDig = 0; macroCDig < 19; ++macroCDig)\
   { /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retU64) = (retU64) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a 10th digit (could overflow)*/\
      if((retU64) <= (ulonglong) 184467440737095516)\
      { /*If: This number is < the max 64 bit size/10*/\
         /*I am not using a temporary variable here
         ` because I do not want this library to call
         ` stdint. I want the user to do that.
         */\
         if((ulong) 10 * (retU64) + *macroRetStr - 48 > 9)\
         { /*If: I have room for 10th digit*/\
            (retU64) = (retU64) * 10 + *macroRetStr - 48;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: This number is < the max 64bit size/10*/\
   } /*If have a 10th digit (could overflow)*/\
   \
   macroRetStr; /*Let the user see were ended at*/\
}) base10/*StrToU64*/

/*-------------------------------------------------------\
| Fun-06: base10StrToU32
|  - Converts a base 10 c-string to an 32 bit unsigned int
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retUI32:
|    o 32 bit unisnged integer to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retU32 to hold the uint32_t integer
| Note:
|  - This function will only convert the first number of
|    digits in an uint32_t unsigned integer or till
|    non-numeric character
\-------------------------------------------------------*/
#define base10StrToU32(inStr, retU32)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   (retU32) = 0;\
   \
   for(macroCDig = 0; macroCDig < 9; ++macroCDig)\
   { /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retU32) = (retU32) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a 10th digit (could overflow)*/\
      if((retU32) <= (ulong) 429496729)\
      { /*If: This number is < the max 32bit size/10*/\
         /*I am not using a temporary variable here because
         ` I do not want this library to call stdint. I
         ` want the user to do that.
         */\
         if((uint) (10 * (retU32) + *macroRetStr - 48 > 9)\
         { /*If: I have room for 10th digit*/\
            (retU32) = (retU32) * 10 + *macroRetStr - 48;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: This number is < the max 32bit size/10*/\
   } /*If have a 10th digit (could overflow)*/\
   \
   macroRetStr; /*Let the user see were ended at*/\
}) base10/*StrToU32*/

/*-------------------------------------------------------\
| Fun-07: base10StrToU16
|  - Converts a base 10 c-string to an 16 bit unsigned int
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retUI16:
|    o 16 bit unisnged integer to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retU16 to hold the uint16_t integer
| Note:
|  - This function will only convert the first number of
|    digits in an uint16_t unsigned integer or till
|    non-numeric character
\-------------------------------------------------------*/
#define base10StrToU16(inStr, retU16)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   (retU16) = 0;\
   \
   for(macroCDig = 0; macroCDig < 3; ++macroCDig)\
   { /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retU16) = (retU16) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a 10th digit (could overflow)*/\
      if((retU16) <= (ulong) 6553)\
      { /*If: This number is < the max 16bit size/10*/\
         /*I am not using a temporary variable here
         ` because I do not want this library to call
         ` stdint. I want the user to do that.
         */\
         if((uint) (10 * (retU16) + *macroRetStr - 48 > 9)\
         { /*If: I have room for 10th digit*/\
            (retU16) = (retU16) * 10 + *macroRetStr - 48;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: This number is < the max 32bit size/10*/\
   } /*If have a 10th digit (could overflow)*/\
   \
   macroRetStr; /*Let the user see were ended at*/\
}) base10/*StrToU16*/

/*-------------------------------------------------------\
| Fun-08: base10StrToU8
|  - Converts a base 10 c-string to an 8 bit unsigned int
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retUI8:
|    o 8 bit unisnged integer to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retU8 to hold the uint8_t integer
| Note:
|  - This function will only convert the first number of
|    digits in an uint8_t unsigned integer or till
|    non-numeric character
\-------------------------------------------------------*/
#define base10StrToU8(inStr, retU8)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   (retU8) = 0;\
   \
   for(macroCDig = 0; macroCDig < 3; ++macroCDig)\
   { /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retU8) = (retU8) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a 10th digit (could overflow)*/\
      if((retU8) <= 25)\
      { /*If: This number is < the max 8bit size/10*/\
         /*I am not using a temporary variable here
         ` because I do not want this library to call
         ` stdint. I want the user to do that.
         */\
         if((uint) (10 * (retU8) + *macroRetStr - 48 > 9)\
         { /*If: I have room for 10th digit*/\
            (retU8) = (retU8) * 10 + *macroRetStr - 48;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: This number is < the max 8bit size/10*/\
   } /*If have a 10th digit (could overflow)*/\
   \
   macroRetStr; /*Let the user see were ended at*/\
}) /*base10StrToU8*/

/*-------------------------------------------------------\
| Fun-09: base10StrToSL
|  - Converts a base 10 c-string to an signed long.
|    Your system determines how many bits are in an long.
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retUL:
|    o an signed long to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retSL to hold the signed long
| Note:
|  - This function will only convert the first number of
|    digits in an long or till non-numeric character
\-------------------------------------------------------*/
#define base10StrToSL(inStr, retSL)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   long macroTmpL = 0;\
   (retSL) = 0;\
   \
   if(*inStr == '-') ++macroRetStr;\
   \
   for(\
      macroCDig = 0;\
      macroCDig < defNumBase10DigInSL();\
      ++macroCDig\
   ){ /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retSL) = (retSL) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a last digit (could overflow)*/\
      if((retSL) <= defMaxSL() / 10)\
      { /*If: number is < the max unsigned int size/10*/\
         macroTmpL = 10 * (retSL) + *macroRetStr - 48;\
         \
         if(macroTmpL > 9)\
         { /*If: I have room for 10th digit*/\
            (retSL) = macroTmpL;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: number is < the max unsigned int size/10*/\
   } /*If have a last digit (could overflow)*/\
   \
   if(*inStr == '-') retSL = -retSL;\
   \
   macroRetStr;\
}) /*base10StrToSL*/

/*-------------------------------------------------------\
| Fun-10: base10StrToSI
|  - Converts a base 10 c-string to an signed integer
|    Your system determines number bits in an integer
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retUL:
|    o an signed long to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retSL to hold the signed integer
| Note:
|  - This function will only convert the first number of
|    digits in an integer or till non-numeric character
\-------------------------------------------------------*/
#define base10StrToSI(inStr, retSI)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   int macroTmpI = 0;\
   (retSI) = 0;\
   \
   if(*inStr == '-') ++macroRetStr;\
   \
   for(\
      macroCDig = 0;\
      macroCDig < defNumBase10DigInSI();\
      ++macroCDig\
   ){ /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retSI) = (retSI) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a last digit (could overflow)*/\
      if((retSI) <= defMaxSI() / 10)\
      { /*If: number is < the max unsigned int size/10*/\
         macroTmpI = 10 * (retSI) + *macroRetStr - 48;\
         \
         if(macroTmpI > 9)\
         { /*If: I have room for 10th digit*/\
            (retSI) = macroTmpI;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: number is < the max unsigned int size/10*/\
   } /*If have a last digit (could overflow)*/\
   \
   if(*inStr == '-') retSI = -retSI;\
   \
   macroRetStr;\
}) /*base10StrToSI*/

/*-------------------------------------------------------\
| Fun-11: base10StrToSS
|  - Converts a base 10 c-string to an signed short
|    Your system determines number bits in an short
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retUL:
|    o an signed long to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retSL to hold the signed short
| Note:
|  - This function will only convert the first number of
|    digits in an short or till non-numeric character
\-------------------------------------------------------*/
#define base10StrToSS(inStr, retSS)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   short macroTmpS = 0;\
   (retSS) = 0;\
   \
   if(*inStr == '-') ++macroRetStr;\
   \
   for(\
      macroCDig = 0;\
      macroCDig < defNumBase10DigInSS();\
      ++macroCDig\
   ){ /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retSS) = (retSS) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a last digit (could overflow)*/\
      if((retSS) <= defMaxSS() / 10)\
      { /*If: number is < the max unsigned int size/10*/\
         macroTmpS = 10 * (retSS) + *macroRetStr - 48;\
         \
         if(macroTmpS > 9)\
         { /*If: I have room for 10th digit*/\
            (retSS) = macroTmpS;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: number is < the max unsigned int size/10*/\
   } /*If have a last digit (could overflow)*/\
   \
   if(*inStr == '-') retSS = -retSS;\
   \
   macroRetStr;\
}) /*base10StrToSL*/

/*-------------------------------------------------------\
| Fun-12: base10StrToSC
|  - Converts a base 10 c-string to an signed char
|    Your system determines number bits in an char
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retUL:
|    o an signed long to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retSL to hold the signed char
| Note:
|  - This function will only convert the first number of
|    digits in an char or till non-numeric character
\-------------------------------------------------------*/
#define base10StrToSC(inStr, retSC)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   char macroTmpC = 0;\
   (retSC) = 0;\
   \
   if(*inStr == '-') ++macroRetStr;\
   \
   for(\
      macroCDig = 0;\
      macroCDig < defNumBase10DigInSC();\
      ++macroCDig\
   ){ /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retSC) = (retSC) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a last digit (could overflow)*/\
      if((retSC) <= defMaxSC() / 10)\
      { /*If: number is < the max unsigned int size/10*/\
         macroTmpC = 10 * (retSC) + *macroRetStr - 48;\
         \
         if(macroTmpC > 9)\
         { /*If: I have room for 10th digit*/\
            (retSC) = macroTmpC;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: number is < the max unsigned int size/10*/\
   } /*If have a last digit (could overflow)*/\
   \
   if(*inStr == '-') retSC = -retSC;\
   \
   macroRetStr;\
}) /*base10StrToSL*/

/*-------------------------------------------------------\
| Fun-13: base10StrToS64
|  - Converts a base 10 c-string to an 64 bit signed int
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retS64:
|    o 64 bit unisnged integer to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUI64 to hold the int64_t integer
| Note:
|  - This function will only convert the first number of
|    digits in an int64_t unsigned integer or till
|    non-numeric character
\-------------------------------------------------------*/
#define base10StrToS64(inStr, retS64)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   (retS64) = 0;\
   \
   if(*inStr == '-') ++macroRetStr;\
   \
   for(macroCDig = 0; macroCDig < 18; ++macroCDig)\
   { /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retS64) = (retS64) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a 10th digit (could overflow)*/\
      if((retS64) <= (longlong) 922337203685477580)\
      { /*If: This number is < the max 64 bit size/10*/\
         /*I am not using a temporary variable here
         ` because I do not want this library to call
         ` stdint. I want the user to do that.
         */\
         if(10 * (retS64) + *macroRetStr - 48 > 9)\
         { /*If: I have room for 10th digit*/\
            (retS64) = (retS64) * 10 + *macroRetStr - 48;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: This number is < the max 64bit size/10*/\
   } /*If have a 10th digit (could overflow)*/\
   \
   if(*inStr == '-') retS64 = -retS64;\
   \
   macroRetStr; /*Let the user see were ended at*/\
}) base10/*StrToS64*/

/*-------------------------------------------------------\
| Fun-14: base10StrToS32
|  - Converts a base 10 c-string to an 32 bit signed int
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retS32:
|    o 32 bit unisnged integer to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUI32 to hold the int32_t integer
| Note:
|  - This function will only convert the first number of
|    digits in an int32_t unsigned integer or till
|    non-numeric character
\-------------------------------------------------------*/
#define base10StrToS32(inStr, retS32)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   (retS32) = 0;\
   \
   if(*inStr == '-') ++macroRetStr;\
   \
   for(macroCDig = 0; macroCDig < 9; ++macroCDig)\
   { /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retS32) = (retS32) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a 10th digit (could overflow)*/\
      if((retS32) <= (ulong) 214748364)\
      { /*If: This number is < the max 32 bit size/10*/\
         /*I am not using a temporary variable here
         ` because I do not want this library to call
         ` stdint. I want the user to do that.
         */\
         if(10 * (retS32) + *macroRetStr - 48 > 9)\
         { /*If: I have room for 10th digit*/\
            (retS32) = (retS32) * 10 + *macroRetStr - 48;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: This number is < the max 32bit size/10*/\
   } /*If have a 10th digit (could overflow)*/\
   \
   if(*inStr == '-') retS32 = -retS32;\
   \
   macroRetStr; /*Let the user see were ended at*/\
}) base10/*StrToS32*/

/*-------------------------------------------------------\
| Fun-15: base10StrToS16
|  - Converts a base 10 c-string to an 16 bit signed int
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retS16:
|    o 16 bit unisnged integer to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUI16 to hold the int16_t integer
| Note:
|  - This function will only convert the first number of
|    digits in an int16_t unsigned integer or till
|    non-numeric character
\-------------------------------------------------------*/
#define base10StrToS16(inStr, retS16)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   (retS16) = 0;\
   \
   if(*inStr == '-') ++macroRetStr;\
   \
   for(macroCDig = 0; macroCDig < 4; ++macroCDig)\
   { /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retS16) = (retS16) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a 10th digit (could overflow)*/\
      if((retS16) <= (ulong) 3276)\
      { /*If: This number is < the max 16 bit size/10*/\
         /*I am not using a temporary variable here
         ` because I do not want this library to call
         ` stdint. I want the user to do that.
         */\
         if(10 * (retS16) + *macroRetStr - 48 > 9)\
         { /*If: I have room for 10th digit*/\
            (retS16) = (retS16) * 10 + *macroRetStr - 48;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: This number is < the max 16bit size/10*/\
   } /*If have a 10th digit (could overflow)*/\
   \
   if(*inStr == '-') retS16 = -retS16;\
   \
   macroRetStr; /*Let the user see were ended at*/\
}) base10/*StrToS16*/

/*-------------------------------------------------------\
| Fun-16: base10StrToS8
|  - Converts a base 10 c-string to an 8 bit signed int
| Input:
|  - inStr:
|    o C-string with base 10 number to convert
|  - retS8:
|    o 8 bit unisnged integer to hold the returned number
| Output:
|  - Returns:
|    o pionter to character after last converted number
|  - Modifies:
|    o retUI8 to hold the int8_t integer
| Note:
|  - This function will only convert the first number of
|    digits in an int8_t unsigned integer or till
|    non-numeric character
\-------------------------------------------------------*/
#define base10StrToS8(inStr, retS8)({\
   char macroCDig = 0;\
   char *macroRetStr = (inStr);\
   (retS8) = 0;\
   \
   if(*inStr == '-') ++macroRetStr;\
   \
   for(macroCDig = 0; macroCDig < 2; ++macroCDig)\
   { /*Loop: Convert digits with no overflow concern*/\
      if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
         (retS8) = (retS8) * 10 + *macroRetStr - 48;\
      else break;\
      \
      ++macroRetStr;\
   } /*Loop: Convert digits with no overflow concern*/\
   \
   if(*macroRetStr - 58 < 0 && *macroRetStr - 47 > 0)\
   { /*If have a 10th digit (could overflow)*/\
      if((retS8) <= 12)\
      { /*If: This number is < the max 8 bit size/10*/\
         /*I am not using a temporary variable here
         ` because I do not want this library to call
         ` stdint. I want the user to do that.
         */\
         if(10 * (retS8) + *macroRetStr - 48 > 9)\
         { /*If: I have room for 10th digit*/\
            (retS8) = (retS8) * 10 + *macroRetStr - 48;\
            ++macroRetStr;\
         } /*If: I have room for 10th digit*/\
      } /*If: This number is < the max 8bit size/10*/\
   } /*If have a 10th digit (could overflow)*/\
   \
   if(*inStr == '-') retS8 = -retS8;\
   \
   macroRetStr; /*Let the user see were ended at*/\
}) /*base10StrToS8*/

#endif
