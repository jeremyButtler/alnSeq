/*#######################################################\
# Name: genMath
#   - Contains math functions (currently max)
# Libraries:
# C Standard Libraries:
\#######################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   - Math functions
'   o header:
'     - definitoins and maximums
'   o fun-01 macroMax:
'     - Maximize function
'   o fun-02 macroIfMax:
'     - Maximize by one score, then select by a second
'       score
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef GEN_MATH_H
#define GEN_MATH_H

/*--------------------------------------------------------\
| Fun-01: macroMax
| Use:
|  - Find the maximum value (branchless)
| Input:
|  o ret:
|    - Value to hold the maximum value
|  o x:
|    - First value to compare, ret is set to x if x >= y
|  o y:
|    - second value to compare, ret is set to y if y > x
| Output:
|  - Sets:
|    - Sets ret to x if x is greater than or equal to y
|    - Sets ret to y if y is greather than x
| From:
|  - https://graphics.stanford.edu/~seander/bithacks.html#IntegerMinOrMax
\--------------------------------------------------------*/
#define macroMax(ret, x, y){\
   (ret) = (x) ^ (((x) ^ (y)) & (-((x) < (y))));\
   /*Logic:
   `  x < y:
   `    if x < y I get 0 (x > y); else I get 1
   `  -(x < y):
   `    If x < y I get -1 (111...)
   `    If x >= y I get 0
   `  x ^ ((x ^ y) & 111...): (only when x < y)
   `    This gives me x
   `  x ^ (x ^ y) & 0: (only when x >= y)
   `    This gives me y
   */\
}

/*--------------------------------------------------------\
| Fun-03: macroIfMax
| Use:
|  - Set a value (ret) to a value based on which value
|    is greater.
| Input:
|  o ret:
|    - This will hold the return value
|  o x:
|    - First value to compare, (if x >= y)
|  o y:
|    - second value to compare, (if y > x)
|  o xRet:
|    - Value to set ret of x is >= y
|  o yRet:
|    - Value to set ret of y is > x
| Output:
|  - Sets:
|    - ret to xRet if x is greater than or equal to y
|    - ret to yRet if y is greater than x
\--------------------------------------------------------*/
#define macroIfMax(ret, x, y, xRet, yRet){\
   (ret) = (xRet) ^ (((xRet) ^ (yRet)) & (-((x) < (y))));\
   /*This follows the same logic as macroMax(ret, x, y), except
   ` instead of setting ret to the highest value, I set
   ` ret to xRet if x is >= to y or yRet if y > x.
   */ \
}

#endif