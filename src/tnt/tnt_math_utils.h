#ifndef MATH_UTILS_H
#define MATH_UTILS_H


#include <cmath> // std::abs(), std::sqrt()
#include <algorithm> // std::min(), std=::max()
#include <limits> //std::numerical_limits()

namespace TNT
{


/**
 * @brief Computes the hypotenuse of real (non-complex) scalars a and b by 
	avoiding underflow/overflow
	using:
	
	$a * sqrt( 1 + (b/a) * (b/a)))$
	
	rather than
	
	$sqrt(a*a + b*b)$
	
 * @tparam T Input data type, must be the "real" part of a number  
 * @param a first number
 * @param b second number 
 * @return hypotenuse of type T 
 */
template <class T>
T hypot(const T &a, const T &b)
{
	if (a== 0)
		return std::abs(b);
	else
	{
		T c = b/a;
		return std::abs(a) * std::sqrt(1 + c*c);
	}
}
} /* TNT namespace */



#endif
/* MATH_UTILS_H */
