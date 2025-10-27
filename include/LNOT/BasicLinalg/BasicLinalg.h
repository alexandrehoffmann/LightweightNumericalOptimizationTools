#ifndef LNOT_BASIC_LINALG_H
#define LNOT_BASIC_LINALG_H

#include <stdbool.h>

#ifdef __cplusplus
	#define LNOT_RESTRICT /* restrict */
#else
	#define LNOT_RESTRICT restrict
#endif // __cplusplus


#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

typedef unsigned int lnot_Size;

#define LNOT_DEFINE_TRIDIAG_LDLT_COMPUTE(Scalar, Suffix)\
	bool lnot_tridiag_ldlt_compute_##Suffix(\
		const Scalar* LNOT_RESTRICT alpha, \
		const Scalar* LNOT_RESTRICT beta, \
		const lnot_Size             size, \
		const Scalar                shift, \
		const Scalar                epsilon, \
		Scalar* LNOT_RESTRICT       invDelta, \
		Scalar* LNOT_RESTRICT       l);\
		\

LNOT_DEFINE_TRIDIAG_LDLT_COMPUTE(float, f)
LNOT_DEFINE_TRIDIAG_LDLT_COMPUTE(double, d)
LNOT_DEFINE_TRIDIAG_LDLT_COMPUTE(long double, ld)

#define LNOT_DEFINE_TRIDIAG_LDLT_SOLVE_LOWER_UNIT(Scalar, Suffix)\
	void lnot_tridiag_ldlt_solveLowerUnit_##Suffix(\
		const Scalar* LNOT_RESTRICT l, \
		const lnot_Size             size, \
		const Scalar                b1, \
		Scalar* LNOT_RESTRICT       x);\
		\

LNOT_DEFINE_TRIDIAG_LDLT_SOLVE_LOWER_UNIT(float, f)
LNOT_DEFINE_TRIDIAG_LDLT_SOLVE_LOWER_UNIT(double, d)
LNOT_DEFINE_TRIDIAG_LDLT_SOLVE_LOWER_UNIT(long double, ld)

#define LNOT_DEFINE_TRIDIAG_LDLT_SOLVE_INPLACE_UPPER(Scalar, Suffix)\
	void lnot_tridiag_ldlt_solveInplaceUpper_##Suffix(\
		const Scalar* LNOT_RESTRICT l, \
		const lnot_Size             size, \
		Scalar* LNOT_RESTRICT       x);\
		\
		
LNOT_DEFINE_TRIDIAG_LDLT_SOLVE_INPLACE_UPPER(float, f)
LNOT_DEFINE_TRIDIAG_LDLT_SOLVE_INPLACE_UPPER(double, d)
LNOT_DEFINE_TRIDIAG_LDLT_SOLVE_INPLACE_UPPER(long double, ld)

#define LNOT_DEFINE_TRIDIAG_LDLT_SOLVE_UNIT(Scalar, Suffix)\
	void lnot_tridiag_ldlt_solveUnit_##Suffix(\
		const Scalar* LNOT_RESTRICT invD, \
		const Scalar* LNOT_RESTRICT l, \
		const lnot_Size             size, \
		const Scalar                b1, \
		Scalar* LNOT_RESTRICT       x);\
		\

LNOT_DEFINE_TRIDIAG_LDLT_SOLVE_UNIT(float, f)
LNOT_DEFINE_TRIDIAG_LDLT_SOLVE_UNIT(double, d)
LNOT_DEFINE_TRIDIAG_LDLT_SOLVE_UNIT(long double, ld)

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // LNOT_BASIC_LINALG_H
