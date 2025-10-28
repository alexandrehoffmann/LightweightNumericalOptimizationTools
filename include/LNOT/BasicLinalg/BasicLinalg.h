#ifndef LNOT_BASIC_LINALG_H
#define LNOT_BASIC_LINALG_H

#include <stdbool.h>
#include <math.h>

#ifdef __cplusplus
	#define LNOT_RESTRICT /* restrict */
#else
	#define LNOT_RESTRICT restrict
#endif // __cplusplus


#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

typedef enum { LNOT_MAT_ROW_MAJOR, LNOT_MAT_COL_MAJOR } lnot_mat_StorageOrder; ///<  @brief Enum for matrix storage order.
typedef enum { LNOT_MAT_LOWER,     LNOT_MAT_UPPER }     lnot_mat_UpLo;         ///<  @brief Enum that specifies whether the upper or lower triangular part of 

typedef unsigned int lnot_Size;

inline float       lnot_abs_f (const float       v) { return fabsf(v); } 
inline double      lnot_abs_d (const double      v) { return fabs (v); } 
inline long double lnot_abs_ld(const long double v) { return fabsl(v); }

inline float       lnot_max_f (const float       a, const float       b) { return fmaxf(a, b); } 
inline double      lnot_max_d (const double      a, const double      b) { return fmax (a, b); } 
inline long double lnot_max_ld(const long double a, const long double b) { return fmaxl(a, b); } 

#define LNOT_DEFINE_SYM_MATRIX_VECTOR_PROD(Scalar, Suffix)\
	void lnot_symMatrixVectorProd_##Suffix(\
		const lnot_mat_StorageOrder layout, \
		const lnot_mat_UpLo uplo, \
		const Scalar alpha, \
		const Scalar* LNOT_RESTRICT A, \
		const Scalar* LNOT_RESTRICT x, \
		const lnot_Size N, \
		Scalar* LNOT_RESTRICT y);\
		\

LNOT_DEFINE_SYM_MATRIX_VECTOR_PROD(float, f)
LNOT_DEFINE_SYM_MATRIX_VECTOR_PROD(double, d)
LNOT_DEFINE_SYM_MATRIX_VECTOR_PROD(long double, ld)

#define LNOT_DEFINE_SYM_RK1_UPDATE(Scalar, Suffix)\
	void lnot_symRk1Update_##Suffix(\
		const lnot_mat_StorageOrder layout, \
		const lnot_mat_UpLo         uplo, \
		const Scalar                alpha, \
		const Scalar* LNOT_RESTRICT x, \
		const lnot_Size             N, \
		Scalar* LNOT_RESTRICT       A);\
		\

LNOT_DEFINE_SYM_RK1_UPDATE(float, f)
LNOT_DEFINE_SYM_RK1_UPDATE(double, d)
LNOT_DEFINE_SYM_RK1_UPDATE(long double, ld)

#define LNOT_DEFINE_SYM_RK2_UPDATE(Scalar, Suffix)\
	void lnot_symRk2Update_##Suffix(\
		const lnot_mat_StorageOrder layout, \
		const lnot_mat_UpLo         uplo, \
		const Scalar                alpha, \
		const Scalar* LNOT_RESTRICT x, \
		const Scalar* LNOT_RESTRICT y, \
		const lnot_Size             N, \
		Scalar* LNOT_RESTRICT       A);\
		\

LNOT_DEFINE_SYM_RK2_UPDATE(float, f)
LNOT_DEFINE_SYM_RK2_UPDATE(double, d)
LNOT_DEFINE_SYM_RK2_UPDATE(long double, ld)

#define LNOT_DEFINE_TRIDIAG_NORM1(Scalar, Suffix)\
	Scalar lnot_tridiag_norm1_##Suffix(\
		const Scalar* LNOT_RESTRICT alpha, \
		const Scalar* LNOT_RESTRICT beta, \
		const lnot_Size N);\
	\

LNOT_DEFINE_TRIDIAG_NORM1(float, f)
LNOT_DEFINE_TRIDIAG_NORM1(double, d)
LNOT_DEFINE_TRIDIAG_NORM1(long double, ld)

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
