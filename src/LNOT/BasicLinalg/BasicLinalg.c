#include <LNOT/BasicLinalg/BasicLinalg.h>

#define LNOT_IMPLEMENT_SYM_MATRIX_VECTOR_PROD(Scalar, Suffix)\
	void lnot_symMatrixVectorProd_##Suffix(\
		const lnot_mat_StorageOrder layout, \
		const lnot_mat_UpLo uplo, \
		const Scalar alpha, \
		const Scalar* LNOT_RESTRICT A, \
		const Scalar* LNOT_RESTRICT x, \
		const lnot_Size N, \
		Scalar* LNOT_RESTRICT y)\
	{\
		const lnot_Size iStride = layout == LNOT_MAT_ROW_MAJOR && uplo == LNOT_MAT_LOWER ? N : 1;\
		const lnot_Size jStride = layout == LNOT_MAT_ROW_MAJOR && uplo == LNOT_MAT_LOWER ? 1 : N; \
	\
		for (lnot_Size i=0; i!=N; ++i)\
		{\
			for (lnot_Size j=0; j!=(i+1); ++j)\
			{\
				y[i] += alpha*A[i*iStride + j*jStride]*x[j];\
			}\
			for (lnot_Size j=i+1; j!=N; ++j)\
			{\
				y[i] += alpha*A[j*iStride + i*jStride]*x[j];\
			}\
		}\
	}\
	\

LNOT_IMPLEMENT_SYM_MATRIX_VECTOR_PROD(float, f)
LNOT_IMPLEMENT_SYM_MATRIX_VECTOR_PROD(double, d)
LNOT_IMPLEMENT_SYM_MATRIX_VECTOR_PROD(long double, ld)

#define LNOT_IMPLEMENT_SYM_RK1_UPDATE(Scalar, Suffix)\
	void lnot_symRk1Update_##Suffix(\
		const lnot_mat_StorageOrder layout, \
		const lnot_mat_UpLo         uplo, \
		const Scalar                alpha, \
		const Scalar* LNOT_RESTRICT x, \
		const lnot_Size             N, \
		Scalar* LNOT_RESTRICT       A)\
	{\
		const lnot_Size iStride = layout == LNOT_MAT_ROW_MAJOR && uplo == LNOT_MAT_LOWER ? N : 1;\
		const lnot_Size jStride = layout == LNOT_MAT_ROW_MAJOR && uplo == LNOT_MAT_LOWER ? 1 : N; \
	\
		for (lnot_Size i=0; i!=N; ++i)\
		{\
			const Scalar alpha_xi = alpha*x[i];\
			\
			for (lnot_Size j=0; j!=(i+1); ++j)\
			{\
				A[i*iStride + j*jStride] += alpha_xi*x[j];\
			}\
		}\
	}\
	\

LNOT_IMPLEMENT_SYM_RK1_UPDATE(float, f)
LNOT_IMPLEMENT_SYM_RK1_UPDATE(double, d)
LNOT_IMPLEMENT_SYM_RK1_UPDATE(long double, ld)

#define LNOT_IMPLEMENT_SYM_RK2_UPDATE(Scalar, Suffix)\
	void lnot_symRk2Update_##Suffix(\
		const lnot_mat_StorageOrder layout, \
		const lnot_mat_UpLo         uplo, \
		const Scalar                alpha, \
		const Scalar* LNOT_RESTRICT x, \
		const Scalar* LNOT_RESTRICT y, \
		const lnot_Size             N, \
		Scalar* LNOT_RESTRICT       A)\
	{\
		const lnot_Size iStride = layout == LNOT_MAT_ROW_MAJOR && uplo == LNOT_MAT_LOWER ? N : 1;\
		const lnot_Size jStride = layout == LNOT_MAT_ROW_MAJOR && uplo == LNOT_MAT_LOWER ? 1 : N; \
	\
		for (lnot_Size i=0; i!=N; ++i)\
		{\
			const Scalar alpha_xi = alpha*x[i];\
			const Scalar alpha_yi = alpha*y[i];\
			\
			for (lnot_Size j=0; j!=(i+1); ++j)\
			{\
				A[i*iStride + j*jStride] += alpha_xi*y[j] + alpha_yi*x[j];\
			}\
		}\
	}\
	\

LNOT_IMPLEMENT_SYM_RK2_UPDATE(float, f)
LNOT_IMPLEMENT_SYM_RK2_UPDATE(double, d)
LNOT_IMPLEMENT_SYM_RK2_UPDATE(long double, ld)

#define LNOT_IMPLEMENT_TRIDIAG_LDLT_COMPUTE(Scalar, Suffix)\
	bool lnot_tridiag_ldlt_compute_##Suffix(\
		const Scalar* restrict alpha, \
		const Scalar* restrict beta, \
		const lnot_Size size, \
		const Scalar shift, \
		const Scalar epsilon, \
		Scalar* restrict invDelta, \
		Scalar* restrict l)\
	{\
		if ((alpha[0] + shift) < epsilon) { return false; }\
	\
		invDelta[0] = 1 / (alpha[0] + shift);\
		if (size == 1) { return true; }\
		l[0] = beta[0]*invDelta[0];\
	\
		for (lnot_Size i=1; i!=(lnot_Size)(size-1);++i)\
		{\
			const Scalar delta_i = alpha[i] + shift - beta[i-1]*l[i-1];\
			if (delta_i < epsilon) { return false; }\
			invDelta[i] = 1 / delta_i;\
			l[i] = beta[i]*invDelta[i];\
		}\
		{\
			const lnot_Size i = (lnot_Size)(size-1);\
			const Scalar delta_i = alpha[i] + shift - beta[i-1]*l[i-1];\
			if (delta_i < epsilon) { return false; }\
			invDelta[i] = 1 / delta_i;\
		}\
	\
		return true;\
	}\
	\

LNOT_IMPLEMENT_TRIDIAG_LDLT_COMPUTE(float, f)
LNOT_IMPLEMENT_TRIDIAG_LDLT_COMPUTE(double, d)
LNOT_IMPLEMENT_TRIDIAG_LDLT_COMPUTE(long double, ld)

#define LNOT_IMPLEMENT_TRIDIAG_LDLT_SOLVE_LOWER_UNIT(Scalar, Suffix)\
	void lnot_tridiag_ldlt_solveLowerUnit_##Suffix(\
		const Scalar* restrict l, \
		const lnot_Size size, \
		const Scalar b1, \
		Scalar* restrict x) \
	{\
		x[0] = b1;\
		for (lnot_Size i=1; i!=size; ++i) { x[i] = -l[i-1]*x[i-1]; }\
	}\
	\

LNOT_IMPLEMENT_TRIDIAG_LDLT_SOLVE_LOWER_UNIT(float, f)
LNOT_IMPLEMENT_TRIDIAG_LDLT_SOLVE_LOWER_UNIT(double, d)
LNOT_IMPLEMENT_TRIDIAG_LDLT_SOLVE_LOWER_UNIT(long double, ld)

#define LNOT_IMPLEMENT_TRIDIAG_LDLT_SOLVE_INPLACE_UPPER(Scalar, Suffix)\
	void lnot_tridiag_ldlt_solveInplaceUpper_##Suffix(\
		const Scalar* restrict l, \
		const lnot_Size size, \
		Scalar* restrict x)\
	{\
		for (lnot_Size i=(lnot_Size)(size-2); i!=(lnot_Size)(-1); --i) { x[i] -= l[i]*x[i+1]; }\
	}\
	\

LNOT_IMPLEMENT_TRIDIAG_LDLT_SOLVE_INPLACE_UPPER(float, f)
LNOT_IMPLEMENT_TRIDIAG_LDLT_SOLVE_INPLACE_UPPER(double, d)
LNOT_IMPLEMENT_TRIDIAG_LDLT_SOLVE_INPLACE_UPPER(long double, ld)

#define LNOT_IMPLEMENT_TRIDIAG_LDLT_SOLVE_UNIT(Scalar, Suffix)\
	void lnot_tridiag_ldlt_solveUnit_##Suffix(\
		const Scalar* restrict invD, \
		const Scalar* restrict l, \
		const lnot_Size size, \
		const Scalar b1, \
		Scalar* restrict x)\
	{\
		lnot_tridiag_ldlt_solveLowerUnit_##Suffix(l, size, b1, x);\
		for (lnot_Size i=0; i!=size; ++i) { x[i] *= invD[i]; }\
		lnot_tridiag_ldlt_solveInplaceUpper_##Suffix(l, size, x);\
	}\
	\
	
LNOT_IMPLEMENT_TRIDIAG_LDLT_SOLVE_UNIT(float, f)
LNOT_IMPLEMENT_TRIDIAG_LDLT_SOLVE_UNIT(double, d)
LNOT_IMPLEMENT_TRIDIAG_LDLT_SOLVE_UNIT(long double, ld)
