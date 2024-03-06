#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fq.h>
#include <flint/fq_poly.h>

#undef WORD
#undef CAT
#undef _CAT

#undef DIGIT
#define DIGIT 64

#define CAT(A, B)			_CAT(A, B)
#define _CAT(A, B)			A ## B

#include <relic/relic.h>

#if (FP_PRIME % 64) > 0
#define RELIC_DIGS	(FP_PRIME/64 + 1)
#else
#define RELIC_DIGS	(FP_PRIME/64)
#endif


#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>

/*
#undef WORD
#define WORD 64
*/

/*! \file bn_ext.h
    \brief An extension of the relic integer-handling functions
    
	These functions handle operations on integers or arrays of integers,
	as well as converstions and macros. The integers are all considered of size
	RELIC_DIGS in digits and reduced to the order of the elliptic curve groups.

    \ingroup utils
*/


// Defines a global mult_type constant which determines the kind of multiplication used:
// 		b -> basic
// 		k -> karatsuba
//		c -> comba

#ifndef MULT_TYPE
#define MULT_TYPE 'b'
#endif

#ifndef MOD_TYPE
#define MOD_TYPE 'b'
#endif


void bn_rand_mod(bn_t a, bn_t mod);

/**
 * Converts an mpz (gmp) integer to a bn (relic) integer 
 * @param[out] out			bn integer
 * @param[in] in 			mpz integer
 */
int mpz2bn(bn_t out, mpz_t in);

/**
 * Converts an fmpz (flint) integer to a bn (relic) integer 
 * @param[out] out			bn integer
 * @param[in] in 			mpz integer
 */
int fmpz2bn(bn_t out, fmpz_t in);

/**
 * Converts a bn (relic) integer to an mpz (gmp integer) 
 * @param[out] out			mpz integer
 * @param[in] in 			bn integer
 */
int bn2mpz(mpz_t out, bn_t in);

/**
 * Converts a bn (relic) integer to an fmpz (flint integer)
 * @param[out] out			mpz integer
 * @param[in] in 			bn integer
 */
int bn2fmpz(fmpz_t out, bn_t in);

/**
 * Converts a string of 4 unsigned bytes to an integer.
 * @param[out] out 			Output integer.
 * @param[in] in 	 		Pointer to the byte array to convert.
 */
int uint8_t2int(int * out, uint8_t * in);

/**
 * Converts a bn integer to an int.
 * @param[out] out 			Output bn integer.
 * @param[in] in 	 		Signed integer
 * @param[in] size 	 		The size of the integer in base 10.
 */
int bn2int(int * out, bn_t in);



