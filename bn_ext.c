#include "bn_ext.h"


   void bn_rand_mod(bn_t a, bn_t mod)
   {
        bn_rand(a,BN_MOD,D_BITS);
        bn_mod_basic(a, a, mod);
        return;
   }

int mpz2bn(bn_t out, mpz_t in)		
{
	int size = mpz_sizeinbase(in, 10) + 2;		// + 2 bytes for minus sign and NULL
	char * temp = malloc(sizeof(char)*size);

	mpz_get_str(temp, 10, in);
	bn_zero(out);
	bn_read_str(out, temp, size, 10);

	free(temp);
	return 0;
}

int fmpz2bn(bn_t out, fmpz_t in)		
{
	mpz_t temp;
	mpz_init(temp);

	fmpz_get_mpz(temp , in);

	mpz2bn(out, temp);

	mpz_clear(temp);

	return 0;
}

int bn2mpz(mpz_t out, bn_t in)
{
	int size = bn_size_str(in, 10);
	char * temp = malloc(sizeof(char)*size);

	bn_write_str(temp, size, in, 10);
	mpz_set_str(out,temp, 10);

	free(temp);
	return 0;
}

int bn2fmpz(fmpz_t out, bn_t in)
{
	mpz_t temp;
	mpz_init(temp);

	bn2mpz(temp, in);

	fmpz_set_mpz(out , temp);

	mpz_clear(temp);
	return 0;
}

int uint8_t2int(int * out, uint8_t * in)
{
	// The byte array is in big endian
	int num = 1;
	if (*(char *)&num == 1)			// Little endian
	{
		*out = (int) in[0];
		*out = (*out << 8);

		*out += (int) in[1];
		*out = (*out << 8);

		*out += (int) in[2];
		*out = (*out << 8);

		*out += (int) in[3];
	}
	else							// Big endian
	{
		*out = (int) in[0];
		*out = (*out >> 8);

		*out += (int) in[1];
		*out = (*out >> 8);

		*out += (int) in[2];
		*out = (*out >> 8);

		*out += (int) in[3];
	}
	

    return 0;
}

int bn2int(int * out, bn_t in)
{
	uint8_t * temp = malloc(sizeof(uint8_t)*4);
	bn_write_bin(temp, 4, in);
	uint8_t2int(out, temp);
	free(temp);

	return 0;
}
