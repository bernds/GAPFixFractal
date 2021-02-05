#include <cstdint>
#include <cassert>
#include <cmath>

#include "fpvec.h"

/* These are some crummy versions operating on the fixed-point values we use on
   the GPU side.  These are just good enough to work with coordinates.  */

double to_double (const vpvec &v)
{
	size_t n = v.size ();
	double v1 = (int32_t)v[n - 1];
	double v2 = v[n - 2];
	unsigned long long div = 1;
	div <<= 32;
	double result = v1;
	double ddiv = div;
	for (int i = 1; i < 10; i++) {
		double tmp = v[n - 1 - i];
		result += tmp / ddiv;
		ddiv *= div;
	}
	return result;
}

void set (vpvec &v, double d)
{
	size_t sz = v.size ();
	vpvec val (sz);
	val[sz - 1] = floor (d);
	while (sz-- > 1) {
		double newv = (d - floor (d)) * 2 * (1 << 31);
		val[sz - 1] = floor (newv);
	}
	v = val;
}

vpvec add (const vpvec &srca, const vpvec &srcb)
{
	vpvec res;
	assert (srca.size () == srcb.size ());
	res.reserve (srca.size ());
	uint32_t carry = 0;
	for (size_t i = 0; i < srca.size (); i++) {
		uint64_t r = (uint64_t)srca[i] + (uint64_t)srcb[i] + carry;
		res.push_back (r);
		carry = r >> 32;
	}
	return res;
}

vpvec sub (const vpvec &srca, const vpvec &srcb)
{
	vpvec res;
	assert (srca.size () == srcb.size ());
	res.reserve (srca.size ());
	uint64_t borrow = 0;
	for (size_t i = 0; i < srca.size (); i++) {
		uint64_t r = (uint64_t)srca[i] - (uint64_t)srcb[i] + borrow;
		res.push_back (r);
		borrow = (int64_t)r >> 32;
	}
	return res;
}

vpvec div1 (const vpvec &src, int32_t d)
{
	size_t len = src.size ();
	uint32_t val = (uint32_t)0 - (src[len - 1] >> 31);
	vpvec result (len);
	for (size_t i = len; i-- > 0;) {
		uint64_t v = ((uint64_t)val << 32) | src[i];
		int64_t iv = v;
		int64_t q = iv / d;
		val = v - q * d;
		result[i] = q;
	}
	return result;
}

vpvec mul1 (const vpvec &src, uint32_t d)
{
	size_t len = src.size ();
	vpvec result;
	uint64_t v = 0;
	for (size_t i = 0; i < len; i++) {
		v += (uint64_t)src[i] * d;
		result.push_back (v);
		v >>= 32;
	}
	return result;
}
