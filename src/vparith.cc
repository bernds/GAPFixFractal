#include <cstdint>
#include <cassert>
#include <cmath>

#include <QString>

#include "fpvec.h"

/* These are some crummy versions operating on the fixed-point values we use on
   the GPU side.  These are just good enough to work with coordinates.  */

double to_double (const uint32_t *v, int n)
{
	double v1 = (int32_t)v[n - 1];
	unsigned long long div = 1;
	div <<= 32;
	double result = v1;
	double ddiv = div;
	for (int i = 1; i < n; i++) {
		double tmp = v[n - 1 - i];
		result += tmp / ddiv;
		ddiv *= div;
	}
	return result;
}

double to_double (const vpvec &v)
{
	return to_double (&v[0], v.size ());
}

void set (vpvec &v, double d)
{
	size_t sz = v.size ();
	vpvec val (sz);
	while (sz-- > 1) {
		double newv_f = floor (d);
		val[sz] = newv_f;
		d -= newv_f;
		d = d * 4 * (1 << 30);
	}
	v = val;
}

QString to_string (vpvec src, int prec)
{
	int sz = src.size ();
	bool neg = false;
	if (src[sz - 1] & 0x80000000) {
		vpvec zero (sz, 0);
		src = sub (zero, src);
		neg = true;
	}
	QString n = QString::number (src[sz - 1]);
	if (neg)
		n = "-" + n;
	QString frac = ".";
	for (int i = 0; i < (prec - 1) * 11; i++) {
		src[sz - 1] = 0;
		bool found = false;
		for (auto val: src)
			if (val != 0) {
				found = true;
				break;
			}
		if (!found)
			break;
		src = mul1 (src, 10);
		frac += QChar ('0' + src[sz - 1]);
	}
	if (frac != ".")
		n += frac;
	return n;
}

vpvec from_string (QString str, int sz)
{
	if (str.isEmpty ())
		throw invalid_decimal_string ();
	bool neg = str[0] == '-';
	if (neg)
		str.remove (0, 1);
	QChar dpoint = '.';
	QChar sep = ',';
	// Try to accomodate German-style numbers with comma instead of decimal point.
	if (str.count (sep) == 1) {
		if (str.count (dpoint) > 1 || str.count (dpoint) == 0)
			std::swap (dpoint, sep);
	}
	int decimal = str.indexOf (dpoint, 0);
	QString fraction = "0";
	QString nonfrac = str;
	if (decimal != -1) {
		nonfrac = str.section (dpoint, 0, 0);
		fraction = str.section (dpoint, 1);
		fraction.remove (sep);
	}
	bool ok;
	int nonfrac_int = nonfrac.toInt (&ok);
	if (!ok)
		throw invalid_decimal_string ();
	vpvec result (sz, 0);
	while (!fraction.isEmpty ()) {
		QString tmp = fraction.right (8);
		fraction.chop (tmp.size ());
		result[sz - 1] = tmp.toInt (&ok);
		if (!ok)
			throw invalid_decimal_string ();
		int d = 1;
		for (int i = 0; i < tmp.size (); i++)
			d *= 10;
		result = div1 (result, d);
	}
	result[sz - 1] = nonfrac_int;
	if (neg) {
		vpvec zero (sz, 0);
		result = sub (zero, result);
	}
	return result;
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
