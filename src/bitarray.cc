#include <string>
#include <exception>
#include <iostream>

#include "bitarray.h"

void bit_array::debug () const
{
	for (unsigned bit = 0; bit < m_n_bits; bit++) {
		uint64_t val = m_bits[bit / 64];
		printf ("%d", (int)((val >> (bit % 64)) & 1));
	}
	putchar ('\n');
}

void bit_array::debug (int linesz) const
{
	int linecnt = 0;
	for (unsigned bit = 0; bit < m_n_bits; bit++) {
		uint64_t val = m_bits[bit / 64];
		printf ("%d", (int)((val >> (bit % 64)) & 1));
		if (++linecnt == linesz)
			linecnt = 0, putchar ('\n');
	}
	if (m_n_bits % linesz != 0)
		putchar ('\n');
}
