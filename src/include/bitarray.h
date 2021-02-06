#ifndef BITARRAY_H
#define BITARRAY_H

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <algorithm>

inline unsigned popcnt (uint64_t val)
{
	unsigned cnt = 0;
	while (val != 0) {
		val &= val - 1;
		cnt++;
	}
	return cnt;
}

class bit_array
{
	unsigned m_n_bits;
	int m_n_elts;
	uint64_t *m_bits;
	uint64_t m_last_mask;

	void calc_mask (unsigned sz)
	{
		uint64_t mask = 0;
		mask = ~mask;
		mask >>= 63 - (sz + 63) % 64;
		m_last_mask = mask;
	}
public:
	bit_array (unsigned sz, bool set = false)
		: m_n_bits (sz), m_n_elts ((sz + 63) / 64), m_bits (new uint64_t[m_n_elts]())
	{
		calc_mask (sz);
		if (set && m_n_elts > 0) {
			memset (m_bits, -1, m_n_elts * sizeof (uint64_t));
			m_bits[m_n_elts - 1] &= m_last_mask;
		}
	}
	bit_array (const bit_array &other)
		: m_n_bits (other.m_n_bits), m_n_elts (other.m_n_elts), m_bits (new uint64_t[m_n_elts]), m_last_mask (other.m_last_mask)
	{
		memcpy (m_bits, other.m_bits, m_n_elts * sizeof (uint64_t));
	}
	bit_array (bit_array &&other) noexcept
		: m_n_bits (other.m_n_bits), m_n_elts (other.m_n_elts), m_bits (other.m_bits), m_last_mask (other.m_last_mask)
	{
		if (&other == this)
			return;
		other.m_n_bits = 0;
		other.m_n_elts = 0;
		other.m_bits = 0;
	}
	bit_array () : bit_array (1)
	{
	}
	~bit_array ()
	{
		delete[] m_bits;
	}
	bit_array &operator= (bit_array other)
	{
		std::swap (m_n_elts, other.m_n_elts);
		std::swap (m_n_bits, other.m_n_bits);
		std::swap (m_bits, other.m_bits);
		std::swap (m_last_mask, other.m_last_mask);
		return *this;
	}
	bool operator== (const bit_array &other) const
	{
		if (other.m_n_bits != m_n_bits)
			return false;
		return memcmp (m_bits, other.m_bits, m_n_elts * sizeof (uint64_t)) == 0;
	}
	bool operator!= (const bit_array &other) const
	{
		return !(*this == other);
	}
	unsigned bitsize () const { return m_n_bits; }
	uint64_t *raw_bits () const { return m_bits; }
	unsigned raw_n_elts () const { return m_n_elts; }
	void grow (unsigned sz)
	{
		if (sz < m_n_bits)
			throw std::logic_error ("grow called with smaller size");
		if (sz == m_n_bits)
			return;
		int new_nelts = (sz + 63) / 64;
		uint64_t *new_bits = new uint64_t[new_nelts]();
		memcpy (new_bits, m_bits, m_n_elts * sizeof (uint64_t));
		calc_mask (sz);
		m_n_bits = sz;
		m_n_elts = new_nelts;
		std::swap (m_bits, new_bits);
		delete[] new_bits;
	}
	void debug () const;
	void debug (int linesz) const;
	void clear ()
	{
		for (int i = 0; i < m_n_elts; i++)
			m_bits[i] = 0;
	}

	bool test_bit (unsigned n) const
	{
		if (n >= m_n_bits)
			return false;
		unsigned off = n / 64;
		unsigned pos = n % 64;
		uint64_t val = m_bits[off];
		val >>= pos;
		return (val & 1) != 0;
	}

	void set_bit (unsigned n)
	{
		if (n >= m_n_bits)
			return;
		unsigned off = n / 64;
		unsigned pos = n % 64;
		uint64_t val = m_bits[off];
		val |= (uint64_t)1 << pos;
		m_bits[off] = val;
	}

	void clear_bit (unsigned n)
	{
		if (n >= m_n_bits)
			return;
		unsigned off = n / 64;
		unsigned pos = n % 64;
		uint64_t val = m_bits[off];
		val &= ~((uint64_t)1 << pos);
		m_bits[off] = val;
	}

	/* IOR with the other bitmap, shifted left (or right, for negative numbers.
	   Apply a shifted mask before the IOR.  */
	bool ior (const bit_array &other, int shift, const bit_array &mask)
	{
		if (shift == 0)
			return ior (other);
		shift = -shift;

		// int sgn = shift / std::abs (shift);

		/* Examples:
		   shift -3: right shift by 3.
		   Dest word 0 receives [word1:0-3 :: word0:3-64]
		   shift 3:  left shift by 3.
		   Dest word 0 receives [word0:0-61 :: 000 ]
		   Dest word 1 receives [word1:0-61 :: word0:62-64 ]
		   shift 64: left shift by 64.
		   Dest word 0 receives 64 bits from src word -1.  */

		shift += other.m_n_elts * 64;
		int wordshift = (shift + 63) / 64 - other.m_n_elts - 1;
		int bitshift = shift % 64;

		uint64_t last = 0;
		if (wordshift >= 0 && wordshift < other.m_n_elts)
			last = other.m_bits[wordshift];
		if (wordshift >= 0 && wordshift < mask.m_n_elts)
			last &= mask.m_bits[wordshift];
		bool changed = false;
		for (int i = 0; i < m_n_elts; i++) {
			wordshift++;
			uint64_t curr = 0;
			if (wordshift >= 0 && wordshift < other.m_n_elts)
				curr = other.m_bits[wordshift];
			if (wordshift >= 0 && wordshift < mask.m_n_elts)
				curr &= mask.m_bits[wordshift];
			uint64_t val = m_bits[i];
			uint64_t val1 = val;
			if (bitshift != 0) {
				val |= curr << (64 - bitshift);
				val |= last >> bitshift;
			} else {
				val |= curr;
			}
			if (i + 1 == m_n_elts)
				val &= m_last_mask;
			changed |= val != val1;
			m_bits[i] = val;
			last = curr;
		}
		return changed;
	}
	bool ior (const bit_array &other, int shift)
	{
		return ior (other, shift, other);
	}

	bool ior (const bit_array &other)
	{
		int limit = std::min (m_n_elts, other.m_n_elts);
		bool changed = false;
		for (int i = 0; i < limit; i++) {
			uint64_t val = m_bits[i];
			uint64_t val1 = val;
			val |= other.m_bits[i];
			if (i + 1 == m_n_elts)
				val &= m_last_mask;
			changed |= val != val1;
			m_bits[i] = val;
		}
		return changed;
	}
	bool and1 (const bit_array &other)
	{
		int limit = std::min (m_n_elts, other.m_n_elts);
		bool changed = false;
		for (int i = 0; i < limit; i++) {
			uint64_t val = m_bits[i];
			uint64_t val1 = val;
			val &= other.m_bits[i];
			changed |= val != val1;
			m_bits[i] = val;
		}
		return changed;
	}
	bool andnot (const bit_array &other)
	{
		unsigned limit = std::min (m_n_elts, other.m_n_elts);
		bool changed = false;
		for (unsigned i = 0; i < limit; i++) {
			uint64_t val = m_bits[i];
			uint64_t val1 = val;
			val &= ~other.m_bits[i];
			changed |= val != val1;
			m_bits[i] = val;
		}
		return changed;
	}
	bool intersect_p (const bit_array &other) const
	{
		unsigned limit = std::min (m_n_elts, other.m_n_elts);
		for (unsigned i = 0; i < limit; i++) {
			if (m_bits[i] & other.m_bits[i])
				return true;
		}
		return false;
	}
	bool subset_of (const bit_array &other) const
	{
		int limit = std::min (m_n_elts, other.m_n_elts);
		int i;
		for (i = 0; i < limit; i++) {
			if (m_bits[i] & ~other.m_bits[i])
				return false;
		}
		for (; i < m_n_elts; i++)
			if (m_bits[i] != 0)
				return false;
		return true;
	}

	unsigned popcnt () const
	{
		unsigned cnt = 0;
		for (int i = 0; i < m_n_elts; i++)
			cnt += ::popcnt (m_bits[i]);

		return cnt;
	}
	unsigned ffs (int test = 0) const
	{
		int elt = test / 64;
		int bitpos = test % 64;
		uint64_t mask = (uint64_t)1 << bitpos;
		while (elt < m_n_elts) {
			uint64_t v = m_bits[elt];
			while (mask != 0) {
				if (v & mask)
					return test;
				mask <<= 1;
				test++;
			}
			mask = 1;
			elt++;
		}
		return m_n_bits;
	}
	unsigned ffz (int test = 0) const
	{
		int elt = test / 64;
		int bitpos = test % 64;
		uint64_t mask = (uint64_t)1 << bitpos;
		while (elt < m_n_elts) {
			uint64_t v = m_bits[elt];
			while (mask != 0) {
				if ((v & mask) == 0)
					return test;
				mask <<= 1;
				test++;
			}
			mask = 1;
			elt++;
		}
		return m_n_bits;
	}
	/* Extract SIZE bits (at most 64) from OFF.  */
	uint64_t extract (unsigned off, unsigned size) const
	{
		int elt = off / 64;
		int bitpos = off % 64;
		if (elt >= m_n_elts)
			return 0;
		uint64_t v = m_bits[elt] >> bitpos;
		if (bitpos + size > 64 && elt + 1 < m_n_elts) {
			v |= m_bits[elt + 1] << (64 - bitpos);
		}
		uint64_t mask = 0;
		mask = ~mask;
		mask >>= 63 - (size + 63) % 64;
		return v & mask;
	}
	/* Insert SIZE bits (at most 64) at OFF.  */
	void insert (unsigned off, unsigned size, uint64_t data)
	{
		int elt = off / 64;
		int bitpos = off % 64;
		if (elt >= m_n_elts)
			return;
		uint64_t mask = ~(uint64_t)0 << bitpos;
		m_bits[elt] &= ~mask;
		m_bits[elt] |= (data << bitpos) & mask;
		if (bitpos + size > 64 && elt + 1 < m_n_elts) {
			uint64_t mask2 = (uint64_t)1 << (size + bitpos - 64);
			mask2--;
			m_bits[elt + 1] &= ~mask2;
			m_bits[elt + 1] |= (data >> (64 - bitpos)) & mask2;
		}
	}
};

#endif

