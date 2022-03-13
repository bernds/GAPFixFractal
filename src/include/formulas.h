#ifndef FORMULAS_H
#define FORMULAS_H

// "Testing" is a placeholder, it can be any temporary formula that's being
// tried out without affecting the order of the others and thereby saves

enum class formula
{
	testing, standard, lambda, tricorn, spider, ship, mix, sqtwice_a, sqtwice_b, celtic, magnet_a,
		facing, facing_b, rings, e90_mix, tails
};

extern const formula formula_table[];

inline bool formula_supports_hybrid (formula f)
{
	return f == formula::tricorn || f == formula::ship || f == formula::celtic;
}

inline bool formula_supports_dem (formula f)
{
	return f == formula::standard || f == formula::lambda || f == formula::mix;
}

inline int n_formula_real_vals (formula, bool dem, bool incolor = false)
{
	int base = dem ? 6 : 4;
	if (incolor)
		return base + 1;
	return base;
}

inline int n_formula_int_vals (formula, bool /* dem */, bool incolor = false)
{
	// Default is only one value in the array: ZPIDX, the index into the zprev array.
	if (incolor)
		// Also return the atom domain period, and keep track of the iteration count.
		return 3;
	return 1;
}

inline int n_formula_extra_doubles (formula, bool /* dem */, bool incolor = false)
{
	// Extra data; the absolute value of the minimum point of the orbit.
	// @@@ Not actually used just yet.
	if (incolor)
		return 1;
	return 0;
}

inline int formula_scratch_space (formula f, int nwords)
{
	if (f == formula::magnet_a || f == formula::facing || f == formula::facing_b || f == formula::rings || f == formula::tails)
		return nwords * 4 + 4;
	return 0;
}

inline bool formula_test_fixpoint (formula)
{
	return false;
}

// Used for shortcuts, can conservatively return false for formulas with zero critpoints.
inline bool formula_zero_critpoint (formula f)
{
	return f == formula::standard || f == formula::tricorn || f == formula::ship || f == formula::celtic;
}

#endif
