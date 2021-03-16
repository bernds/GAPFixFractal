#ifndef FORMULAS_H
#define FORMULAS_H

// "Testing" is a placeholder, it can be any temporary formula that's being
// tried out without affecting the order of the others and thereby saves

enum class formula
{
	testing, standard, lambda, tricorn, spider, ship, mix, sqtwice_a, sqtwice_b, celtic, magnet_a,
		facing, facing_b, rings
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

/* The number of real values that the kernel needs, and which must be preserved across
   invocations.  */
inline int n_formula_real_vals (formula, bool dem)
{
	/* Z, C (as computed from the coords), and possibly ZDER.  */
	return dem ? 6 : 4;
}

/* The number of integer values that the kernel needs, and which must be preserved across
   invocations.  */
inline int n_formula_int_vals (formula, bool /* dem */)
{
	/* ZPIDX - the index into the zprev array.  */
	return 1;
}

inline int formula_scratch_space (formula f, int nwords)
{
	if (f == formula::magnet_a || f == formula::facing || f == formula::facing_b || f == formula::rings)
		return nwords * 4 + 4;
	return 0;
}

inline bool formula_test_fixpoint (formula)
{
	return false;
}
#endif
