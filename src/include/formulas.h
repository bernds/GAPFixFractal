#ifndef FORMULAS_H
#define FORMULAS_H

// "Testing" is a placeholder, it can be any temporary formula that's being
// tried out without affecting the order of the others and thereby saves

enum class formula
{
	testing, standard, lambda, tricorn, spider, ship, mix, sqtwice_a, sqtwice_b, celtic, magnet_a,
	facing
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

inline int n_formula_cplx_vals (formula, bool dem)
{
	return dem ? 3 : 2;
}

inline int formula_scratch_space (formula f, int nwords)
{
	if (f == formula::magnet_a || f == formula::facing)
		return nwords * 4 + 4;
	return 0;
}

inline bool formula_test_fixpoint (formula)
{
	return false;
}
#endif
