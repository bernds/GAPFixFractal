#ifndef FORMULAS_H
#define FORMULAS_H

// "Testing" is a placeholder, it can be any temporary formula that's being
// tried out without affecting the order of the others and thereby saves

enum class formula
{
	testing, standard, lambda, tricorn, spider, ship, mix, sqtwice_a, sqtwice_b, celtic
};

extern const formula formula_table[];

inline bool formula_supports_hybrid (formula f)
{
	return f == formula::tricorn || f == formula::ship || f == formula::celtic;
}

#endif
