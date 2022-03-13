#include <QString>

#include <vector>
#include <cstring>
#include <memory>
#include <array>

#include <iostream>

#include "formulas.h"
#include "genkernel.h"

using std::string;
using std::vector;
using std::array;
using std::unique_ptr;
using std::make_unique;

class expr;

class generator
{
	int m_tmp = 0;
	int m_lbl_tmp = 0;
	vector<unique_ptr<expr>> m_cleanup;

public:
	QString m_code, m_regs;

	QString gen_reg (const QString &type, QString nm = "reg", bool declare = true)
	{
		QString r = QString ("%%1%2").arg (nm).arg (m_tmp++);
		if (declare)
			m_regs += QString ("\t.reg.%1\t%2;\n").arg (type).arg (r);
		return r;
	}
	QString gen_label (QString nm = "lbl")
	{
		QString r = QString ("%1%2").arg (nm).arg (m_lbl_tmp++);
		return r;
	}
	void append_code (const QString &c)
	{
		m_code += c;
	}
	void append_move (const QString &type, const QString &dst, const QString &src)
	{
		m_code += QString ("\tmov.%1\t%2, %3;\n").arg (type, dst, src);
	}
	QString code ()
	{
		QString t = m_regs + m_code;
		m_code = QString ();
		m_regs = QString ();
		return t;
	}
	template<class T, class ... ARGS>
	T *make (ARGS &&... args)
	{
		T *v = new T (std::forward<ARGS>(args)...);
		m_cleanup.emplace_back ((expr *)v);
		return v;
	}
};

static generator *codegen;

template<class T, class ... ARGS>
T *make (ARGS &&... args)
{
	return codegen->make<T> (std::forward<ARGS>(args)...);
}

class expr
{
protected:
	vector<QString> m_values;
	QString m_dstreg = "reg";
	size_t m_len;

public:
	expr (int len) : m_len (len)
	{
	}
	virtual ~expr () { }
	int length ()
	{
		return m_len;
	}
	void set_destreg (const QString &s)
	{
		m_dstreg = s;
	}
	QString get_destreg (QString type)
	{
		return codegen->gen_reg (type, m_dstreg + QString ("_%1_").arg (m_values.size ()));
	}
	void calculate_full ()
	{
		while (m_values.size () < m_len)
			m_values.push_back (next_piece ());
	}
	virtual void require_carry ()
	{
		calculate_full ();
	}
	QString get_piece (size_t v)
	{
		while (v >= m_values.size ())
			m_values.push_back (next_piece ());

		return m_values[v];
	}
	QString get_piece_high (size_t v)
	{
		if (v >= m_len)
			return "0";
		return get_piece (m_len - v - 1);
	}
	virtual uint32_t nonzero_bits (size_t v)
	{
		if (v >= m_len)
			return 0;
		return ~(uint32_t)0;
	}
	uint32_t nonzero_bits_high (size_t v)
	{
		if (v >= m_len)
			return 0;
		return nonzero_bits (m_len - v - 1);
	}
	virtual QString next_piece () = 0;
};

template<int N, uint32_t F=0>
class const_expr : public expr
{
public:
	const_expr (int len) : expr (len)
	{
		while (len-- > 2)
			m_values.push_back ("0");
		m_values.push_back (QString ("%1").arg (F));
		m_values.push_back (QString ("%1").arg (N));
	}
	QString next_piece () override
	{
		abort ();
	}
};

class reg_expr : public expr
{
public:
	reg_expr (int len, QString nm = "reg", bool declare = true) : expr (len)
	{
		while (len-- > 0)
			m_values.push_back (codegen->gen_reg ("u32", nm, declare));
	}
	QString next_piece () override
	{
		abort ();
	}
};

class round_cst_expr : public expr
{
public:
	round_cst_expr (int len, uint32_t val) : expr (len)
	{
		m_values.push_back (QString ("%1").arg (val));
	}
	QString next_piece () override
	{
		return "0";
	}
};

class ldg_expr : public expr
{
	QString m_addr;
public:
	ldg_expr (const QString &addr, int len)
		: expr (len), m_addr (addr)
	{
	}
	void require_carry () override
	{
	}
	QString next_piece () override
	{
		int l = m_values.size ();
		QString r = get_destreg ("u32");
		codegen->append_code (QString ("\tld.global.u32\t%1, [%2 + %3];\n")
				      .arg (r, m_addr).arg (l * 4));
		return r;
	}
};

class ldc_expr : public expr
{
	QString m_addr;
public:
	ldc_expr (const QString &addr, int len)
		: expr (len)
	{
		m_addr = codegen->gen_reg ("u32");
		codegen->append_move ("u32", m_addr, addr);
	}
	void require_carry () override
	{
	}
	QString next_piece () override
	{
		int l = m_values.size ();
		QString r = get_destreg ("u32");
		codegen->append_code (QString ("\tld.const.u32\t%1, [%2 + %3];\n")
				      .arg (r, m_addr).arg (l * 4));
		return r;
	}
};

class lds_expr : public expr
{
	QString m_addr, m_stride;
public:
	lds_expr (const QString &addr, const QString &stride, int len)
		: expr (len), m_addr (codegen->gen_reg ("u32")), m_stride (stride)
	{
		codegen->append_move ("u32", m_addr, addr);
	}
	void require_carry () override
	{
	}
	QString next_piece () override
	{
		QString r = get_destreg ("u32");
		codegen->append_code (QString ("\tld.shared.u32\t%1, [%2];\n")
				      .arg (r, m_addr));
		codegen->append_code (QString ("\tadd.u32\t%1, %1, %2;\n")
				      .arg (m_addr, m_stride));
		return r;
	}
};

class from_float_expr : public expr
{
public:
	from_float_expr (const QString &f64) : expr (2)
	{
		QString r0 = codegen->gen_reg ("u32");
		QString r1 = codegen->gen_reg ("u32");
		QString tmpf = codegen->gen_reg ("f64");
		m_values.push_back (r0);
		m_values.push_back (r1);
		codegen->append_code (QString ("\tcvt.rzi.u32.f64\t%1, %2;\n").arg (r1, f64));
		codegen->append_code (QString ("\tcvt.rn.f64.u32\t%1, %2;\n").arg (tmpf, r1));
		codegen->append_code (QString ("\tsub.f64\t%1, %2, %1;\n").arg (tmpf, f64));
		codegen->append_code (QString ("\tmul.f64\t%1, %1, 4294967296.0;\n").arg (tmpf));
		codegen->append_code (QString ("\tcvt.rzi.u32.f64\t%1, %2;\n").arg (r0, tmpf));
	}
	QString next_piece () override
	{
		abort ();
	}
};

class addsub_expr : public expr
{
	QString m_op, m_cond;
	expr *m_a, *m_b;
	bool nonzero_added = false;
public:
	addsub_expr (const QString &op, expr *srca, expr *srcb, QString cond = QString ()) :
		expr (std::max (srca->length (), srcb->length ())),
		m_op (op), m_cond (cond.isEmpty () ? "" : "@" + cond), m_a (srca), m_b (srcb)
	{
		srca->require_carry ();
		srcb->require_carry ();
	}
	QString next_piece () override
	{
		int idx = m_values.size ();
		int rev = m_len - idx - 1;
		QString a = m_a->get_piece_high (rev);
		QString b = m_b->get_piece_high (rev);
		if (a == "0" && !nonzero_added && m_op == "add")
			return b;
		if (b == "0" && !nonzero_added)
			return a;

		QString dst = get_destreg ("u32");
		codegen->append_code (QString ("%1\t%2%3.cc.u32\t%4,%5,%6;\n")
				      .arg (m_cond, m_op, nonzero_added ? "c" : "", dst, a, b));
		nonzero_added = true;
		return dst;
	}
	uint32_t nonzero_bits (size_t v) override
	{
		if (m_op != "add")
			return expr::nonzero_bits (v);
		uint32_t anz = m_a->nonzero_bits (v);
		uint32_t bnz = m_b->nonzero_bits (v);
		uint32_t result = ~(uint32_t) 0;
		uint32_t testbits = (uint32_t)1 << 30;
		uint32_t test = anz | bnz;
		for (int i = 0; i < 31; i++) {
			if (test & testbits)
				return result;
			result >>= 1;
		}
		if (v == 0)
			return 1;
		anz = m_a->nonzero_bits (v - 1);
		bnz = m_b->nonzero_bits (v - 1);
		uint32_t testbit = (uint32_t)1 << 31;
		if (((anz | bnz) & testbit) == 0)
			return 0;
		return 1;
	}

};

template<int N>
class arshift_expr : public expr
{
	expr *m_op;
public:
	arshift_expr (expr *op) :
		expr (op->length ()), m_op (op)
	{
	}
	void require_carry () override
	{
		m_op->require_carry ();
	}
	QString next_piece () override
	{
		int idx = m_values.size ();
		QString a = m_op->get_piece (idx);
		QString b;
		if (idx + 1 == m_op->length ()) {
			b = codegen->gen_reg ("u32");
			codegen->append_code (QString ("\tshr.s32\t%1, %2, 31;\n").arg (b, a));
		} else
			b = m_op->get_piece (idx + 1);
		QString dst = codegen->gen_reg ("u32");
		codegen->append_code (QString ("\tshf.r.clamp.b32\t%1,%2,%3,%4;\n")
				      .arg (dst, a, b).arg (N));
		return dst;
	}
};

class lshift_expr : public expr
{
	expr *m_op;
	int m_shift;
public:
	lshift_expr (expr *op, int shift) :
		expr (op->length ()), m_op (op), m_shift (shift)
	{
	}
	void require_carry () override
	{
		m_op->require_carry ();
	}
	QString next_piece () override
	{
		int idx = m_values.size ();
		QString a = m_op->get_piece (idx);
		QString b;
		QString dst = codegen->gen_reg ("u32");
		if (idx == 0) {
			codegen->append_code (QString ("\tshl.b32\t%1, %2, %3;\n").arg (dst, a).arg (m_shift));
		} else {
			b = m_op->get_piece (idx - 1);
			codegen->append_code (QString ("\tshf.l.clamp.b32\t%1,%2,%3,%4;\n")
					      .arg (dst, b, a).arg (m_shift));
		}
		return dst;
	}
};

class padlow_expr : public expr
{
	expr *m_op;
public:
	padlow_expr (expr *src, int len)
		: expr (len), m_op (src)
	{
		assert (len > m_op->length ());
	}
	void require_carry () override
	{
		m_op->require_carry ();
	}
	QString next_piece () override
	{
		int idx = m_values.size ();
		if (idx < length () - m_op->length ())
			return "0";
		return m_op->get_piece (idx + m_op->length () - m_len);
	}
	uint32_t nonzero_bits (size_t v) override
	{
		if ((int)v < length () - m_op->length ())
			return 0;
		return m_op->nonzero_bits (v + m_op->length () - m_len);
	}
};

class padhigh_expr : public expr
{
	expr *m_op;
public:
	padhigh_expr (expr *src, int len)
		: expr (len), m_op (src)
	{
		assert (len >= m_op->length ());
	}
	void require_carry () override
	{
		m_op->require_carry ();
	}
	QString next_piece () override
	{
		int idx = m_values.size ();
		if (idx >= m_op->length ())
			return "0";
		return m_op->get_piece (idx);
	}
	uint32_t nonzero_bits (size_t v) override
	{
		if (v >= m_op->length ())
			return 0;
		return m_op->nonzero_bits (v);
	}
};

class trunc_expr : public expr
{
	expr *m_op;
public:
	trunc_expr (expr *src, int len)
		: expr (len), m_op (src)
	{
		assert (len <= m_op->length ());
	}
	void require_carry () override
	{
		m_op->require_carry ();
	}
	QString next_piece () override
	{
		return m_op->get_piece (m_op->length () - m_len + m_values.size ());
	}
	uint32_t nonzero_bits (size_t v) override
	{
		return m_op->nonzero_bits (m_op->length () - m_len + v);
	}
};

class concat_expr : public expr
{
	expr *m_low, *m_high;
public:
	concat_expr (expr *op_low, expr *op_high)
		: expr (op_low->length () + op_high->length ()), m_low (op_low), m_high (op_high)
	{
	}
	void require_carry () override
	{
		m_low->require_carry ();
		m_high->require_carry ();
	}
	QString next_piece () override
	{
		int idx = m_values.size ();
		if (idx < m_low->length ())
			return m_low->get_piece (idx);
		return m_high->get_piece (idx - m_low->length ());
	}
	uint32_t nonzero_bits (size_t v) override
	{
		if (v < m_low->length ())
			return m_low->nonzero_bits (v);
		return m_high->nonzero_bits (v - m_low->length ());
	}
};

class lowpart_expr : public expr
{
	expr *m_op;
	int m_off;

public:
	lowpart_expr (expr *src, int off, int len)
		: expr (len), m_op (src), m_off (off)
	{
		assert (len < m_op->length ());
	}
	void require_carry () override
	{
		m_op->require_carry ();
	}
	QString next_piece () override
	{
		int idx = m_values.size ();
		int rev = m_len - idx - 1;
		return m_op->get_piece_high (rev + m_off);
	}
};

void gen_cond_negate (expr *dest, expr *src, QString pred)
{
	int len = dest->length ();
	assert (len == src->length ());
	QString lab = codegen->gen_label ("flip");
	for (size_t i = 0; i < len; i++) {
		QString sr = src->get_piece (i);
		QString dr = dest->get_piece (i);
		codegen->append_move ("u32", dr, sr);
	}
	codegen->append_code (QString ("@!%1\tbra\t%2;\n").arg (pred, lab));
	for (size_t i = 0; i < len; i++) {
		QString dr = dest->get_piece (i);
		codegen->append_code (QString ("\tsub%1.cc.u32\t%2, 0, %2;\n").arg (i == 0 ? "" : "c", dr));
	}
	codegen->append_code (lab + ":\n");
}

class abs_expr : public expr
{
	expr *m_op;
	QString m_pred, m_pred_reg;
public:
	abs_expr (expr *src, QString pred = QString ())
		: expr (src->length ()), m_op (src), m_pred_reg (pred)
	{
		src->require_carry ();
	}
	QString next_piece () override
	{
		assert (m_values.size () == 0);
		m_op->calculate_full ();

		QString p = get_pred ();
		for (size_t i = 0; i < length (); i++) {
			QString dst = get_destreg ("u32");
			m_values.push_back (dst);
		}
		gen_cond_negate (this, m_op, p);
		return m_values[0];
	}
	QString get_pred ()
	{
		if (m_pred.isEmpty ()) {
			m_pred = m_pred_reg.isEmpty () ? codegen->gen_reg ("pred", "isn") : m_pred_reg;
			QString high = m_op->get_piece_high (0);
			codegen->append_code (QString ("\tsetp.lt.s32\t%1, %2, 0;\n").arg (m_pred, high));
		}
		return m_pred;
	}
};

class mult_sign_fixup_expr : public expr
{
	expr *m_op;
	QString m_sgn1, m_sgn2;

public:
	mult_sign_fixup_expr (expr *op, abs_expr *sgn1, abs_expr *sgn2)
		: expr (op->length ()), m_op (op), m_sgn1 (sgn1->get_pred ()), m_sgn2 (sgn2->get_pred ())
	{
	}
	mult_sign_fixup_expr (expr *op, const QString &pred1, const QString &pred2)
		: expr (op->length ()), m_op (op), m_sgn1 (pred1), m_sgn2 (pred2)
	{
	}
	QString next_piece () override
	{
		assert (m_values.size () == 0);
		m_op->calculate_full ();

		QString pred = codegen->gen_reg ("pred", "flip");
		QString lab = codegen->gen_label ("flip");
		codegen->append_code (QString ("\txor.pred\t%1, %2, %3;\n").arg (pred, m_sgn1, m_sgn2));
		for (size_t i = 0; i < length (); i++) {
			QString dst = get_destreg ("u32");
			m_values.push_back (dst);
		}
		gen_cond_negate (this, m_op, pred);
		return m_values[0];
	}
};


class mult_expr : public expr
{
	int m_parts_len;
	expr *m_a, *m_b;
	QString m_carry, m_carry2;
	vector<QString> m_preds_a, m_preds_b;
	bool m_skip_ones;
public:
	mult_expr (expr *srca, expr *srcb, int discard_top = 1, bool skip_ones = false) :
		expr (2 * std::max (srca->length (), srcb->length ()) - discard_top), m_a (srca), m_b (srcb),
		m_skip_ones (skip_ones)
	{
		m_parts_len = std::max (srca->length (), srcb->length ());
		srca->require_carry ();
		srcb->require_carry ();
		m_carry = "0";
		m_carry2 = "0";
	}
	QString next_piece () override
	{
		bool sqr = m_a == m_b;
		int idx = m_values.size ();
		int size = m_parts_len;
		int n_prods = idx < size ? idx + 1 : size * 2 - idx - 1;
		int first_i = idx < size ? 0 : idx - size + 1;

		QString res_in = m_carry;
		QString carry_in = m_carry2;

		QString res;
		if (n_prods == 0) {
			res = get_destreg ("u32");
			codegen->append_move ("u32", res, res_in);
			return res;
		}
		int n_steps = sqr ? (n_prods + 1) / 2 : n_prods;
		bool use_inreg = !sqr || n_prods == 1;
		if (!use_inreg) {
			for (int j = 0; j < n_steps; j++) {
				int idx1 = first_i + j;
				int idx2 = first_i + n_prods - j - 1;
				QString t1 = m_a->get_piece_high (m_parts_len - idx1 - 1);
				QString t2 = m_b->get_piece_high (m_parts_len - idx2 - 1);
				uint32_t nza = m_a->nonzero_bits_high (m_parts_len - idx1 - 1);
				uint32_t nzb = m_b->nonzero_bits_high (m_parts_len - idx2 - 1);
				if (t1 != "0" && t2 != "0" && ((nza != 1 && nzb != 1) || !m_skip_ones))
					break;
				if (j + 1 == n_prods / 2) {
					use_inreg = true;
					break;
				}
			}
		}
		QString current_r = use_inreg ? res_in : "0";
		QString current_c = use_inreg ? carry_in : "0";
		QString current_c2 = "0";
		for (int j = 0; j < n_steps; j++) {
			int idx1 = first_i + j;
			int idx2 = first_i + n_prods - j - 1;
			QString t1 = m_a->get_piece_high (m_parts_len - idx1 - 1);
			QString t2 = m_b->get_piece_high (m_parts_len - idx2 - 1);
			uint32_t nza = m_a->nonzero_bits_high (m_parts_len - idx1 - 1);
			uint32_t nzb = m_b->nonzero_bits_high (m_parts_len - idx2 - 1);
			if (t1 != "0" && t2 != "0" && ((nza != 1 && nzb != 1) || !m_skip_ones)) {
				if (res.isEmpty ()) {
					res = get_destreg ("u32");
					m_carry = codegen->gen_reg ("u32", "carry");
					m_carry2 = codegen->gen_reg ("u32", "carry2_");
				}
				if (current_r == "0") {
					codegen->append_code (QString ("\tmul.lo.u32 %1, %2, %3;\n").arg (res, t1, t2));
					codegen->append_code (QString ("\tmul.hi.u32 %1, %2, %3;\n").arg (m_carry, t1, t2));
					codegen->append_code (QString ("\tmov.u32 %1, 0;\n").arg (m_carry2));
				} else {
					codegen->append_code (QString ("\tmad.lo.cc.u32 %1, %2, %3, %4;\n")
							      .arg (res, t1, t2, current_r));
					codegen->append_code (QString ("\tmadc.hi.cc.u32 %1, %2, %3, %4;\n")
							      .arg (m_carry, t1, t2, current_c));
					codegen->append_code (QString ("\taddc.u32 %1, %2, 0;\n")
							      .arg (m_carry2, current_c2));
				}
				current_r = res;
				current_c = m_carry;
				current_c2 = m_carry2;
			}
			if (!use_inreg && j + 1 == n_prods / 2) {
				codegen->append_code (QString ("\tadd.cc.u32 %1, %1, %1;\n").arg (res));
				codegen->append_code (QString ("\taddc.cc.u32 %1, %1, %1;\n").arg (m_carry));
				codegen->append_code (QString ("\taddc.u32 %1, %1, %1;\n").arg (m_carry2));
			}
		}
		if (current_r == "0") {
			if (!use_inreg) {
				current_r = res_in;
				current_c = carry_in;
				current_c2 = "0";
			}
		} else if (!use_inreg) {
			codegen->append_code (QString ("\tadd.cc.u32 %1, %1, %2;\n").arg (res, res_in));
			codegen->append_code (QString ("\taddc.cc.u32 %1, %1, %2;\n").arg (m_carry, carry_in));
			codegen->append_code (QString ("\taddc.u32 %1, %1, 0;\n").arg (m_carry2));
		}
		/* Rounding */
		if (0 && idx + 2 == size) {
			codegen->append_code (QString ("\tadd.cc.u32 %1, %1, 0x80000000;\n").arg (res));
			codegen->append_code (QString ("\taddc.cc.u32 %1, %1, 0;\n").arg (m_carry));
			codegen->append_code (QString ("\taddc.u32 %1, %1, 0;\n").arg (m_carry2));
		}
		codegen->append_code ("\n");
		m_carry = current_c;
		m_carry2 = current_c2;
		return current_r;
	}
};

/* An alternative way to do a plain non-Karatsuba multiplication.
   This one uses the nonzero_bits mechanism to identify cases where one of the words involved
   is either zero or one: in that case, a conditional add is all that's required, rather than
   a multiplication.
   This case arises in the Karatsuba algorithm, where (al + ah) is multiplied by (bl + bh),
   and the sums are one word wider than their addends, but the high word is always either zero
   or one.  */
class mult_special_expr : public expr
{
	expr *m_a, *m_b, *m_tmpm;

public:
	mult_special_expr (expr *srca, expr *srcb, int discard_top = 1) :
		expr (2 * std::max (srca->length (), srcb->length ()) - discard_top), m_a (srca), m_b (srcb)
	{
		srca->require_carry ();
		srcb->require_carry ();
		m_tmpm = make<mult_expr> (m_a, m_b, discard_top, true);

	}
	QString next_piece () override
	{
		// Generate everything in one go.
		if (m_values.size () > 0)
			abort ();
		m_tmpm->set_destreg (m_dstreg);
		bool any_found = false;
		size_t alen = m_a->length ();
		size_t blen = m_b->length ();
		for (size_t i = 0; i < alen; i++)
			any_found |= m_a->nonzero_bits (i) == 1;
		for (size_t i = 0; i < blen; i++)
			any_found |= m_b->nonzero_bits (i) == 1;
		for (size_t i = 0; i < m_tmpm->length (); i++) {
			QString p = m_tmpm->get_piece (i);
			if (p == "0" && any_found) {
				p = codegen->gen_reg ("u32", "mresz");
				codegen->append_move ("u32", p, "0");
			}
			m_values.push_back (p);
		}
		if (!any_found)
			return m_values[0];
		for (size_t i = 0; i < alen; i++)
			if (m_a->nonzero_bits_high (i) == 1) {
				QString pred = codegen->gen_reg ("pred", "one");
				QString lab = codegen->gen_label ("onelab");
				codegen->append_code (QString ("\tsetp.ne.u32\t%1, %2, 0;\n")
						      .arg (pred, m_a->get_piece_high (i)));
				codegen->append_code (QString ("@!%1\tbra\t%2;\n").arg (pred, lab));
				bool usecc = false;
				for (size_t j = 0; j < blen + 1; j++) {
					size_t idx = (alen - i - 1) + j;
					if (idx >= m_values.size ())
						break;
					QString src = j < blen ? m_b->get_piece (j) : "0";
					if (src != "0" || usecc)
						codegen->append_code (QString ("\tadd%1.cc.u32\t%2, %2, %3;\n")
								      .arg (usecc ? "c" : "", m_values[idx], src));
					usecc = src != "0" || usecc;
				}
				codegen->append_code (lab + ":\n");
			} 
		for (size_t i = 0; i < blen; i++)
			if (m_b->nonzero_bits_high (i) == 1) {
				QString pred = codegen->gen_reg ("pred", "one");
				QString lab = codegen->gen_label ("onelab");
				codegen->append_code (QString ("\tsetp.ne.u32\t%1, %2, 0;\n")
						      .arg (pred, m_b->get_piece_high (i)));
				codegen->append_code (QString ("@!%1\tbra\t%2;\n").arg (pred, lab));
				bool usecc = false;
				for (size_t j = 0; j < alen + 1; j++) {
					size_t idx = (blen - i - 1) + j;
					if (idx >= m_values.size ())
						break;
					QString src = j < alen ? m_a->get_piece (j) : "0";
					// Avoid duplicate pieces.
					if (m_a->nonzero_bits (j) == 1)
						src = "0";
					if (src != "0" || usecc)
						codegen->append_code (QString ("\tadd%1.cc.u32\t%2, %2, %3;\n")
								      .arg (usecc ? "c" : "", m_values[idx], src));
					usecc = src != "0" || usecc;
				}
				codegen->append_code (lab + ":\n");
			}
		return m_values[0];
	}
};

void gen_store (const QString &dst, expr *ex)
{
	int len = ex->length ();
	QString addr = dst; // codegen->gen_reg ("u64", "addr");
	// codegen->append_move ("u64", addr, dst);
	for (int i = 0; i < len; i++) {
		QString val = ex->get_piece (i);
		codegen->append_code (QString ("\tst.global.u32\t[%1 + %2], %3;\n")
				      .arg (addr).arg (i * 4).arg (val));
	}
}

void gen_store (reg_expr *dst, expr *ex)
{
	int len = ex->length ();
	for (int i = 0; i < len; i++) {
		QString dstp = dst->get_piece (i);
		QString val = ex->get_piece (i);
		codegen->append_move ("u32", dstp, val);
	}
}

// Store v into scratch, normalized so that 1/2 <= stored < 1
// We store one extra word so as to not lose precision. scratch needs to be
// sized appropriately.
QString emit_normalize (expr *v, const QString &scratch)
{
	QString predreg = codegen->gen_reg ("pred", "nzfound");
	QString ptr = codegen->gen_reg ("u32", "ptr");
	QString endptr = codegen->gen_reg ("u32", "endptr");
	QString amt1 = codegen->gen_reg ("u32", "amt1");
	QString amt = codegen->gen_reg ("u32", "amt");
	QString lastp = codegen->gen_reg ("u32", "lastp");
	QString store = codegen->gen_reg ("u32", "store");

	int len = v->length ();

	codegen->append_move ("u32", ptr, scratch);
	codegen->append_move ("u32", amt, "0");
	codegen->append_move ("u32", amt1, "0");
	codegen->append_move ("u32", lastp, "0");
	codegen->append_code (QString ("\tadd.u32\t%1, %2, %3;\n").arg (endptr, ptr).arg (len * 4));
	codegen->append_code (QString ("\tsetp.ne.u32\t%1, 0, 0;\n").arg (predreg));

	for (int i = 0; i < len + 1; i++) {
		QString p = i >= len ? "0" : v->get_piece_high (i);
		QString lab = codegen->gen_label ("pastfirst");
		codegen->append_code (QString ("@%1\tbra\t\t%2;\n").arg (predreg, lab));
		codegen->append_code (QString ("\tsetp.ne.u32\t%1, %2, 0;\n").arg (predreg, p));
		codegen->append_code (QString ("@!%1\tadd.u32\t\t%2, %2, 32;\n").arg (predreg, amt1));
		codegen->append_code (QString ("@%1\tbfind.shiftamt.u32\t%2, %3;\n").arg (predreg, amt, p));
		codegen->append_code (lab + ":\n");
		codegen->append_code (QString("\tshf.l.clamp.b32\t%1,%2,%3,%4;\n")
				      .arg (store, p, lastp, amt));
		codegen->append_code (QString("@%1\tst.shared.u32\t[%2], %3;\n")
				      .arg (predreg, endptr, store));
		codegen->append_code (QString("@!%1\tst.shared.u32\t[%2], 0;\n")
				      .arg (predreg, ptr));
		codegen->append_code (QString("@!%1\tadd.u32\t%2, %2, 4;\n")
				      .arg (predreg, ptr));
		codegen->append_code (QString("@%1\tsub.u32\t%2, %2, 4;\n")
				      .arg (predreg, endptr));
		codegen->append_move ("u32", lastp, p);
	}
	codegen->append_code (QString ("\tadd.u32\t\t%1, %1, %2;\n").arg (amt, amt1));
	return amt;
}

expr *emit_normalize_amt (expr *v, const QString &amt, const QString &scratch, int len)
{
	QString amtreg = codegen->gen_reg ("u32", "amtcp");
	codegen->append_move ("u32", amtreg, amt);
	QString predreg = codegen->gen_reg ("pred", "skip");
	QString predreg2 = codegen->gen_reg ("pred", "skip");
	QString ptr = codegen->gen_reg ("u32", "ptr");
	QString endptr = codegen->gen_reg ("u32", "endptr");
	QString lastp = codegen->gen_reg ("u32", "lastp");
	QString store = codegen->gen_reg ("u32", "store");
	QString discarded = codegen->gen_reg ("u32", "discarded");

	int vlen = v->length ();

	codegen->append_move ("u32", ptr, scratch);
	codegen->append_move ("u32", lastp, "0");
	codegen->append_move ("u32", discarded, "0");
	codegen->append_code (QString ("\tadd.u32\t%1, %2, %3;\n").arg (endptr, ptr).arg (len * 4 - 4));

	QString endlab = codegen->gen_label ("stend");
	for (int i = 0; i < vlen; i++) {
		QString p = v->get_piece_high (i);
		codegen->append_code (QString ("\tsetp.lt.u32\t%1, %2, 32;\n").arg (predreg, amtreg));
		codegen->append_code (QString("@!%1\tsub.u32\t%2, %2, 32;\n")
				      .arg (predreg, amtreg));
		codegen->append_code (QString("\tshf.l.wrap.b32\t%1,%2,%3,%4;\n")
				      .arg (store, p, lastp, amtreg));
		codegen->append_code (QString("@%1\tst.shared.u32\t[%2], %3;\n")
				      .arg (predreg, endptr, store));
		codegen->append_code (QString("@!%1\tor.b32\t%2, %2, %3;\n")
				      .arg (predreg, discarded, store));
		if (i + 1 >= len) {
			codegen->append_code (QString ("\tsetp.eq.u32\t%1, %2, %3;\n").arg (predreg2, ptr, endptr));
			codegen->append_code (QString ("@%1\tbra\t%2;\n").arg (predreg2, endlab));
		}

		codegen->append_code (QString("@%1\tsub.u32\t%2, %2, 4;\n")
				      .arg (predreg, endptr));
		codegen->append_move ("u32", lastp, p);
	}
	codegen->append_code (endlab + ":\n");
	auto v1 = make<lds_expr> (scratch, "4", len);
	auto reg = make<reg_expr> (len, "namt");
	gen_store (reg, v1);
	QString allok = codegen->gen_label ("allok");
	codegen->append_code (QString ("\tsetp.eq.u32\t%1, %2, 0;\n").arg (predreg2, discarded));
	codegen->append_code (QString ("@%1\tbra\t%2;\n").arg (predreg2, allok));
	codegen->append_move ("u32", reg->get_piece_high (0), "0x7FFFFFFF");
	codegen->append_code (allok + ":\n");
	return reg;
}

constexpr int karatsuba_cutoff = 12;
/*
    n := half length
    F := xh * yh
    G := xl * yl
    H := (xh + xl) (yh + yl)
    K := H − F − G == xh * yl + xl * yh
    result = F << 2n + K << n + G == (F << 2n) + (H << n) - (F << n) - (G << n) + G
 */
expr *gen_mult_karatsuba_1 (expr *a, expr *b)
{
	if (a->length () < karatsuba_cutoff) {
		return make<mult_special_expr> (a, b, 0);
	}
	int half = (a->length () + 1) / 2;
	int full = 2 * half;

	auto ah = make<trunc_expr> (a, half);
	auto bh = a == b ? ah : make<trunc_expr> (b, half);
	auto al = make<lowpart_expr> (a, half, half);
	auto bl = a == b ? al : make<lowpart_expr> (b, half, half);

	auto H1 = make<addsub_expr> ("add", make<padhigh_expr> (ah, half + 1),
					    make<padhigh_expr> (al, half + 1));
	auto H2 = a == b ? H1 : make<addsub_expr> ("add", make<padhigh_expr> (bh, half + 1),
							  make<padhigh_expr> (bl, half + 1));
	auto H = make<lowpart_expr> (gen_mult_karatsuba_1 (H1, H2), 1, full + 1);
	auto F = gen_mult_karatsuba_1 (ah, bh);
	auto G = gen_mult_karatsuba_1 (al, bl);
	assert (F->length () == full);
	assert (G->length () == full);
	assert (H->length () == full + 1);
	F->set_destreg ("F");
	G->set_destreg ("G");
	H->set_destreg ("H");
	auto sub1 = make<addsub_expr> ("sub", H, make<padhigh_expr> (G, full + 1));
	auto sub2 = make<addsub_expr> ("sub", sub1, make<padhigh_expr> (F, full + 1));
	sub2->set_destreg ("K");
	auto highlow = make<concat_expr> (G, F);
	auto Ke = make<padlow_expr> (make<padhigh_expr> (sub2, full + half), 2 * full);
	auto sum = make<addsub_expr> ("add", highlow, Ke);
	sum->set_destreg ("msum");
	if (half * 2 > a->length ())
		return make<trunc_expr> (sum, 2 * a->length ());
	return sum;
}

expr *gen_mult_unsigned (expr *a, expr *b, bool truncate = true)
{
	expr *raw_mult;
	if (a->length () < karatsuba_cutoff)
		raw_mult = make<mult_expr> (a, b);
	else {
		auto me = gen_mult_karatsuba_1 (a, b);
		raw_mult = make<lowpart_expr> (me, 1, me->length () - 1);
	}
	if (truncate)
		return make<trunc_expr> (raw_mult, a->length ());
	return raw_mult;
}

struct real_val;

/* A structure for information about a real-valued reg, and optionally some of its
   related values: abs_val and neg_pred can hold a sign/magnitude representation,
   which is often useful, and square_val can be used to keep the square of the
   value around.  */
class real_reg
{
	reg_expr *m_reg;
	reg_expr *m_abs_val;
	reg_expr *m_square_val;
	QString m_neg_pred;
	friend class real_val;

public:
	real_reg () = default;
	real_reg (int size, const QString &name, bool keep_abs = false, bool keep_square = false)
		: m_reg (make<reg_expr> (size, name)),
		  m_abs_val (keep_abs || keep_square ? make<reg_expr> (size, name + "a") : nullptr),
		  m_square_val (keep_square ? make<reg_expr> (size, name + "sq") : nullptr),
		  m_neg_pred (keep_abs ? codegen->gen_reg ("pred", name + "n") : nullptr)
	{
	}
	void store (expr *val)
	{
		gen_store (m_reg, val);
		if (!m_neg_pred.isEmpty () || m_square_val != nullptr) {
			auto absexp = make<abs_expr> (val, m_neg_pred);
			if (!m_neg_pred.isEmpty ())
				gen_store (m_abs_val, absexp);
			if (m_square_val != nullptr)
				gen_store (m_square_val, gen_mult_unsigned (absexp, absexp));
		}
	}
	operator expr *() const
	{
		return m_reg;
	}
	expr *ex () const
	{
		return m_reg;
	}
	bool have_squared () const
	{
		return m_square_val != nullptr;
	}
	expr *squared ()
	{
		if (m_square_val == nullptr) {
			if (m_abs_val != nullptr)
				return gen_mult_unsigned (m_abs_val, m_abs_val);
			auto absv = abs_val ();
			return gen_mult_unsigned (absv, absv);
		}
		return m_square_val;
	}
	expr *abs_val () const
	{
		if (m_abs_val == nullptr)
			return make<abs_expr> (m_reg, m_neg_pred);
		return m_abs_val;
	}
	QString neg_pred () const
	{
		if (!m_neg_pred.isEmpty ())
			return m_neg_pred;

		QString pred = codegen->gen_reg ("pred", "neg");
		codegen->append_code (QString ("\tsetp.lt.s32\t%1, %2, 0;\n").arg (pred, m_reg->get_piece_high (0)));
		return pred;
	}
};

class real_val
{
	expr *m_val {};
	expr *m_abs_val {};
	expr *m_square_val {};
	QString m_neg_pred;

public:
	real_val (const real_reg &reg) : m_val (reg.m_reg), m_abs_val (reg.m_abs_val), m_square_val (reg.m_square_val), m_neg_pred (reg.m_neg_pred)
	{
	}
	real_val (expr *expr) : m_val (expr)
	{
	}
	real_val (expr *expr, const QString &pred) : m_abs_val (expr), m_neg_pred (pred)
	{
	}
	real_val () = default;
	operator expr *() const
	{
		return ex ();
	}
	expr *ex () const
	{
		if (m_val != nullptr)
			return m_val;
		expr *dst = make<reg_expr> (m_abs_val->length (), "sgn");
		gen_cond_negate (dst, m_abs_val, m_neg_pred);
		return dst;
	}
	bool have_squared () const
	{
		return m_square_val != nullptr;
	}
	expr *squared ()
	{
		if (m_square_val == nullptr) {
			auto absv = abs_val ();
			m_square_val = gen_mult_unsigned (absv, absv);
		}
		return m_square_val;
	}
	expr *abs_val ()
	{
		if (m_abs_val == nullptr)
			m_abs_val = make<abs_expr> (m_val, m_neg_pred);
		return m_abs_val;
	}
	QString neg_pred ()
	{
		if (m_neg_pred.isEmpty ()) {
			m_neg_pred = codegen->gen_reg ("pred", "neg");
			codegen->append_code (QString ("\tsetp.lt.s32\t%1, %2, 0;\n")
					      .arg (m_neg_pred, m_val->get_piece_high (0)));
		}
		return m_neg_pred;
	}
};

template<class R1, class R2>
real_val gen_mult (R1 &a, R2 &b, bool abs_required = true, bool truncate = true)
{
	expr *aop = a.abs_val ();
	expr *bop = b.abs_val ();
	auto prod = gen_mult_unsigned (aop, bop, truncate);
	if (!abs_required || &a == &b)
		return prod;

	QString pred = codegen->gen_reg ("pred", "sfix");
	codegen->append_code (QString ("\txor.pred\t%1, %2, %3;\n").arg (pred, a.neg_pred (), b.neg_pred ()));
	return real_val (prod, pred);
}

QString convert_to_float (expr *v)
{
	QString dstreg = codegen->gen_reg ("f64", "convdst");
	int len = v->length ();
	QString conv = codegen->gen_reg ("f64", "fltconv");
	QString mult = codegen->gen_reg ("f64", "multiplier");
	// 2^-32.
	codegen->append_move ("f64", conv, "2.3283064365386963e-10");

	codegen->append_move ("f64", mult, "1.0");
	bool first = true;
	while (len-- > 0) {
		QString ftmp = codegen->gen_reg ("f64", "tmp");
		codegen->append_code (QString ("\tcvt.rn.f64.%1\t%2, %3;\n")
				      .arg (first ? "s32" : "u32", ftmp, v->get_piece (len)));
		if (!first)
			codegen->append_code (QString ("\tmul.f64\t%1, %1, %2;\n").arg (ftmp, mult));
		if (first)
			codegen->append_code (QString ("\tmov.f64\t%1, %2;\n").arg (dstreg, ftmp));
		else
			codegen->append_code (QString ("\tadd.f64\t%1, %1, %2;\n").arg (dstreg, ftmp));
		if (len > 0)
			codegen->append_code (QString ("\tmul.f64\t%1, %1, %2;\n").arg (mult, conv));
		first = false;
	}
	return dstreg;
}

real_val emit_newton_iteration (real_val &d, real_val &xi, int dlen, const QString &name)
{
	// Compute Xi + Xi (1 - D * Xi)
	auto mult = gen_mult_unsigned (d, xi);
	auto const1 = make<const_expr<1>> (mult->length ());
	real_val sub (make<addsub_expr> ("sub", const1, mult));
	real_val subt (make<trunc_expr> (sub, dlen));
	auto result = make<addsub_expr> ("add", xi, gen_mult (subt, xi));
	result->set_destreg (name);
	return real_val { result };
}

expr *invert (const QString &addr, int len)
{
	real_val d (make<lds_expr> (addr, "4", len));
	int est_n = len > 2 ? 3 : 2;
	auto d2 = make<trunc_expr> (d, est_n);
	QString flt = convert_to_float (d2);
	codegen->append_code (QString ("\trcp.rz.f64\t%1, %1;\n").arg (flt));
	real_val estimate (make<padlow_expr> (make<from_float_expr> (flt), len));
	for (int nwords = 1; nwords + 1 < len; nwords *= 2) {
		estimate = emit_newton_iteration (d, estimate, len, QString ("est%1").arg (nwords));
	}
	if (estimate.ex ()->length () > len)
		return make<trunc_expr> (estimate, len);
	return estimate;
}

std::pair<expr *, QString> gen_inverse (expr *d, const QString &scratch)
{
	QString shift_amt = emit_normalize (d, scratch);
	QString zpred = codegen->gen_reg ("pred", "divz");
	int len = d->length ();
	// QString lab = codegen->gen_label ("divz");
	// codegen->append_code (QString ("\tsetp.ge.u32\t%1, %2, %3;\n").arg (zpred, shift_amt).arg (len * 32));
	// codegen->append_code (QString ("@%1\tbra\t\t%2;\n").arg (zpred, lab));
	auto inverse = invert (scratch, len + 1);
	inverse->calculate_full ();
	// codegen->append_code (QString ("%1:\n").arg (lab));
	return { inverse, shift_amt };
}

expr *gen_div (expr *n, expr *d, const QString &scratch)
{
	auto na = make<abs_expr> (n);
	auto da = make<abs_expr> (d);
	auto [ inverse, shift_amt ] = gen_inverse (da, scratch);
	int len = da->length ();
	auto nanorm = emit_normalize_amt (na, shift_amt, scratch, len + 1);
	nanorm->set_destreg ("nanorm");
	auto q1 = make<mult_expr> (nanorm, inverse);
	q1->set_destreg ("q1_");
	auto q1r = make<addsub_expr> ("add", q1, make<round_cst_expr> (1 + len, 0x80000000));
	q1r->set_destreg ("q1r_");
#if 0
	/* Now we have an estimate of the result of the division. Perform one more iteration:
	   q2 = q1 + inv(n - d*q1).  */
	auto dapad = make<padlow_expr> (da, len + 1);
	auto dq1 = gen_mult (dapad, make<trunc_expr> (q1, len + 1));
	auto nmdq1 = make<addsub_expr> ("sub", nanorm, dq1);
	auto addend = make<mult_expr> (nmdq1, inverse);
	auto q2 = make<addsub_expr> ("add", q1, addend);
	auto q2r = make<addsub_expr> ("add", q2, make<round_cst_expr> (1 + len, 0x80000000));
#endif
	// auto pos_result = emit_normalize_amt (q1, shift_amt, scratch, len + 1);
	return make<mult_sign_fixup_expr> (make<trunc_expr> (q1r, len), na, da);
}

void gen_coord_muladd (QString &result, int dstsize, int srcsize)
{
	result += ".func coord_muladd (.reg.u64 %dst, .reg.u64 %src_step, .reg.u64 %src_base, .reg.u32 %srcc)\n{\n";
	result += "\t.reg.pred\t%neg;\n";
	result += "\tsetp.lt.s32\t%neg, %srcc, 0;\n";
	result += "@%neg\tneg.s32\t%srcc, %srcc;\n";

	QString last = "0";
	QString last_carry = "0";
	vector<QString> pieces;
	for (int i = 0; i < srcsize; i++) {
		QString t1 = QString ("%vala%1").arg (i);
		QString res = QString ("%res%1").arg (i);
		QString carry = QString ("%carry%1").arg (i);
		result += QString ("\t.reg.u32 %1, %2, %3;\n").arg (t1, res, carry);
		result += QString ("\tld.global.u32 %1, [%src_step + %2];\n").arg (t1).arg (i * 4);
		result += QString ("\tmad.lo.cc.u32 %1, %2, %srcc, %3;\n").arg (res, t1, last);
		result += QString ("\tmadc.hi.u32 %1, %2, %srcc, 0;\n").arg (carry, t1);
		last = carry;
		pieces.push_back (res);
		result += "\n";
	}
	for (int i = 0; i < srcsize; i++) {
		QString t2 = QString ("%valb%1").arg (i);
		result += QString ("\t.reg.u32\t%1;\n").arg (t2);
		result += QString ("\tld.global.u32 %1, [%src_base + %2];\n").arg (t2).arg (i * 4);
		auto res = pieces[i];
		result += QString ("@%neg\tsub%3.cc.u32\t%1, %1, %2;\n@!%neg\tadd%3.cc.u32\t%1, %1, %2;\n").arg (t2, res, i == 0 ? "" : "c");
		if (i >= srcsize - dstsize)
			result += QString ("\tst.global.u32 [%dst + %1], %2;\n").arg ((i - srcsize + dstsize) * 4).arg (t2);
	}
	result += "}\n";
}

void gen_coord_mul (QString &result, int dstsize, int srcsize)
{
	result += ".func coord_mul (.reg.u64 %dst, .reg.u64 %src_step, .reg.u32 %srcc)\n{\n";
	result += "\t.reg.pred\t%neg;\n";
	result += "\tsetp.lt.s32\t%neg, %srcc, 0;\n";
	result += "@%neg\tneg.s32\t%srcc, %srcc;\n";

	QString last = "0";
	QString last_carry = "0";
	vector<QString> pieces;
	for (int i = 0; i < srcsize; i++) {
		QString t1 = QString ("%vala%1").arg (i);
		QString res = QString ("%res%1").arg (i);
		QString carry = QString ("%carry%1").arg (i);
		result += QString ("\t.reg.u32 %1, %2, %3;\n").arg (t1, res, carry);
		result += QString ("\tld.global.u32 %1, [%src_step + %2];\n").arg (t1).arg (i * 4);
		result += QString ("\tmad.lo.cc.u32 %1, %2, %srcc, %3;\n").arg (res, t1, last);
		result += QString ("\tmadc.hi.u32 %1, %2, %srcc, 0;\n").arg (carry, t1);
		last = carry;
		pieces.push_back (res);
		result += "\n";
	}
	for (int i = 0; i < srcsize; i++) {
		auto res = pieces[i];
		result += QString ("@%neg\tsub%2.cc.u32\t%1, 0, %1;\n").arg (res, i == 0 ? "" : "c");
		if (i >= srcsize - dstsize)
			result += QString ("\tst.global.u32 [%dst + %1], %2;\n").arg ((i - srcsize + dstsize) * 4).arg (res);
	}
	result += "}\n";
}

void gen_mul_func (QString &result, int size, bool sqr)
{
#if 0
	/* ??? This would be nice, but the multiple return values cause the following:
	   ptxas fatal   : 'Compile only(-c)' requires ABI
	   It's unclear how to avoid this.  */
	generator cg;
	codegen = &cg;
	auto dst = make<reg_expr> (size, "dst", false);
	auto srca = make<reg_expr> (size, "srca", false);
	auto srcb = sqr ? srca : make<reg_expr> (size, "srcb", false);
	result += ".func (";
	for (int i = 0; i < dst->length (); i++) {
		if (i > 0)
			result += ", ";
		result += QString (".reg.u32 %1").arg (dst->get_piece (i));
	}
	if (sqr)
		result += ") sqr_real (";
	else
		result += ") mul_real (";
	for (int i = 0; i < srca->length (); i++) {
		if (i > 0)
			result += ", ";
		result += QString (".reg.u32 %1").arg (srca->get_piece (i));
	}
	if (!sqr)
		for (int i = 0; i < srcb->length (); i++) {
			result += ", ";
			result += QString (".reg.u32 %1").arg (srcb->get_piece (i));
		}
	result += ")\n{\n";

	gen_store (dst, gen_mult (srca, srcb));
	result += cg.code ();
	codegen = nullptr;
	result += "}\n";
#else
	result += "\t.reg.u32 %res_in, %carry_in, %carry, %carry2;\n";
	result += "\tmov.u32 %carry, 0;\n";
	result += "\tmov.u32 %carry2, 0;\n";

	for (int i = 0; i < size; i++) {
		QString t1 = QString ("%vala%1").arg (i);
		QString t2 = QString ("%valb%1").arg (i);
		result += QString ("\t.reg.u32 %1;\n").arg (sqr ? t1 : t1 + ", " + t2);
		result += QString ("\tld.global.u32 %1, [%srca + %2];\n").arg (t1).arg (i * 4);
		if (!sqr)
			result += QString ("\tld.global.u32 %1, [%srcb + %2];\n").arg (t2).arg (i * 4);
	}
	result += "\t.reg.pred %sgna;\n";
	result += QString ("\tsetp.lt.s32 %sgna, %vala%1, 0;\n").arg (size - 1);
	for (int i = 0; i < size; i++) {
		QString t1 = QString ("%vala%1").arg (i);
		result += QString ("\t@%sgna sub%1.cc.u32 %2, 0, %3;\n").arg (i == 0 ? "" : "c").arg (t1).arg (t1);
	}
	if (!sqr) {
		result += "\t.reg.pred %sgnb;\n";
		result += QString ("\tsetp.lt.s32 %sgnb, %valb%1, 0;\n").arg (size - 1);
		for (int i = 0; i < size; i++) {
			QString t1 = QString ("%valb%1").arg (i);
			result += QString ("\t@%sgnb sub%1.cc.u32 %2, 0, %3;\n").arg (i == 0 ? "" : "c").arg (t1).arg (t1);
		}
	}

	int count = 0;
	for (int i = 0; i < size * 2 - 1; i++) {
		int n_prods = i < size ? i + 1 : size * 2 - i - 1;
		int first_i = i < size ? 0 : i - size + 1;
		QString res = QString ("%res%1").arg (i);
		result += QString ("\t.reg.u32 %1;\n").arg (res);
		result += QString ("\tmov.u32 %res_in, %carry;\n");
		result += "\tmov.u32 %carry_in, %carry2;\n";
		QString initial_r = sqr ? "0" : "%res_in";
		QString initial_c = sqr ? "0" : "%carry_in";
		int n_steps = sqr ? (n_prods + 1) / 2 : n_prods;
		for (int j = 0; j < n_steps; j++) {
			int idx1 = first_i + j;
			int idx2 = first_i + n_prods - j - 1;
			QString t1 = QString ("%vala%1").arg (idx1);
			QString t2 = QString (sqr ? "%vala%1" : "%valb%1").arg (idx2);
			result += QString ("\tmad.lo.cc.u32 %1, %2, %3, %4;\n").arg (res).arg (t1).arg (t2).arg (j == 0 ? initial_r : res);
			result += QString ("\tmadc.hi.cc.u32 %carry, %1, %2, %3;\n").arg (t1).arg (t2).arg (j == 0 ? initial_c : "%carry");
			result += QString ("\taddc.u32 %carry2, %1, 0;\n").arg (j == 0 ? "0" : "%carry2");
			count += 2;
			if (sqr && j + 1 == n_prods / 2) {
				result += QString ("\tadd.cc.u32 %1, %2, %3;\n").arg (res).arg (res).arg (res);
				result += "\taddc.cc.u32 %carry, %carry, %carry;\n";
				result += "\taddc.u32 %carry2, %carry2, %carry2;\n";
			}

		}
		if (sqr) {
			result += QString ("\tadd.cc.u32 %1, %2, %res_in;\n").arg (res, res);
			result += "\taddc.cc.u32 %carry, %carry, %carry_in;\n";
			result += "\taddc.u32 %carry2, %carry2, 0;\n";
		}
		/* Rounding */
		if (i + 2 == size) {
			result += QString ("\tadd.cc.u32 %1, %2, 0x80000000;\n").arg (res, res);
			result += "\taddc.cc.u32 %carry, %carry, 0;\n";
			result += "\taddc.u32 %carry2, %carry2, 0;\n";
		}
		if (i + 1 >= size && sqr)
			result += QString ("\tst.global.u32 [%dst + %1], %2;\n").arg ((i - size + 1) * 4).arg (res);
		result += "\n";
	}
	if (!sqr) {
		result += "\t.reg.pred %sgnr;\n";
		result += "\txor.pred %sgnr, %sgna,%sgnb;\n";
		for (int i = 0; i < size; i++) {
			QString t1 = QString ("%res%1").arg (i + size - 1);
			result += QString ("\t@%sgnr sub%1.cc.u32 %2, 0, %3;\n").arg (QString (i == 0 ? "" : "c"), t1, t1);
			result += QString ("\tst.global.u32 [%dst + %1], %2;\n").arg (i * 4).arg (t1);
		}
	}
	result += "}\n";
#endif
}

template<typename... ARGS>
void gen_kernel_header (QString &result, const QString &name, ARGS... in_args)
{
	vector<QString> args = { in_args... };
	result += QString (".visible .entry %1 (").arg (name);
	bool first = true;
	for (size_t i = 0; i < args.size (); i += 2) {
		if (!first)
			result += ",\n\t\t\t";
		result += QString (".param.%1 %in_%2").arg (args[i], args[i + 1]);
		first = false;
	}
	result += ")\n{\n";
	for (size_t i = 0; i < args.size (); i += 2)
		result += QString ("\t.reg.%1 %%2;\n").arg (args[i], args[i + 1]);

	for (size_t i = 0; i < args.size (); i += 2)
		result += QString ("\tld.param.%1 %%2, [%in_%3];\n").arg (args[i]).arg (args[i + 1]).arg (args[i + 1]);
}

struct cplx_val;
struct cplx_reg
{
	real_reg re;
	real_reg im;
	cplx_reg (const real_reg &r, const real_reg &i) : re (r), im (i)
	{
	}
	cplx_reg () = default;
	void store (const cplx_val &v);
};

struct cplx_val
{
	real_val re;
	real_val im;
	cplx_val (const real_val &r, const real_val &i) : re (r), im (i)
	{
	}
	cplx_val (expr *r, expr *i) : re (r), im (i)
	{
	}
	cplx_val (const cplx_reg &reg) : re (reg.re), im (reg.im)
	{
	}
	void set_re (expr *v)
	{
		re = real_val (v);
	}
	void set_im (expr *v)
	{
		im = real_val (v);
	}
	cplx_val () = default;
};

void cplx_reg::store (const cplx_val &v)
{
	re.store (v.re);
	im.store (v.im);
}

constexpr int mult_alg_cutoff = 3;
cplx_val emit_complex_mult (cplx_val &va, cplx_val &vb)
{
	if (va.re.ex ()->length () < mult_alg_cutoff) {
		// (a+ib) * (c+id) = ac - bd + i(bc + ad)
		real_val rp1_expr = gen_mult (va.re, vb.re);
		real_val rp2_expr = gen_mult (va.im, vb.im);

		real_val ip1_expr = gen_mult (va.re, vb.im);
		real_val ip2_expr = gen_mult (va.im, vb.re);

		auto re = make<addsub_expr> ("sub", rp1_expr, rp2_expr);
		auto im = make<addsub_expr> ("add", ip1_expr, ip2_expr);
		return { re, im };
	}
	// k1 = c(a + b)
	// k2 = a(d − c)
	// k3 = b(c + d)
	// re = k1 − k3
	// im = k1 + k2.
	real_val sum1 (make<addsub_expr> ("add", va.re, va.im));
	real_val sum2 (make<addsub_expr> ("sub", vb.im, vb.re));
	real_val sum3 (make<addsub_expr> ("add", vb.re, vb.im));
	real_val k1 = gen_mult (vb.re, sum1);
	real_val k2 = gen_mult (va.re, sum2);
	real_val k3 = gen_mult (va.im, sum3);
	auto re = make<addsub_expr> ("sub", k1, k3);
	auto im = make<addsub_expr> ("add", k1, k2);
	return { re, im };
}

template<class C>
cplx_val emit_complex_sqr (C &v)
{
	if (!v.re.have_squared () && !v.im.have_squared ()) {
		cplx_val vc (v);
		return emit_complex_mult (vc, vc);
	}

	auto sqre = v.re.squared ();
	auto sqim = v.im.squared ();

	// zr <- z2r - z2i
	// zi <- zr * zi * 2

	auto tmp = gen_mult (v.re, v.im);

	auto re = make<addsub_expr> ("sub", sqre, sqim);
	auto im = make<addsub_expr> ("add", tmp, tmp);
	return { re, im };
}

template<class C>
cplx_val gen_mult_int (C &v, int factor)
{
	cplx_val result;
	bool first = true;
	for (int i = 0; i < 10; i++)
		if (factor & (1 << i)) {
			expr *shfr = i == 0 ? v.re.ex () : make<lshift_expr> (v.re, i);
			expr *shfi = i == 0 ? v.im.ex () : make<lshift_expr> (v.im, i);
			if (first) {
				result = { shfr, shfi };
				first = false;
			} else {
				auto re = make<addsub_expr> ("add", result.re, shfr);
				auto im = make<addsub_expr> ("add", result.im, shfi);
				result = { re, im };
			}
		}
	return result;
}

expr *gen_mult_by_inverse (expr *v, expr *inverse,
				      const QString &scratch, const QString &shift_amt)
{
	auto va = make<abs_expr> (v);
	int len = inverse->length () - 1;
	auto vanorm = emit_normalize_amt (va, shift_amt, scratch, len + 1);
	vanorm->set_destreg ("vanorm");
	auto q1 = make<mult_expr> (vanorm, inverse);
	q1->set_destreg ("q1_");
	auto q1r = make<addsub_expr> ("add", q1, make<round_cst_expr> (1 + len, 0x80000000));
	q1r->set_destreg ("q1r_");
	// auto pos_result = emit_normalize_amt (q1, shift_amt, scratch, len + 1);
	auto result = make<reg_expr> (len);
	gen_cond_negate (result, make<trunc_expr> (q1r, len), va->get_pred ());
	return result;
}

cplx_val emit_complex_div (cplx_val &n, cplx_val &d)
{
	auto rp1_expr = gen_mult (n.re, d.re, true, false);
	auto rp2_expr = gen_mult (n.im, d.im, true, false);

	auto ip2_expr = gen_mult (n.re, d.im, true, false);
	auto ip1_expr = gen_mult (n.im, d.re, true, false);

	auto nre = make<addsub_expr> ("add", rp1_expr, rp2_expr);
	auto nim = make<addsub_expr> ("sub", ip1_expr, ip2_expr);

	auto d2re = d.re.squared ();
	auto d2im = d.im.squared ();
	auto d_sum = make<addsub_expr> ("add", d2re, d2im);
	auto [ inv, shift_amt ] = gen_inverse (d_sum, "%scratch");
	auto rre = gen_mult_by_inverse (nre, inv, "%scratch", shift_amt);
	rre->calculate_full ();
	auto rim = gen_mult_by_inverse (nim, inv, "%scratch", shift_amt);
	rim->calculate_full ();
	return { rre, rim };
}

template<class C>
void build_powers (array<cplx_val, 20> &powers, C &z, int power)
{
	powers[0] = z;
	powers[1] = emit_complex_sqr (z);
	int pidx = 2;
	for (int p = 4; p <= power; p *= 2, pidx++) {
		powers[pidx] = emit_complex_sqr (powers[pidx - 1]);
	}
}

cplx_val get_power (array<cplx_val, 20> &powers, int power)
{
	cplx_val v;
	bool first = true;
	for (int i = 0; i < 10; i++)
		if (power & (1 << i)) {
			if (first) {
				v = powers[i];
				first = false;
			} else {
				v = emit_complex_mult (powers[i], v);
			}
		}
	return v;
}

cplx_val cplx_ldc (const QString &name, int size, int stepsize)
{
	auto re = make<ldc_expr> (QString ("%1 + %2")
					 .arg (name).arg (stepsize * 4 - size * 4), size);
	auto im = make<ldc_expr> (QString ("%1 + %2")
					 .arg (name).arg (2 * stepsize * 4 - size * 4), size);
	return { re, im };
}

cplx_val emit_complex_add (const cplx_val &v1, const cplx_val &v2)
{
	auto newzr = make<addsub_expr> ("add", v1.re, v2.re);
	auto newzi = make<addsub_expr> ("add", v1.im, v2.im);
	return { newzr, newzi };
}

cplx_val emit_complex_sub (const cplx_val &v1, const cplx_val &v2)
{
	auto newzr = make<addsub_expr> ("sub", v1.re, v2.re);
	auto newzi = make<addsub_expr> ("sub", v1.im, v2.im);
	return { newzr, newzi };
}

cplx_val mul_by_power (cplx_val &v, array<cplx_val, 20> &powers, int p)
{
	if (p == 0)
		return v;
	cplx_val pwr = get_power (powers, p);
	return emit_complex_mult (pwr, v);
}

void gen_inner_standard (int size, int power, bool julia, bool dem,
			 cplx_reg zreg, cplx_val cval, cplx_reg zder)
{
	array<cplx_val, 20> powers;
	build_powers (powers, zreg, power);
	cplx_val pwr = get_power (powers, power);

	cplx_val newz = emit_complex_add (pwr, cval);
	if (dem) {
		cplx_val zdval (zder);
		cplx_val nrd = mul_by_power (zdval, powers, power - 1);
		nrd = gen_mult_int (nrd, power);
		if (!julia) {
			nrd.set_re (make<addsub_expr> ("add", nrd.re,
							      make<ldg_expr> ("%dem_step", size)));
		}
		zder.store (nrd);
	}
	zreg.store (newz);
}

// Burning ship, because why not set zre and zim to their absolute values before
// continuing.
void gen_inner_ship (int power, bool /* julia */, cplx_reg zreg, cplx_val cval)
{
	auto zar = make<abs_expr> (zreg.re);
	auto zai = make<abs_expr> (zreg.im);
	zreg.store ({ zar, zai});
	array<cplx_val, 20> powers;
	build_powers (powers, zreg, power);
	cplx_val pwr = get_power (powers, power);
	cplx_val newz = emit_complex_add (pwr, cval);

	zreg.store (newz);
}

/* Tricorn, also known as Mandelbar.
   Like the default function except the conjugate is taken after the final step
   before assigning Z.  */
void gen_inner_tricorn (int size, int power, bool /* julia */, cplx_reg zreg, cplx_val cval)
{
	array<cplx_val, 20> powers;
	build_powers (powers, zreg, power);
	cplx_val pwr = get_power (powers, power);
	cplx_val newz = emit_complex_add (pwr, cval);

	auto const0 = make<const_expr<0>> (size);
	newz.set_im (make<addsub_expr> ("sub", const0, newz.im));

	zreg.store (newz);
}

void gen_inner_celtic (int power, bool /* julia */, cplx_reg zreg, cplx_val cval)
{
	array<cplx_val, 20> powers;
	build_powers (powers, zreg, power);
	cplx_val pwr = get_power (powers, power);
	pwr.set_re (make<abs_expr> (pwr.re));
	cplx_val newz = emit_complex_add (pwr, cval);

	zreg.store (newz);
}

/* c*(z^N - Nz)
   The critical point is 1.
   The "actual" fractal named Lambda in other programs is c*(z^2 - z), with critical point 1/2.
   Seems to look identical except for scaling.  */

void gen_inner_lambda (int size, int power, bool julia, bool dem,
		       cplx_reg zreg, cplx_val cval, cplx_reg zder)
{
	array<cplx_val, 20> powers;
	build_powers (powers, zreg, power);
	cplx_val pwr = get_power (powers, power);
	cplx_val f = gen_mult_int (zreg, power);

	cplx_val sum = emit_complex_sub (pwr, f);
	cplx_val newz = emit_complex_mult (sum, cval);
	if (dem) {
		cplx_val dpwr = get_power (powers, power - 1);
		cplx_val zdval (zder);
		cplx_val m1 = emit_complex_mult (dpwr, zdval);
		cplx_val nrd = gen_mult_int (m1, power);
		cplx_val f = gen_mult_int (zder, power);
		cplx_val sum2 = emit_complex_sub (nrd, f);
		cplx_val m = emit_complex_mult (sum2, cval);
		if (!julia) {
			real_val step (make<ldg_expr> ("%dem_step", size));
			expr *scaled_re = gen_mult (sum.re, step);
			expr *scaled_im = gen_mult (sum.im, step);
			cplx_val scaled { scaled_re, scaled_im };
			zder.store (emit_complex_add (scaled, m));
		} else {
			zder.store (m);
		}
	}

	zreg.store (newz);
}


void gen_inner_magnet_a (int size, int power, bool /* julia */, cplx_reg zreg, cplx_val creg)
{
	cplx_val pwr = emit_complex_sqr (zreg);
	cplx_val newz1 = emit_complex_add (pwr, creg);
	auto const1 = make<const_expr<1>> (size);
	auto const2 = make<const_expr<2>> (size);
	newz1.set_re (make<addsub_expr> ("sub", newz1.re, const1));
	cplx_val dnewz1 = gen_mult_int (zreg, 2);
	cplx_val dnewz2 = emit_complex_add (dnewz1, creg);
	dnewz2.set_re (make<addsub_expr> ("sub", dnewz2.re, const2));
	auto d = emit_complex_div (newz1, dnewz2);
	auto newz = emit_complex_sqr (d);
	zreg.store (newz);
}

// (z + 1/z)^n
// Generates Mandelbrot shapes facing each other in a ring.
void gen_inner_facing (int size, int power, bool /* julia */, cplx_reg zreg, cplx_val creg)
{
	cplx_val zval = zreg;
	auto q = emit_complex_div (creg, zval);
	auto sum = emit_complex_add (zreg, q);
	array<cplx_val, 20> powers;
	build_powers (powers, sum, power);
	cplx_val pwr = get_power (powers, power);
	zreg.store (pwr);
}

// c(z^n + z^-n), critical point 1. Called "Cczcpaczcp" in Saturn.
// Similar ring shape to the other "Facing" formula.
void gen_inner_facing_b (int size, int power, bool /* julia */, cplx_reg zreg, cplx_val creg)
{
	array<cplx_val, 20> powers;
	build_powers (powers, zreg, power);
	cplx_val pwr = get_power (powers, power);

	auto q = emit_complex_div (creg, pwr);
	auto prod = emit_complex_mult (pwr, creg);
	auto sum = emit_complex_add (prod, q);

	zreg.store (sum);
}

// Another variant of "Cczcpaczcp". This one is c(z^-n - (n/2)z^-2), where n is one more than the chosen power
void gen_inner_rings (int size, int power, bool /* julia */, cplx_reg zreg, cplx_val creg)
{
	/* Compute the division first and then the powers - that's better for maintaining precision.  */
	auto const0 = make<const_expr<0>> (size);
	auto const1 = make<const_expr<1>> (size);
	cplx_val one { const1, const0 };
	cplx_val zval { zreg };
	auto zinv = emit_complex_div (one, zval);

	array<cplx_val, 20> powers;
	build_powers (powers, zinv, power + 1);
	cplx_val pwr = get_power (powers, power + 1);
	cplx_val pwr2 = get_power (powers, 2);

	if (power & 1) {
		cplx_val p2mp = gen_mult_int (pwr2, (power + 1) / 2);
		auto sum = emit_complex_sub (pwr, p2mp);
		auto newz = emit_complex_mult (sum, creg);
		zreg.store (newz);
	} else {
		cplx_val p2m2 = gen_mult_int (pwr2, power + 1);
		auto p2mp2r = make<arshift_expr<1>> (p2m2.re);
		auto p2mp2i = make<arshift_expr<1>> (p2m2.im);
		cplx_val p2mp2 { p2mp2r, p2mp2i };
		auto sum = emit_complex_sub (pwr, p2mp2);
		auto newz = emit_complex_mult (sum, creg);
		zreg.store (newz);
	}
}

void gen_inner_tails (int size, int stepsize, int power, bool julia, cplx_reg zreg, cplx_val cval)
{
	auto const1 = make<const_expr<1>> (size);
	auto const0 = make<const_expr<0>> (size);

	array<cplx_val, 20> powers;
	build_powers (powers, zreg, power);
	cplx_val pwr = get_power (powers, power);

	cplx_val newz1 = emit_complex_add (pwr, cval);
	cplx_val one { const1, const0 };
	cplx_val newz1b = emit_complex_div (one, newz1);
	cplx_val sum = emit_complex_add (newz1, newz1b);
	cplx_val scaled { make<arshift_expr<1>> (sum.re), make<arshift_expr<1>> (sum.im) };
	auto pq = cplx_ldc ("const_param_q", size, stepsize);
	cplx_val newz = emit_complex_sub (scaled, pq);
	zreg.store (newz);
}

void gen_inner_spider (int size, int power, cplx_reg zreg, cplx_reg creg)
{
	array<cplx_val, 20> powers;
	build_powers (powers, zreg, power);
	cplx_val pwr = get_power (powers, power);

	cplx_val newz = emit_complex_add (pwr, creg);

	auto prhalf = make<arshift_expr<1>> (creg.re);
	auto pihalf = make<arshift_expr<1>> (creg.im);
	cplx_val np = emit_complex_add ({ prhalf, pihalf}, newz);
	zreg.store (newz);
	creg.store (np);
}

// Similar to lambda, but computes c(z^N - Nz - 2) + q

void gen_inner_mix (int size, int stepsize, int power, bool julia, bool dem,
		    cplx_reg zreg, cplx_val cval, cplx_reg zder)
{
	array<cplx_val, 20> powers;
	build_powers (powers, zreg, power);
	cplx_val pwr = get_power (powers, power);
	cplx_val f = gen_mult_int (zreg, power);
	cplx_val sum1 = emit_complex_sub (pwr, f);
	cplx_val sum = sum1;
	auto const2 = make<const_expr<2>> (size);
	sum.set_re (make<addsub_expr> ("sub", sum.re, const2));
	cplx_val newz1 = emit_complex_mult (sum, cval);
	auto pq = cplx_ldc ("const_param_q", size, stepsize);
	cplx_val newz = emit_complex_sub (newz1, pq);
	if (dem) {
		cplx_val dpwr = get_power (powers, power - 1);
		cplx_val zdval (zder);
		cplx_val m1 = emit_complex_mult (dpwr, zdval);
		cplx_val nrd = gen_mult_int (m1, power);
		if (!julia) {
			cplx_val sum2 = emit_complex_sub (nrd, gen_mult_int (zder, power));
			cplx_val m = emit_complex_mult (sum2, cval);
			real_val step (make<ldg_expr> ("%dem_step", size));
			expr *scaled_re = gen_mult (sum.re, step);
			expr *scaled_im = gen_mult (sum.im, step);
			cplx_val scaled { scaled_re, scaled_im };
			zder.store (emit_complex_add (scaled, m));
		} else {
			cplx_val sum2 = emit_complex_sub (nrd, gen_mult_int (zder, power));
			cplx_val m = emit_complex_mult (sum2, cval);
			zder.store (m);
		}
	}

	zreg.store (newz);
}

// "Mixture" from the Saturn program: (c(2z - z^2) + q)^p

void gen_inner_e90mix (int size, int stepsize, int power, bool julia, bool dem,
		       cplx_reg zreg, cplx_val cval, cplx_reg zder)
{
	cplx_val zval { zreg };
	cplx_val twice = emit_complex_add (zval, zval);
	cplx_val pwr = emit_complex_mult (zval, zval);
	cplx_val diff = emit_complex_sub (twice, pwr);
	cplx_val cdiff = emit_complex_mult (cval, diff);
	auto pq = cplx_ldc ("const_param_q", size, stepsize);
	cplx_val sum = emit_complex_add (cdiff, pq);
	array<cplx_val, 20> powers;
	build_powers (powers, sum, power);
	cplx_val newz = get_power (powers, power);
	if (dem) {
		// Needs working out.
	}

	zreg.store (newz);
}

static void gen_inner (formula f, int size, int stepsize, int power, bool julia, bool dem,
		       cplx_reg zreg, cplx_reg creg, cplx_val cval, cplx_reg zder)
{
	switch (f) {
	default:
	case formula::standard:
		gen_inner_standard (size, power, julia, dem, zreg, cval, zder);
		break;
	case formula::tricorn:
		gen_inner_tricorn (size, power, julia, zreg, cval);
		break;
	case formula::ship:
		gen_inner_ship (power, julia, zreg, cval);
		break;
	case formula::celtic:
		gen_inner_celtic (power, julia, zreg, cval);
		break;
	case formula::spider:
		gen_inner_spider (size, power, zreg, creg);
		break;
	case formula::lambda:
		gen_inner_lambda (size, power, julia, dem, zreg, cval, zder);
		break;
	case formula::mix:
		gen_inner_mix (size, stepsize, power, julia, dem, zreg, cval, zder);
		break;
	case formula::e90_mix:
		gen_inner_e90mix (size, stepsize, power, julia, dem, zreg, cval, zder);
		break;
	case formula::magnet_a:
		gen_inner_magnet_a (size, power, julia, zreg, cval);
		break;
	case formula::facing:
		gen_inner_facing (size, power, julia, zreg, cval);
		break;
	case formula::facing_b:
		gen_inner_facing_b (size, power, julia, zreg, cval);
		break;
	case formula::rings:
		gen_inner_rings (size, power, julia, zreg, cval);
		break;
	case formula::tails:
		gen_inner_tails (size, stepsize, power, julia, zreg, cval);
		break;
	}
}

static void gen_inner_hybrid (QString &result, generator &cg,
			      formula f, int size, int stepsize, int power, bool julia,
			      cplx_reg zreg, cplx_reg creg, cplx_val cval)
{
	result += "\t.reg.u32 %hval;\n";
	result += "\t.reg.pred %hpred;\n";
	result += "\tand.b32\t%hval, %hybrid_code, %hybrid_mask;\n";
	result += "\tshl.b32\t%hybrid_code, %hybrid_code, 1;\n";
	result += "\tsetp.ne.u32\t%hpred, %hval, 0;\n";
	result += "@%hpred\tbra.uni\tstditer;\n";

	gen_inner (f, size, stepsize, power, julia, false, zreg, creg, cval, cplx_reg ());
	result += cg.code ();

	result += "\tbra\tloopend;\n";
	result += "stditer:\n";
	result += "\tadd.u32\t%hybrid_code, %hybrid_code, 1;\n";
	gen_inner_standard (size, power, julia, false, zreg, cval, cplx_reg ());
	result += cg.code ();
	result += "loopend:\n";
}

static void gen_store_zprev (QString &result, real_reg &zr, real_reg &zi)
{
	QString fzre = convert_to_float (zr);
	QString fzim = convert_to_float (zi);
	result += codegen->code ();

	result += "\tsub.u32\t\t%zpidx, %zpidx, 1;\n";
	result += QString ("\tand.b32\t\t%zpidx, %zpidx, %nprev_mask;\n");

	result += "\tcvt.u64.u32\t\t%prevaddr, %zpidx;\n";
	result += "\tmul.lo.u64\t\t%prevaddr, %prevaddr, 16;\n";
	result += "\tadd.u64\t\t%prevaddr, %prevaddr, %ar_zprev;\n";
	result += QString ("\tst.f64\t[%prevaddr], %1;\n").arg (fzre);
	result += QString ("\tst.f64\t[%prevaddr + 8], %1;\n").arg (fzim);
}

static void gen_compare_and_branch (QString &result, real_val &r1, real_val &r2, const QString &predreg, const QString &tmpreg, const QString &neq_label, const QString distance)
{
	int len = r1.ex ()->length ();
	for (int i = 0; i < len; i++) {
		QString p1 = r1.ex ()->get_piece_high (i);
		QString p2 = r2.ex ()->get_piece_high (i);
		if (i + 1 < len) {
			result += QString ("\tsetp.ne.u32\t%1, %2, %3;\n").arg (predreg, p1, p2);
		} else {
			result += QString ("\tsub.u32\t%1, %2, %3;\n").arg (tmpreg, p1, p2);
			result += QString ("\tabs.s32\t%1, %1;\n").arg (tmpreg);
			result += QString ("\tsetp.ge.u32\t%1, %2, %3;\n").arg (predreg, tmpreg, distance);
		}
		result += QString ("@%1\tbra\t\t%2;\n").arg (predreg, neq_label);
	}
}

static void gen_test_and_branch (QString &result, real_val &r1, const QString &predreg, const QString &tmpreg, const QString &neq_label, const QString distance, int maxlen = 0)
{
	int len = maxlen != 0 ? maxlen : r1.ex ()->length ();
	for (int i = 0; i < len; i++) {
		QString p1 = r1.ex ()->get_piece_high (i);
		if (i + 1 < len) {
			result += QString ("\tsetp.ne.u32\t%1, %2, 0;\n").arg (predreg, p1);
		} else {
			result += QString ("\tsetp.ge.u32\t%1, %2, %3;\n").arg (predreg, p1, distance);
		}
		result += QString ("@%1\tbra\t\t%2;\n").arg (predreg, neq_label);
	}
}

static void gen_compare_lt_and_branch (QString &result, real_val &r1, real_val &r2, const QString &predreg, const QString &lt_label, const QString &nlt_label)
{
	int len = r1.ex ()->length ();
	for (int i = 0; i < len; i++) {
		QString p1 = r1.ex ()->get_piece_high (i);
		QString p2 = r2.ex ()->get_piece_high (i);
		result += QString ("\tsetp.lt.u32\t%1, %2, %3;\n").arg (predreg, p1, p2);
		result += QString ("@%1\tbra\t\t%2;\n").arg (predreg, lt_label);
		if (i + 1 < len) {
			result += QString ("\tsetp.ne.u32\t%1, %2, %3;\n").arg (predreg, p1, p2);
			result += QString ("@%1\tbra\t\t%2;\n").arg (predreg, nlt_label);
		} else
			result += QString ("bra\t\t%1;\n").arg (nlt_label);
	}
}

void gen_kernel (formula f, QString &result, int size, int stepsize, int power, int pwrb, bool incolor,
		 bool julia, bool dem, bool hybrid = false)
{
	QString nm = julia ? "iter_julia" : "iter_mandel";
	if (dem)
		nm += "_dem";
	if (hybrid)
		nm += "_hybrid";
	gen_kernel_header (result, nm,
			   "u64", "ar_z", "u64", "ar_zprev", "u64", "ar_intvals",
			   "u64", "ar_coords", "u64", "ar_step",
			   "u32", "maxidx", "u64", "ar_result", "u32", "count", "u32", "init",
			   "u32", "hybrid_code", "u32", "hybrid_mask");

	int n_rvals = n_formula_real_vals (f, dem, incolor);
	int n_ivals = n_formula_int_vals (f, dem, incolor);
	int n_dvals = n_formula_extra_doubles (f, dem, incolor);

	generator cg;
	codegen = &cg;

	real_reg cr (size, "cr", true);
	real_reg ci (size, "ci", true);
	real_reg zr (size, "zr", true, true);
	real_reg zi (size, "zi", true, true);
	real_reg zderr = dem ? real_reg (size, "zderr", true) : real_reg ();
	real_reg zderi = dem ? real_reg (size, "zderi", true) : real_reg ();

	cplx_reg zreg (zr, zi);
	cplx_reg zder (zderr, zderi);
	cplx_reg creg (cr, ci);
	cplx_val cval;

	result += cg.code ();

	QString kernel_init = R"(
	.reg.s32 %idx, %tidx, %ctaidx, %ntidx, %n_prev, %nprev_mask, %npoff;
	ld.const.u32	%n_prev, [const_nprev];
	sub.u32		%nprev_mask, %n_prev, 1;

	.reg.u32 %iter, %scratch, %scratchoff;
	mov.u32         %ntidx, %ntid.x;
	mov.u32         %ctaidx, %ctaid.x;
	mov.u32         %tidx, %tid.x;
	mad.lo.s32      %idx, %ctaidx, %ntidx, %tidx;
	.reg.pred	%p1;
	setp.ge.s32     %p1, %idx, %maxidx;
@%p1	exit;
	.reg.pred	%pinit;
	setp.lt.s32     %pinit, %idx, %init;

	.reg.u64	%addroff, %iaddroff, %ar_zpidx, %ar_extra;
	add.u32		%npoff, %n_prev, %5;
	add.u32		%npoff, %npoff, %n_prev;
	mul.lo.u32	%npoff, %npoff, %idx;
	mul.lo.u32	%npoff, %npoff, 8;
	mul.lo.u32	%idx, %idx, 4;
	cvt.u64.u32	%addroff, %idx;
	add.u64		%ar_result, %ar_result, %addroff;
	add.u64		%ar_coords, %ar_coords, %addroff;
	cvt.u64.u32	%iaddroff, %idx;
	mul.lo.u64	%iaddroff, %iaddroff, %4;
	add.u64		%ar_zpidx, %ar_intvals, %iaddroff;
	mul.lo.u64	%addroff, %addroff, %2;
	add.u64		%ar_z, %ar_z, %addroff;
	cvt.u64.u32	%addroff, %npoff;
	add.u64		%ar_zprev, %ar_zprev, %addroff;
	cvt.u64.u32	%addroff, %n_prev;
	mul.lo.u64	%addroff, %addroff, 16;
	add.u64		%ar_extra, %ar_zprev, %addroff;

	mov.u32		%scratch, shared;
	mul.lo.u32	%scratchoff, %tidx, %3;
	add.u32		%scratch, %scratch, %scratchoff;

	.reg.u64	%ar_zim, %ar_t, %ar_tim;
	add.u64		%ar_zim, %ar_z, %1;
	add.u64		%ar_t, %ar_zim, %1;
	add.u64		%ar_tim, %ar_t, %1;
)";
	result += kernel_init.arg (size * 4).arg (size * n_rvals).arg (4 * size + 4).arg (n_ivals).arg (n_dvals);

	QString last_elt = "%ar_tim";
	if (dem) {
		result += QString (R"(
	.reg.u64	%ar_zder, %ar_zderim, %dem_step;
	add.u64		%ar_zder, %3, %1;
	add.u64		%ar_zderim, %ar_zder, %1;
	add.u64		%dem_step, %ar_step, %2;
)").arg (size * 4).arg (stepsize * 4 - size * 4).arg (last_elt);
		last_elt = "%ar_zderim";
	}
	if (incolor) {
		result += QString (R"(
	.reg.u64	%ar_zmin, %ar_iters, %ar_period;
	add.u64		%ar_zmin, %2, %1;
	add.u64		%ar_period, %ar_zpidx, 4;
	add.u64		%ar_iters, %ar_period, 4;
)").arg (size * 4).arg (last_elt);
		last_elt = "%ar_zmin";
	}

	result += R"(
	.reg.s32	%xpos, %ypos;
	.reg.s16	%xplow;
	ld.global.s32	%xpos, [%ar_coords];
	shr.s32		%ypos, %xpos, 16;
	cvt.s16.s32	%xplow, %xpos;
	cvt.s32.s16	%xpos, %xplow;

@%pinit bra		notfirst;
	call		coord_mul, (%ar_t, %ar_step, %xpos);
	call		coord_mul, (%ar_tim, %ar_step, %ypos);
)";

	real_val incoord_x (make<ldg_expr> ("%ar_t", size));
	real_val incoord_y (make<ldg_expr> ("%ar_tim", size));
	real_val originx (make<ldc_expr> (QString ("const_origin_x + %1")
					  .arg (stepsize * 4 - size * 4), size));
	real_val originy (make<ldc_expr> (QString ("const_origin_y + %1")
					  .arg (stepsize * 4 - size * 4), size));
	real_val mat00 (make<ldc_expr> (QString ("const_matrix00 + %1")
					.arg (stepsize * 4 - size * 4), size));
	real_val mat01 (make<ldc_expr> (QString ("const_matrix01 + %1")
					.arg (stepsize * 4 - size * 4), size));
	real_val mat10 (make<ldc_expr> (QString ("const_matrix10 + %1")
					.arg (stepsize * 4 - size * 4), size));
	real_val mat11 (make<ldc_expr> (QString ("const_matrix11 + %1")
					.arg (stepsize * 4 - size * 4), size));

	auto coord_x1 = make<addsub_expr> ("add", gen_mult (incoord_x, mat00), gen_mult (incoord_y, mat01));
	auto coord_y1 = make<addsub_expr> ("add", gen_mult (incoord_x, mat10), gen_mult (incoord_y, mat11));
	auto coord_x = make<addsub_expr> ("add", originx, coord_x1);
	auto coord_y = make<addsub_expr> ("add", originy, coord_y1);

	coord_x->calculate_full ();
	coord_y->calculate_full ();
	gen_store ("%ar_t", coord_x);
	gen_store ("%ar_tim", coord_y);

	auto const0 = make<const_expr<0>> (size);
	auto const1 = make<const_expr<1>> (size);
	if (dem) {
		if (julia)
			gen_store ("%ar_zder", make<ldg_expr> ("%dem_step", size));
		else
			gen_store ("%ar_zder", const0);
		gen_store ("%ar_zderim", const0);
	}

	expr *critr = nullptr;
	expr *criti = nullptr;
	if (!julia && formula_zero_critpoint (f)) {
		if (f == formula::facing) {
			gen_store ("%ar_z", coord_x);
			gen_store ("%ar_zim", coord_y);
		} else {
			gen_store ("%ar_z", const0);
			gen_store ("%ar_zim", const0);
		}
	} else if (!julia) {
		critr = make<ldc_expr> (QString ("const_critpoint + %1")
					.arg (stepsize * 4 - size * 4), size);
		criti = make<ldc_expr> (QString ("const_critpoint + %1")
					.arg (2 * stepsize * 4 - size * 4), size);
		if (f == formula::facing) {
			gen_store ("%ar_z", coord_x);
			gen_store ("%ar_zim", coord_y);
		} else {
			gen_store ("%ar_z", critr);
			gen_store ("%ar_zim", criti);
		}
	} else {
		if (f == formula::spider) {
			cplx_val parm_p = cplx_ldc ("const_param_p", size, stepsize);
			gen_store ("%ar_t", parm_p.re);
			gen_store ("%ar_tim", parm_p.im);
		}
		gen_store ("%ar_z", coord_x);
		gen_store ("%ar_zim", coord_y);
	}
	auto do_setup_loads = [ f, &cr, &ci, &zr, &zi, &zderr, &zderi, &cval, &creg, size, stepsize, dem, julia ] ()
	{
		cr.store (make<ldg_expr> ("%ar_t", size));
		ci.store (make<ldg_expr> ("%ar_tim", size));
		zr.store (make<ldg_expr> ("%ar_z", size));
		zi.store (make<ldg_expr> ("%ar_zim", size));
		if (dem) {
			zderr.store (make<ldg_expr> ("%ar_zder", size));
			zderi.store (make<ldg_expr> ("%ar_zderim", size));
		}
		cval = creg;
		if (julia)
			cval = cplx_ldc ("const_param_p", size, stepsize);
		if (f == formula::facing) {
			cval = emit_complex_sqr (cval);
		}
	};
	auto do_finish_stores = [ f, &cr, &ci, &zr, &zi, &zderr, &zderi, dem ] ()
	{
		gen_store ("%ar_z", zr);
		gen_store ("%ar_zim", zi);

		if (f == formula::spider) {
			gen_store ("%ar_t", cr);
			gen_store ("%ar_tim", ci);
		}
		if (dem) {
			gen_store ("%ar_zder", zderr);
			gen_store ("%ar_zderim", zderi);
		}
	};
	if (incolor) {
		result += "\tst.global.u32\t[%ar_iters], 0;\n";
		result += "\tst.global.u32\t[%ar_period], 0;\n";
		gen_store ("%ar_zmin", make<const_expr<100000>> (size));
	}
	result += cg.code ();
	result += "\tst.global.u32\t[%ar_zpidx], 0;\n";

	if (pwrb > 1 && julia) {
		do_setup_loads ();
		gen_inner (f, size, stepsize, pwrb, julia, dem, zreg, creg, cval, zder);
		do_finish_stores ();
		result += cg.code ();
	}

	result += "notfirst:\n";
	result += "\t.reg.u32\t\t%zpidx;\n";
	result += "\tld.global.u32\t%zpidx, [%ar_zpidx];\n";

	do_setup_loads ();

	result += cg.code ();

	bool fixpoints = formula_test_fixpoint (f);
	real_reg zlastr = fixpoints ? real_reg (size, "zlastr", true) : real_reg ();
	real_reg zlasti = fixpoints ? real_reg (size, "zlasti", true) : real_reg ();
	cplx_reg zlastreg (zlastr, zlasti);

	real_reg zmin = incolor ? real_reg (size, "zmin") : real_reg ();
	if (incolor) {
		result += "\t.reg.u32\t\t%iters;\n";
		result += "\tld.global.u32\t%iters, [%ar_iters];\n";
		zmin.store (make<ldg_expr> ("%ar_zmin", size));
	}

	result += cg.code ();

	result += R"(

	.reg.u32 %niter;
	mov.u32	%niter, 0;
	.reg.pred	%using_zprev;
	setp.gt.u32	%using_zprev, %n_prev, 1;
loop:
	.reg.u64	%prevaddr;
	add.u32		%niter, %niter, 1;
@!%using_zprev bra	skip_zprev;
)";

	gen_store_zprev (result, zr, zi);
	result += cg.code ();
	result += "skip_zprev:\n";

	if (fixpoints)
		zlastreg.store (zreg);

	if (hybrid)
		gen_inner_hybrid (result, cg, f, size, stepsize, power, julia, zreg, creg, cval);
	else
		gen_inner (f, size, stepsize, power, julia, dem, zreg, creg, cval, zder);

	result += cg.code ();

	if (incolor) {
		result += "\t.reg.pred\t%early, %uneq;\n";
		result += "\t.reg.u32\t%thisperiod, %dist_tmp;\n";
		result += "\tadd.u32\t\t%iters, %iters, 1;\n";
		real_val zroff { critr == nullptr ? zr : make<addsub_expr> ("sub", zr, critr) };
		real_val zioff { criti == nullptr ? zi : make<addsub_expr> ("sub", zi, criti) };
		real_val dist { make<addsub_expr> ("add", zroff.squared (), zioff.squared ()) };
		dist.ex ()->calculate_full ();
		result += cg.code ();
		real_val ziv { zi };
		real_val zminv { zmin };
		gen_compare_lt_and_branch (result, dist, zminv, "%uneq", "min_found", "not_new_min");
		result += "min_found:\n";
		zmin.store (dist);
		QString fmin = convert_to_float (dist);
		result += cg.code ();
		result += QString ("\tst.f64\t[%ar_extra], %1;\n").arg (fmin);
		result += "\tst.global.u32\t[%ar_period], %iters;\n";
		result += "not_new_min:\n";
	}

	if (f == formula::magnet_a) {
		/* This is a hack for Magnet A, which is currently the only formula that needs this
		   check: both infinity and (1,0) are attractors.
		   The first attempt was to test the difference between zlast and z, but that caused
		   artifacts.  */
		result += "\t.reg.u32\t%diff;\n";
		result += "\t.reg.pred\t%fp_neq;\n";
		auto constm1 = make<const_expr<-1>> (size);
		real_val zrm1 (make<addsub_expr> ("add", zr, constm1));
		real_val dist (make<addsub_expr> ("add", gen_mult (zrm1, zrm1), zi.squared ()));
		dist.ex ()->calculate_full ();
		result += cg.code ();
		gen_test_and_branch (result, dist, "%fp_neq", "%diff", "skip2", "16", 2);
		result += "\tbra\t\tbailout;\n";
	} else if (fixpoints) {
		result += "\t.reg.u32\t%diff;\n";
		result += "\t.reg.pred\t%fp_neq;\n";
		real_val zrdelta (make<addsub_expr> ("sub", zr, zlastr));
		real_val zidelta (make<addsub_expr> ("sub", zi, zlasti));
		real_val dist (make<addsub_expr> ("add", zrdelta.squared (), zidelta.squared ()));
		dist.ex ()->calculate_full ();
		result += cg.code ();
		gen_test_and_branch (result, dist, "%fp_neq", "%diff", "skip2", "16");
		result += "\tbra\t\tbailout;\n";
	}
	QString z2r_high = zr.squared ()->get_piece_high (0);
	QString z2i_high = zi.squared ()->get_piece_high (0);
	QString loop_end = R"(
skip2:
	.reg.u32 %sqsum;
	add.u32		%sqsum, %1, %2;
	.reg.pred	%cont;
	setp.lt.u32	%cont, %sqsum, %3;
@%cont	bra		skip;
bailout:
)";
	result += loop_end.arg (z2r_high, z2i_high, dem ? "100" : "10000");

	result += "\tst.global.u32\t[%ar_result], %niter;\n";

	gen_store_zprev (result, zr, zi);
	result += cg.code ();

	result += R"(
	bra		store_and_exit;
skip:
	sub.u32		%count, %count, 1;
	.reg.pred	%again;
	setp.gt.u32	%again, %count, 0;
@%again	bra		loop;

	st.global.u32	[%ar_result], 0;

store_and_exit:
	st.global.u32	[%ar_zpidx], %zpidx;
)";

	if (incolor) {
		result += "\tst.global.u32\t[%ar_iters], %iters;\n";
		gen_store ("%ar_zmin", zmin);
	}
	do_finish_stores ();
	result += cg.code ();
	result += "}\n";
}

char *gen_mprec_funcs (formula f, int size, int stepsize, int power, int pwrb, bool incolor)
{
	QString result;
	result += QString (R"(	.version	6.2
	.target	sm_61
	.address_size 64
	.extern .shared .align 4 .b8 shared[];
	.const .align 4 .u32 const_origin_x[%1];
	.const .align 4 .u32 const_origin_y[%1];
	.const .align 4 .u32 const_param_p[%1];
	.const .align 4 .u32 const_param_q[%1];
	.const .align 4 .u32 const_matrix00[%1];
	.const .align 4 .u32 const_matrix01[%1];
	.const .align 4 .u32 const_matrix10[%1];
	.const .align 4 .u32 const_matrix11[%1];
	.const .align 4 .u32 const_critpoint[%1];
	.const .align 4 .u32 const_nprev;
)").arg (stepsize * 2);

#if 0
	gen_mul_func (result, size, false);
	gen_mul_func (result, size, true);
#endif
	gen_coord_muladd (result, size, stepsize);
	gen_coord_mul (result, size, stepsize);

	gen_kernel (f, result, size, stepsize, power, pwrb, incolor, true, false);
	gen_kernel (f, result, size, stepsize, power, pwrb, incolor, false, false);
	if (formula_supports_dem (f)) {
		gen_kernel (f, result, size, stepsize, power, pwrb, incolor, true, true);
		gen_kernel (f, result, size, stepsize, power, pwrb, incolor, false, true);
	} else if (formula_supports_hybrid (f)) {
		gen_kernel (f, result, size, stepsize, power, pwrb, incolor, true, false, true);
		gen_kernel (f, result, size, stepsize, power, pwrb, incolor, false, false, true);
	}

	QByteArray a = result.toLatin1 ();
#if 0
	std::cerr << a.constData ();
#endif
	return strdup (a.constData ());
}
