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
using std::shared_ptr;
using std::make_shared;

class generator
{
	int m_tmp = 0;
	int m_lbl_tmp = 0;
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
};

generator *codegen;

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

class addsub_expr : public expr
{
	QString m_op, m_cond;
	shared_ptr<expr> m_a, m_b;
	bool nonzero_added = false;
public:
	addsub_expr (const QString &op, shared_ptr<expr> srca, shared_ptr<expr> srcb, QString cond = QString ()) :
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
	shared_ptr<expr> m_op;
public:
	arshift_expr (shared_ptr<expr> op) :
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
	shared_ptr<expr> m_op;
	int m_shift;
public:
	lshift_expr (shared_ptr<expr> op, int shift) :
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
	shared_ptr<expr> m_op;
public:
	padlow_expr (shared_ptr<expr> src, int len)
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
	shared_ptr<expr> m_op;
public:
	padhigh_expr (shared_ptr<expr> src, int len)
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
	shared_ptr<expr> m_op;
public:
	trunc_expr (shared_ptr<expr> src, int len)
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
	shared_ptr<expr> m_low, m_high;
public:
	concat_expr (shared_ptr<expr> op_low, shared_ptr<expr> op_high)
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
	shared_ptr<expr> m_op;
	int m_off;

public:
	lowpart_expr (shared_ptr<expr> src, int off, int len)
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
	shared_ptr<expr> m_op;
	QString m_pred, m_pred_reg;
public:
	abs_expr (shared_ptr<expr> src, QString pred = QString ())
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
		gen_cond_negate (&*this, &*m_op, p);
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
	shared_ptr<expr> m_op;
	QString m_sgn1, m_sgn2;

public:
	mult_sign_fixup_expr (shared_ptr<expr> op, shared_ptr<abs_expr> sgn1, shared_ptr<abs_expr> sgn2)
		: expr (op->length ()), m_op (op), m_sgn1 (sgn1->get_pred ()), m_sgn2 (sgn2->get_pred ())
	{
	}
	mult_sign_fixup_expr (shared_ptr<expr> op, const QString &pred1, const QString &pred2)
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
		gen_cond_negate (this, &*m_op, pred);
		return m_values[0];
	}
};


class mult_expr : public expr
{
	int m_parts_len;
	shared_ptr<expr> m_a, m_b;
	QString m_carry, m_carry2;
	vector<QString> m_preds_a, m_preds_b;
	bool m_skip_ones;
public:
	mult_expr (shared_ptr<expr> srca, shared_ptr<expr> srcb, int discard_top = 1, bool skip_ones = false) :
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
	shared_ptr<expr> m_a, m_b, m_tmpm;

public:
	mult_special_expr (shared_ptr<expr> srca, shared_ptr<expr> srcb, int discard_top = 1) :
		expr (2 * std::max (srca->length (), srcb->length ()) - discard_top), m_a (srca), m_b (srcb)
	{
		srca->require_carry ();
		srcb->require_carry ();
		m_tmpm = make_shared<mult_expr> (m_a, m_b, discard_top, true);

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

void gen_store (const QString &dst, shared_ptr<expr> ex)
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

void gen_store (reg_expr *dst, shared_ptr<expr> ex)
{
	int len = ex->length ();
	for (int i = 0; i < len; i++) {
		QString dstp = dst->get_piece (i);
		QString val = ex->get_piece (i);
		codegen->append_move ("u32", dstp, val);
	}
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
shared_ptr<expr> gen_mult_karatsuba_1 (shared_ptr<expr> a, shared_ptr<expr> b)
{
	if (a->length () < karatsuba_cutoff) {
		return make_shared<mult_special_expr> (a, b, 0);
	}
	int half = (a->length () + 1) / 2;
	int full = 2 * half;

	auto ah = make_shared<trunc_expr> (a, half);
	auto bh = a == b ? ah : make_shared<trunc_expr> (b, half);
	auto al = make_shared<lowpart_expr> (a, half, half);
	auto bl = a == b ? al : make_shared<lowpart_expr> (b, half, half);

	auto H1 = make_shared<addsub_expr> ("add", make_shared<padhigh_expr> (ah, half + 1),
					    make_shared<padhigh_expr> (al, half + 1));
	auto H2 = a == b ? H1 : make_shared<addsub_expr> ("add", make_shared<padhigh_expr> (bh, half + 1),
							  make_shared<padhigh_expr> (bl, half + 1));
	auto H = make_shared<lowpart_expr> (gen_mult_karatsuba_1 (H1, H2), 1, full + 1);
	auto F = gen_mult_karatsuba_1 (ah, bh);
	auto G = gen_mult_karatsuba_1 (al, bl);
	assert (F->length () == full);
	assert (G->length () == full);
	assert (H->length () == full + 1);
	F->set_destreg ("F");
	G->set_destreg ("G");
	H->set_destreg ("H");
	auto sub1 = make_shared<addsub_expr> ("sub", H, make_shared<padhigh_expr> (G, full + 1));
	auto sub2 = make_shared<addsub_expr> ("sub", sub1, make_shared<padhigh_expr> (F, full + 1));
	sub2->set_destreg ("K");
	auto highlow = make_shared<concat_expr> (G, F);
	auto Ke = make_shared<padlow_expr> (make_shared<padhigh_expr> (sub2, full + half), 2 * full);
	auto sum = make_shared<addsub_expr> ("add", highlow, Ke);
	sum->set_destreg ("msum");
	if (half * 2 > a->length ())
		return make_shared<trunc_expr> (sum, 2 * a->length ());
	return sum;
}

shared_ptr<expr> gen_mult_unsigned (shared_ptr<expr> a, shared_ptr<expr> b)
{
	if (a->length () < karatsuba_cutoff)
		return make_shared<trunc_expr> (make_shared<mult_expr> (a, b), a->length ());

	auto me = gen_mult_karatsuba_1 (a, b);
	return make_shared<trunc_expr> (make_shared<lowpart_expr> (me, 1, me->length () - 1), a->length ());
}

struct real_val;

/* A structure for information about a real-valued reg, and optionally some of its
   related values: abs_val and neg_pred can hold a sign/magnitude representation,
   which is often useful, and square_val can be used to keep the square of the
   value around.  */
class real_reg
{
	shared_ptr<reg_expr> m_reg;
	shared_ptr<reg_expr> m_abs_val;
	shared_ptr<reg_expr> m_square_val;
	QString m_neg_pred;
	friend class real_val;

public:
	real_reg () = default;
	real_reg (int size, const QString &name, bool keep_abs = false, bool keep_square = false)
		: m_reg (make_shared<reg_expr> (size, name)),
		  m_abs_val (keep_abs || keep_square ? make_shared<reg_expr> (size, name + "a") : nullptr),
		  m_square_val (keep_square ? make_shared<reg_expr> (size, name + "sq") : nullptr),
		  m_neg_pred (keep_abs ? codegen->gen_reg ("pred", name + "n") : nullptr)
	{
	}
	void store (shared_ptr<expr> val)
	{
		gen_store (&*m_reg, val);
		if (!m_neg_pred.isEmpty () || m_square_val != nullptr) {
			auto absexp = make_shared<abs_expr> (val, m_neg_pred);
			if (!m_neg_pred.isEmpty ())
				gen_store (&*m_abs_val, absexp);
			if (m_square_val != nullptr)
				gen_store (&*m_square_val, gen_mult_unsigned (absexp, absexp));
		}
	}
	operator shared_ptr<expr> () const
	{
		return m_reg;
	}
	shared_ptr<expr> ex () const
	{
		return m_reg;
	}
	shared_ptr<expr> squared ()
	{
		if (m_square_val == nullptr) {
			if (m_abs_val != nullptr)
				return gen_mult_unsigned (m_abs_val, m_abs_val);
			auto absv = abs_val ();
			return gen_mult_unsigned (absv, absv);
		}
		return m_square_val;
	}
	shared_ptr<expr> abs_val ()
	{
		if (m_abs_val == nullptr)
			return make_shared<abs_expr> (m_reg, m_neg_pred);
		return m_abs_val;
	}
	QString neg_pred ()
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
	shared_ptr<expr> m_val;
	shared_ptr<expr> m_abs_val;
	shared_ptr<expr> m_square_val;
	QString m_neg_pred;

public:
	real_val (const real_reg &reg) : m_val (reg.m_reg), m_abs_val (reg.m_abs_val), m_square_val (reg.m_square_val), m_neg_pred (reg.m_neg_pred)
	{
	}
	real_val (shared_ptr<expr> expr) : m_val (expr)
	{
	}
	real_val (shared_ptr<expr> expr, const QString &pred) : m_abs_val (expr), m_neg_pred (pred)
	{
	}
	real_val () = default;
	operator shared_ptr<expr> () const
	{
		return ex ();
	}
	shared_ptr<expr> ex () const
	{
		if (m_val != nullptr)
			return m_val;
		shared_ptr<expr> dst = make_shared<reg_expr> (m_abs_val->length (), "sgn");
		gen_cond_negate (&*dst, &*m_abs_val, m_neg_pred);
		return dst;
	}
	shared_ptr<expr> squared ()
	{
		if (m_square_val == nullptr) {
			auto absv = abs_val ();
			m_square_val = gen_mult_unsigned (absv, absv);
		}
		return m_square_val;
	}
	shared_ptr<expr> abs_val ()
	{
		if (m_abs_val == nullptr)
			m_abs_val = make_shared<abs_expr> (m_val, m_neg_pred);
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
real_val gen_mult (R1 &a, R2 &b, bool abs_required = true)
{
	shared_ptr<expr> aop = a.abs_val ();
	shared_ptr<expr> bop = b.abs_val ();
	auto prod = gen_mult_unsigned (aop, bop);
	if (!abs_required || &a == &b)
		return prod;

	QString pred = codegen->gen_reg ("pred", "sfix");
	codegen->append_code (QString ("\txor.pred\t%1, %2, %3;\n").arg (pred, a.neg_pred (), b.neg_pred ()));
	return real_val (prod, pred);
}

QString convert_to_float (shared_ptr<expr> v)
{
	QString dstreg = codegen->gen_reg ("f64", "convdst");
	int len = v->length ();
	QString conv = codegen->gen_reg ("f64", "fltconv");
	QString mult = codegen->gen_reg ("f64", "multiplier");
	codegen->append_move ("f64", conv, "1.0");
	codegen->append_code (QString ("\tdiv.rn.f64\t%1, %1, 4294967296.0;\n").arg (conv));

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
	auto dst = make_shared<reg_expr> (size, "dst", false);
	auto srca = make_shared<reg_expr> (size, "srca", false);
	auto srcb = sqr ? srca : make_shared<reg_expr> (size, "srcb", false);
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

	gen_store (&*dst, gen_mult (srca, srcb));
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
	cplx_val (shared_ptr<expr> r, shared_ptr<expr> i) : re (r), im (i)
	{
	}
	cplx_val (const cplx_reg &reg) : re (reg.re), im (reg.im)
	{
	}
	void set_re (shared_ptr<expr> v)
	{
		re = real_val (v);
	}
	void set_im (shared_ptr<expr> v)
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

template<class C>
cplx_val emit_complex_sqr (C &v)
{
	auto sqre = v.re.squared ();
	auto sqim = v.im.squared ();

	// zr <- z2r - z2i
	// zi <- zr * zi * 2

	auto tmp = gen_mult (v.re, v.im);

	auto re = make_shared<addsub_expr> ("sub", sqre, sqim);
	auto im = make_shared<addsub_expr> ("add", tmp, tmp);
	return { re, im };
}

cplx_val emit_complex_mult (cplx_val &va, cplx_val &vb)
{
	// (a+ib) * (c+id) = ac - bd + i(bc + ad)
	real_val rp1_expr = gen_mult (va.re, vb.re);
	real_val rp2_expr = gen_mult (va.im, vb.im);

	real_val ip1_expr = gen_mult (va.re, vb.im);
	real_val ip2_expr = gen_mult (va.im, vb.re);

	auto re = make_shared<addsub_expr> ("sub", rp1_expr, rp2_expr);
	auto im = make_shared<addsub_expr> ("add", ip1_expr, ip2_expr);
	return { re, im };
}

template<class C>
cplx_val gen_mult_int (C &v, int factor)
{
	cplx_val result;
	bool first = true;
	for (int i = 0; i < 10; i++)
		if (factor & (1 << i)) {
			shared_ptr<expr> shfr = i == 0 ? v.re.ex () : make_shared<lshift_expr> (v.re, i);
			shared_ptr<expr> shfi = i == 0 ? v.im.ex () : make_shared<lshift_expr> (v.im, i);
			if (first) {
				result = { shfr, shfi };
				first = false;
			} else {
				auto re = make_shared<addsub_expr> ("add", result.re, shfr);
				auto im = make_shared<addsub_expr> ("add", result.im, shfi);
				result = { re, im };
			}
		}
	return result;
}

template<class C>
void build_powers (array<cplx_val, 20> &powers, C &z, int power)
{
	powers[0] = z;
	powers[1] = emit_complex_sqr (z);
	int pidx = 2;
	for (int p = 4; p <= power; p *= p, pidx++) {
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
	auto re = make_shared<ldc_expr> (QString ("%1 + %2")
					 .arg (name).arg (stepsize * 4 - size * 4), size);
	auto im = make_shared<ldc_expr> (QString ("%1 + %2")
					 .arg (name).arg (2 * stepsize * 4 - size * 4), size);
	return { re, im };
}

cplx_val emit_complex_add (const cplx_val &v1, const cplx_val &v2)
{
	shared_ptr<addsub_expr> newzr, newzi;
	newzr = make_shared<addsub_expr> ("add", v1.re, v2.re);
	newzi = make_shared<addsub_expr> ("add", v1.im, v2.im);
	return { newzr, newzi };
}

cplx_val emit_complex_sub (const cplx_val &v1, const cplx_val &v2)
{
	shared_ptr<addsub_expr> newzr, newzi;
	newzr = make_shared<addsub_expr> ("sub", v1.re, v2.re);
	newzi = make_shared<addsub_expr> ("sub", v1.im, v2.im);
	return { newzr, newzi };
}

void gen_inner_standard (int size, int power, bool julia, bool dem,
			 cplx_reg zreg, cplx_val cval, cplx_reg zder)
{
	array<cplx_val, 20> powers;
	build_powers (powers, zreg, power);
	cplx_val pwr = get_power (powers, power);

	cplx_val newz = emit_complex_add (pwr, cval);
	if (dem) {
		cplx_val dpwr = get_power (powers, power - 1);
		cplx_val zdval (zder);
		cplx_val nrd = emit_complex_mult (dpwr, zdval);
		nrd = gen_mult_int (nrd, power);
		if (!julia) {
			nrd.set_re (make_shared<addsub_expr> ("add", nrd.re,
							      make_shared<ldg_expr> ("%dem_step", size)));
		}
		zder.store (nrd);
	}
	zreg.store (newz);
}

// Burning ship, because why not set zre and zim to their absolute values before
// continuing.
void gen_inner_ship (int power, bool /* julia */, cplx_reg zreg, cplx_val cval)
{
	auto zar = make_shared<abs_expr> (zreg.re);
	auto zai = make_shared<abs_expr> (zreg.im);
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

	auto const0 = make_shared<const_expr<0>> (size);
	newz.set_im (make_shared<addsub_expr> ("sub", const0, newz.im));

	zreg.store (newz);
}

void gen_inner_celtic (int power, bool /* julia */, cplx_reg zreg, cplx_val cval)
{
	array<cplx_val, 20> powers;
	build_powers (powers, zreg, power);
	cplx_val pwr = get_power (powers, power);
	pwr.set_re (make_shared<abs_expr> (pwr.re));
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
			real_val step (make_shared<ldg_expr> ("%dem_step", size));
			shared_ptr<expr> scaled_re = gen_mult (sum.re, step);
			shared_ptr<expr> scaled_im = gen_mult (sum.im, step);
			cplx_val scaled { scaled_re, scaled_im };
			zder.store (emit_complex_add (scaled, m));
		} else {
			zder.store (m);
		}
	}

	zreg.store (newz);
}

void gen_inner_sqtwice_a (int power, bool /* julia */, cplx_reg zreg, cplx_val cval)
{
	array<cplx_val, 20> powers;
	build_powers (powers, zreg, power);
	cplx_val pwr = get_power (powers, power);
	cplx_val newz1 = emit_complex_add (pwr, cval);
	cplx_val newz = emit_complex_sqr (newz1);

	zreg.store (newz);
}


void gen_inner_sqtwice_b (int power, bool /* julia */, cplx_reg zreg, cplx_val cval)
{
	cplx_val pwr = emit_complex_sqr (zreg);
	cplx_val newz1 = emit_complex_add (pwr, cval);

	array<cplx_val, 20> powers_b;
	build_powers (powers_b, newz1, power);
	cplx_val newz = get_power (powers_b, power);

	zreg.store (newz);
}

void gen_inner_spider (int size, int power, cplx_reg zreg, cplx_reg creg)
{
	array<cplx_val, 20> powers;
	build_powers (powers, zreg, power);
	cplx_val pwr = get_power (powers, power);

	cplx_val newz = emit_complex_add (pwr, creg);

	auto prhalf = make_shared<arshift_expr<1>> (creg.re);
	auto pihalf = make_shared<arshift_expr<1>> (creg.im);
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
	auto const2 = make_shared<const_expr<2>> (size);
	sum.set_re (make_shared<addsub_expr> ("sub", sum.re, const2));
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
			real_val step (make_shared<ldg_expr> ("%dem_step", size));
			shared_ptr<expr> scaled_re = gen_mult (sum.re, step);
			shared_ptr<expr> scaled_im = gen_mult (sum.im, step);
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
	case formula::sqtwice_a:
		gen_inner_sqtwice_a (power, julia, zreg, cval);
		break;
	case formula::sqtwice_b:
		gen_inner_sqtwice_b (power, julia, zreg, cval);
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

static void gen_store_zprev (QString &result, real_reg &zr, real_reg &zi, int n_prev)
{
	QString fzre = convert_to_float (zr);
	QString fzim = convert_to_float (zi);
	result += codegen->code ();

	result += "\tsub.u32\t\t%zpidx, %zpidx, 1;\n";
	result += QString ("\tand.b32\t\t%zpidx, %zpidx, %1;\n").arg (n_prev - 1);

	result += "\tcvt.u64.u32\t\t%prevaddr, %zpidx;\n";
	result += "\tmul.lo.u64\t\t%prevaddr, %prevaddr, 16;\n";
	result += "\tadd.u64\t\t%prevaddr, %prevaddr, %ar_zprev;\n";
	result += QString ("\tst.f64\t[%prevaddr], %1;\n").arg (fzre);
	result += QString ("\tst.f64\t[%prevaddr + 8], %1;\n").arg (fzim);
}

void gen_kernel (formula f, QString &result, int size, int stepsize, int power, int n_prev,
		 bool julia, bool dem, bool hybrid = false)
{
	assert ((n_prev & (n_prev - 1)) == 0);
	QString nm = julia ? "iter_julia" : "iter_mandel";
	if (dem)
		nm += "_dem";
	if (hybrid)
		nm += "_hybrid";
	gen_kernel_header (result, nm,
			   "u64", "ar_z", "u64", "ar_zprev", "u64", "ar_zpidx",
			   "u64", "ar_coords", "u64", "ar_step",
			   "u32", "maxidx", "u64", "ar_result", "u32", "count", "u32", "init",
			   "u32", "hybrid_code", "u32", "hybrid_mask");

	QString kernel_init = R"(
	.reg.s32 %idx, %tidx, %ctaidx, %ntidx;

	.reg.u32 %iter;
	mov.u32         %ntidx, %ntid.x;
	mov.u32         %ctaidx, %ctaid.x;
	mov.u32         %tidx, %tid.x;
	mad.lo.s32      %idx, %ctaidx, %ntidx, %tidx;
	.reg.pred	%p1;
	setp.ge.s32     %p1, %idx, %maxidx;
@%p1	exit;
	.reg.pred	%pinit;
	setp.lt.s32     %pinit, %idx, %init;

	.reg.u64	%addroff, %addridx;
	cvt.u64.u32	%addridx, %idx;
	mul.lo.u32	%idx, %idx, 4;
	cvt.u64.u32	%addroff, %idx;
	add.u64		%ar_zpidx, %ar_zpidx, %addroff;
	add.u64		%ar_result, %ar_result, %addroff;
	add.u64		%ar_coords, %ar_coords, %addroff;
	mul.lo.u64	%addroff, %addroff, %2;
	add.u64		%ar_z, %ar_z, %addroff;
	mul.lo.u64	%addroff, %addridx, %3;
	add.u64		%ar_zprev, %ar_zprev, %addroff;

	.reg.u64	%ar_zim, %ar_t, %ar_tim;
	add.u64		%ar_zim, %ar_z, %1;
	add.u64		%ar_t, %ar_zim, %1;
	add.u64		%ar_tim, %ar_t, %1;
)";
	result += kernel_init.arg (size * 4).arg (size * 2 * n_formula_cplx_vals (f, dem)).arg (n_prev * 8 * 2);;

	if (dem) {
		result += QString (R"(
	.reg.u64	%ar_zder, %ar_zderim, %dem_step;
	add.u64		%ar_zder, %ar_tim, %1;
	add.u64		%ar_zderim, %ar_zder, %1;
	add.u64		%dem_step, %ar_step, %2;
)").arg (size * 4).arg (stepsize * 4 - size * 4);
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
	generator cg;
	codegen = &cg;

	real_val incoord_x (make_shared<ldg_expr> ("%ar_t", size));
	real_val incoord_y (make_shared<ldg_expr> ("%ar_tim", size));
	real_val originx (make_shared<ldc_expr> (QString ("const_origin_x + %1")
						 .arg (stepsize * 4 - size * 4), size));
	real_val originy (make_shared<ldc_expr> (QString ("const_origin_y + %1")
						 .arg (stepsize * 4 - size * 4), size));
	real_val mat00 (make_shared<ldc_expr> (QString ("const_matrix00 + %1")
					       .arg (stepsize * 4 - size * 4), size));
	real_val mat01 (make_shared<ldc_expr> (QString ("const_matrix01 + %1")
					       .arg (stepsize * 4 - size * 4), size));
	real_val mat10 (make_shared<ldc_expr> (QString ("const_matrix10 + %1")
					       .arg (stepsize * 4 - size * 4), size));
	real_val mat11 (make_shared<ldc_expr> (QString ("const_matrix11 + %1")
					       .arg (stepsize * 4 - size * 4), size));

	auto coord_x1 = make_shared<addsub_expr> ("add", gen_mult (incoord_x, mat00), gen_mult (incoord_y, mat01));
	auto coord_y1 = make_shared<addsub_expr> ("add", gen_mult (incoord_x, mat10), gen_mult (incoord_y, mat11));
	auto coord_x = make_shared<addsub_expr> ("add", originx, coord_x1);
	auto coord_y = make_shared<addsub_expr> ("add", originy, coord_y1);

	coord_x->calculate_full ();
	coord_y->calculate_full ();
	gen_store ("%ar_t", coord_x);
	gen_store ("%ar_tim", coord_y);

	auto const0 = make_shared<const_expr<0>> (size);
	auto const1 = make_shared<const_expr<1>> (size);
	if (dem) {
		if (julia)
			gen_store ("%ar_zder", make_shared<ldg_expr> ("%dem_step", size));
		else
			gen_store ("%ar_zder", const0);
		gen_store ("%ar_zderim", const0);
	}

	if (!julia) {
		auto critr = make_shared<ldc_expr> (QString ("const_critpoint + %1")
						    .arg (stepsize * 4 - size * 4), size);
		auto criti = make_shared<ldc_expr> (QString ("const_critpoint + %1")
						    .arg (2 * stepsize * 4 - size * 4), size);
		gen_store ("%ar_z", critr);
		gen_store ("%ar_zim", criti);
	} else {
		if (f == formula::spider) {
			cplx_val parm_p = cplx_ldc ("const_param_p", size, stepsize);
			gen_store ("%ar_t", parm_p.re);
			gen_store ("%ar_tim", parm_p.im);
		}
		gen_store ("%ar_z", coord_x);
		gen_store ("%ar_zim", coord_y);
	}
	result += cg.code ();
	result += "\tst.global.u32\t[%ar_zpidx], 0;\n";

	result += "notfirst:\n";
	result += "\t.reg.u32\t\t%zpidx;\n";
	result += "\tld.global.u32\t%zpidx, [%ar_zpidx];\n";

	real_reg cr (size, "cr", true);
	real_reg ci (size, "ci", true);
	real_reg zr (size, "zr", true, true);
	real_reg zi (size, "zi", true, true);

	cr.store (make_shared<ldg_expr> ("%ar_t", size));
	ci.store (make_shared<ldg_expr> ("%ar_tim", size));
	zr.store (make_shared<ldg_expr> ("%ar_z", size));
	zi.store (make_shared<ldg_expr> ("%ar_zim", size));
	real_reg zderr = dem ? real_reg (size, "zderr", true) : real_reg ();
	real_reg zderi = dem ? real_reg (size, "zderi", true) : real_reg ();
	if (dem) {
		zderr.store (make_shared<ldg_expr> ("%ar_zder", size));
		zderi.store (make_shared<ldg_expr> ("%ar_zderim", size));
	}

	result += cg.code ();

	cplx_reg zreg (zr, zi);
	cplx_reg zder (zderr, zderi);
	cplx_reg creg (cr, ci);
	cplx_val cval = creg;
	if (julia)
		cval = cplx_ldc ("const_param_p", size, stepsize);

	result += cg.code ();

	result += R"(

	.reg.u32 %niter;
	mov.u32	%niter, 0;

loop:
	.reg.u64	%prevaddr;
	add.u32		%niter, %niter, 1;

)";

	if (n_prev > 1)
		gen_store_zprev (result, zr, zi, n_prev);

	if (hybrid)
		gen_inner_hybrid (result, cg, f, size, stepsize, power, julia, zreg, creg, cval);
	else
		gen_inner (f, size, stepsize, power, julia, dem, zreg, creg, cval, zder);

	result += cg.code ();

	QString z2r_high = zr.squared ()->get_piece_high (0);
	QString z2i_high = zi.squared ()->get_piece_high (0);
	QString loop_end = R"(
	.reg.u32 %sqsum;
	add.u32		%sqsum, %1, %2;
	.reg.pred	%cont;
	setp.lt.u32	%cont, %sqsum, %3;
@%cont	bra		skip;

	st.global.u32	[%ar_result], %niter;
)";
	result += loop_end.arg (z2r_high, z2i_high, dem ? "100" : "10000");

	gen_store_zprev (result, zr, zi, n_prev);
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
	result += cg.code ();
	result += "}\n";
}

char *gen_mprec_funcs (formula f, int size, int stepsize, int power, int n_prev)
{
	QString result;
	result += QString (R"(	.version	6.2
	.target	sm_61
	.address_size 64
	// .extern .shared .align 4 .b8 shared[];
	.const .align 4 .u32 const_origin_x[%1];
	.const .align 4 .u32 const_origin_y[%1];
	.const .align 4 .u32 const_param_p[%1];
	.const .align 4 .u32 const_param_q[%1];
	.const .align 4 .u32 const_matrix00[%1];
	.const .align 4 .u32 const_matrix01[%1];
	.const .align 4 .u32 const_matrix10[%1];
	.const .align 4 .u32 const_matrix11[%1];
	.const .align 4 .u32 const_critpoint[%1];
)").arg (stepsize * 2);

#if 0
	gen_mul_func (result, size, false);
	gen_mul_func (result, size, true);
#endif
	gen_coord_muladd (result, size, stepsize);
	gen_coord_mul (result, size, stepsize);

	gen_kernel (f, result, size, stepsize, power, n_prev, true, false);
	gen_kernel (f, result, size, stepsize, power, n_prev, false, false);
	if (formula_supports_dem (f)) {
		gen_kernel (f, result, size, stepsize, power, n_prev, true, true);
		gen_kernel (f, result, size, stepsize, power, n_prev, false, true);
	} else if (formula_supports_hybrid (f)) {
		gen_kernel (f, result, size, stepsize, power, n_prev, true, false, true);
		gen_kernel (f, result, size, stepsize, power, n_prev, false, false, true);
	}

	// std::cerr << result;
	QByteArray a = result.toLatin1 ();
#if 0
	std::cerr << a.constData ();
#endif
	return strdup (a.constData ());
}
