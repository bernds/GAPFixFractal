#include <cmath>
#include <QColor>

#include "colors.h"

// Found on stackoverflow as "Fritsch-Carlson method", intended to ensure
// monotonicity (i.e. not overshooting the valid range in the case of
// interpolating colours.
static double choose_delta (double x0, double x1, double x2)
{
	double d0 = x1 - x0;
	double d1 = x2 - x1;
	double m = 0;
	if (d0 * d1 > 0) {
		m = 3 / (3/d0 + 3/d1);
	}
	return m;
}

static double cubic (double x0, double x1, double x2, double x3, double step)
{
	double m1 = choose_delta (x0, x1, x2);
	double m2 = choose_delta (x1, x2, x3);
	auto p = [=] (double t) -> double
	{
		double t2 = t*t;
		double t3 = t2*t;
		return (2*t3 - 3*t2 + 1) * x1 + (t3 - 2*t2 + t) * m1 + (-2*t3 + 3*t2) * x2 + (t3 - t2) * m2;
	};
	return p (step);
}

static double cubic_vector_offs (QVector<double> &v, int offs, double step)
{
	size_t sz = v.size ();
	double x0 = v[(offs + sz - 1) % sz];
	double x1 = v[offs % sz];
	double x2 = v[(offs + 1) % sz];
	double x3 = v[(offs + 2) % sz];
	return cubic (x0, x1, x2, x3, step);
}

static double srgb_to_linear (int v)
{
	double s = v / 255.;
	constexpr double a = 0.055;
	return pow ((s + a) / (1 + a), 2.4);
}

static int linear_to_srgb (double v)
{
	constexpr double a = 0.055;
	double nv = (1 + a) * pow (v, 1 / 2.4) - a;
	return nv * 255;
}

QVector<uint32_t> interpolate_colors (const QVector<uint32_t> &src, int steps,
				      bool narrow_blacks, bool narrow_whites, int nfactor)
{
	QVector<double> rvals;
	QVector<double> gvals;
	QVector<double> bvals;
	int count = src.size ();
	for (int i = 0; i < count; i++) {
		QColor c = QColor::fromRgb (src[i]);
		rvals.push_back (srgb_to_linear (c.red ()));
		gvals.push_back (srgb_to_linear (c.green ()));
		bvals.push_back (srgb_to_linear (c.blue ()));
	}
	QVector<uint32_t> new_colors;
	for (int i = 0; i < count; i++) {
		new_colors.push_back (src[i] | 0xFF000000);
		uint32_t thisc = src[i];
		uint32_t nextc = src[(i + 1) % count];

		int inc = 1;
		if (((thisc == 0 || nextc == 0) && narrow_blacks)
		    || ((thisc == 0xFFFFFF || nextc == 0xFFFFFF) && narrow_whites))
			inc = nfactor;
		for (int step = inc; step < steps; step += inc) {
			double nr = linear_to_srgb (cubic_vector_offs (rvals, i, (double)step / steps));
			double ng = linear_to_srgb (cubic_vector_offs (gvals, i, (double)step / steps));
			double nb = linear_to_srgb (cubic_vector_offs (bvals, i, (double)step / steps));
			new_colors.push_back (QColor::fromRgb (nr, ng, nb).rgb ());
		}
	}
	return new_colors;
}
