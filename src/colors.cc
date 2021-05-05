#include <cmath>
#include <QColor>

#include "colors.h"

// This is an experiment to apply the traditional colour star to hue rotations.
// This one has purple vs yellow, red vs green as opposite colour pairs.

// Color	HSV	trad
// Red		0	0
// Orange	30	60
// Yellow	60	120
// Yellow green	75	150
// Green	120	180
// Blue green	195	210
// Blue		240	240
// Purple	285	300

static double h2t_coeffs[8][4] = {
{-0.00027368755185889614, 0.01807898495841258, 1.703949247920629, 0, }, // 0:30
{+0.0001631720954003085, -0.021238383294915838, 2.8834702955204814, -11.795210475998525, }, // 30:60
{-0.0009738303938267054, 0.18342206476594666, -9.396156588131268, 233.79732719703648, }, // 60:75
{+0.0002908351895108052, -0.10112769148499322, 11.945075130689222, -299.7334657734758, }, // 75:120
{-0.000018930242946646284, 0.010387864199689322, -1.4367915514726821, 235.5412015130004, }, // 120:195
{+0.0001197745190042626, -0.07075442154159238, 14.38595416807725, -792.9372702577451, }, // 195:240
{-0.00023915102473215028, 0.1876719699486249, -47.636379789574896, 4168.849446354427, }, // 240:285
{+0.000177536815158701, -0.16859613315805297, 53.90002959582829, -5477.109445258877, } // 285:360
};

double hue_to_trad (double h)
{
	if (h == -1)
		return h;
	int cf = h <= 60 ? (int)h / 30 : h <= 75 ? 2 : h <= 120 ? 3 : h <= 195 ? 4 : h <= 240 ? 5 : h <= 285 ? 6 : 7;
	return h * h * h * h2t_coeffs[cf][0] + h * h * h2t_coeffs[cf][1] + h * h2t_coeffs[cf][2] + h2t_coeffs[cf][3];
}
static double t2h_coeffs[8][4] = {
{+0.00007766128116582096, -0.010295806121969555, 0.8381677551212178, 0, }, // 0:60
{-0.00004511286843011965, 0.011803540805299755, -0.4877930605149408, 26.519216312723174, }, // 60:120
{+0.00026325770923385326, -0.09920986715373051, 12.83381589456869, -506.3451418906221, }, // 120:150
{+0.00009062868691139928, -0.021526807108626197, 1.1813568878030445, 76.27780844766022, }, // 150:180
{-0.0006257724568794504, 0.3653298105384326, -68.45283428866755, 4254.3292790358955, }, // 180:210
{+0.00019023891838418008, -0.14875735587765457, 39.50547065871077, -3302.7520672805863, }, // 210:240
{+0.0001338597622888346, -0.10816436348900582, 29.763152485435068, -2523.3666134185305, }, // 240:300
{-0.00013756543452045356, 0.1361183136393535, -43.52165065307273, 4805.11370043225, } // 300:360
};
double trad_to_hue (double h)
{
	if (h == -1)
		return h;
#if 0
	if (h >= 360)
		h -= 360;
	int cf = h <= 120 ? (int)h / 60 : h <= 150 ? 2 : h <= 240 ? ((int)h - 150) / 30 + 3 : h <= 300 ? 6 : 7;
	return h * h * h * t2h_coeffs[cf][0] + h * h * t2h_coeffs[cf][1] + h * t2h_coeffs[cf][2] + t2h_coeffs[cf][3];
#else
	double left = 0;
	double right = 360;
	for (int i = 0; i < 9; i++) {
		double mid = left + (right - left) / 2;
		if (hue_to_trad (mid) < h)
			left = mid;
		else
			right = mid;
	}
	return left;
#endif
}

class trad_tester
{
public:
	trad_tester ()
	{
		for (int i = 0; i < 360; i++) {
			double vnew1 = hue_to_trad (i);
			double vnew = trad_to_hue (vnew1);
			double diff = std::abs (vnew - i);
			if (diff > 1)
				printf ("hue %d: to %d to %d, diff %f\n", i, (int)vnew1, (int)vnew, diff);
		}
	}
};
static trad_tester test;

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

static bool is_white (uint32_t col)
{
	return (col & 0xFFFFFF) == 0xFFFFFF;
}

static bool is_black (uint32_t col)
{
	return (col & 0xFFFFFF) == 0;
}

QVector<uint32_t> interpolate_colors (const QVector<uint32_t> &src, int steps, int hue_shift, bool trad,
				      bool narrow_blacks, bool narrow_whites, int nfactor, bool circular)
{
	QVector<double> rvals;
	QVector<double> gvals;
	QVector<double> bvals;
	int count = src.size ();
	for (int i = 0; i < count; i++) {
		QColor c = QColor::fromRgb (src[i]);
		if (hue_shift != 0) {
			int h = c.hslHue ();
			if (h != -1) {
				int s = c.hslSaturation ();
				int v = c.lightness ();
				double ht = (trad ? hue_to_trad (h) : h) + hue_shift;
				while (ht > 360)
					ht -= 360;
				h = trad ? trad_to_hue (ht) : ht;
				c = QColor::fromHsl (h, s, v);
			}
		}
		rvals.push_back (srgb_to_linear (c.red ()));
		gvals.push_back (srgb_to_linear (c.green ()));
		bvals.push_back (srgb_to_linear (c.blue ()));
	}
	QVector<uint32_t> new_colors;
	int i_count = circular ? count : count - 1;
	for (int i = 0; i < i_count; i++) {
		uint32_t thisc = src[i];
		uint32_t nextc = src[(i + 1) % count];

		int inc = 1;
		if (((is_black (thisc) || is_black (nextc)) && narrow_blacks)
		    || ((is_white (thisc) || is_white (nextc)) && narrow_whites))
			inc = nfactor;
		for (int step = 0; step < steps; step += inc) {
			double nr = linear_to_srgb (cubic_vector_offs (rvals, i, (double)step / steps));
			double ng = linear_to_srgb (cubic_vector_offs (gvals, i, (double)step / steps));
			double nb = linear_to_srgb (cubic_vector_offs (bvals, i, (double)step / steps));
			new_colors.push_back (QColor::fromRgb (nr, ng, nb).rgb ());
		}
	}
	if (!circular)
		new_colors.push_back (src.back ());
	return new_colors;
}
