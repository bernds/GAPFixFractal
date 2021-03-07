#include <cmath>
#include <cfloat>
#include <QVector>
#include <QSemaphore>

#include "renderer.h"

static inline uint32_t color_merge (uint32_t c1, uint32_t c2, double m1)
{
	double m2 = 1 - m1;
	int r = ((c1 >> 16) & 255) * m1 + ((c2 >> 16) & 255) * m2;
	int g = ((c1 >> 8) & 255) * m1 + ((c2 >> 8) & 255) * m2;
	int b = ((c1 >> 0) & 255) * m1 + ((c2 >> 0) & 255) * m2;
	return (r << 16) | (g << 8) | b;
}

static QVector<uint32_t> cols_bw = { 0, 0xFFFFFF };
static QVector<uint32_t> cols_bcw = {
	0x000000, 0x801100, 0xffffff, 0x000000,
	0xcd845a, 0xffffff, 0x000000, 0x2020e0,
	0xffffff, 0x000000, 0x00bbcc, 0xffffff,
	0x000000, 0x004488, 0xffffff, 0x000000,
	0x6644cc, 0xffffff, 0x000000, 0x5555ee,
	0xffffff, 0x000000, 0x954928, 0xffffff,
	0x000000, 0xeb8322, 0xffffff, 0x000000,
	0x302010, 0xffffff, 0x000000, 0x007000,
	0xffffff, 0x000000, 0x88aa66, 0xffffff,
	0x000000, 0x447711, 0xffffff, 0x000000,
	0xcc3311, 0xffffff, 0x000000, 0x772277,
	0xffffff, 0x000000, 0x4477aa, 0xffffff
};

static QVector<uint32_t> cols_bcwc = {
	0x000000, 0x112550, 0xffffff, 0x8af4fe,
	0x000000, 0xc7e9fe, 0xffffff, 0xc4ebeb,
	0x000000, 0x55cccc, 0xffffff, 0x41758e,
	0x000000, 0x2b7d8e, 0xffffff, 0x32a8cc,
	0x000000, 0x79dceb, 0xffffff, 0x36b9c2,
	0x000000, 0x172e4a, 0xffffff, 0x15202d,
	0x000000, 0x802908, 0xffffff, 0xbb4411,
	0x000000, 0xe4c1bc, 0xffffff, 0xdec4be,
	0x000000, 0xaa5533, 0xffffff, 0x72341e,
	0x000000, 0x87160d, 0xffffff, 0xebc8c1,
	0x000000, 0xc8596c, 0xffffff, 0xdc7d98,
	0x000000, 0xd27c95, 0xffffff, 0xe0bbca,
	0x000000, 0xb0006e, 0xffffff, 0x770046,
	0x000000, 0x7213a6, 0xffffff, 0xaa22e8,
	0x000000, 0xdebdf6, 0xffffff, 0xdbe1f6,
	0x000000, 0xa2b3e8, 0xffffff, 0x6c7aa6,
	0x000000, 0x592940, 0xffffff, 0x884466,
	0x000000, 0xd2c1c8, 0xffffff, 0xdec1d2,
	0x000000, 0xaa4488, 0xffffff, 0x722959,
	0x000000, 0x131e65, 0xffffff, 0x223399,
	0x000000, 0xbdbed7, 0xffffff, 0xbec4d7,
	0x000000, 0x335599, 0xffffff, 0x1e3465,
	0x000000, 0x5929ab, 0xffffff, 0x8844ee,
	0x000000, 0xd2c1f8, 0xffffff, 0xcdd2d7,
	0x000000, 0x778899, 0xffffff, 0x4c5965,
	0x000000, 0x801e8e, 0xffffff, 0xbb33cc,
	0x000000, 0xe4beeb, 0xffffff, 0xc4ebe4,
	0x000000, 0x55ccbb, 0xffffff, 0x348e80,
	0x000000
};

static QVector<uint32_t> cols_a = {
	0x080202, 0x77321c, 0xe68300, 0xfdbd63,
	0x080202, 0xcc4422, 0x53120c, 0x000200,
	0x3b5829, 0x0d2312, 0x000000, 0x046008,
	0xd6ff6d, 0xffffff, 0x000000, 0x77bbee,
	0x4477aa, 0x113388, 0x000020, 0xffed93,
	0x300000
};
static QVector<uint32_t> cols_b = {
	0x000000, 0xb41500, 0x200404, 0xee8833, 0xffee00, 0xffffff, 0xb1c9fd, 0x130542,
	0x2233aa, 0x92c7ff, 0x200404, 0xf73e06, 0x0f0715, 0x331166, 0x7739a6, 0xffffff,
	0xeecc33, 0xe68300, 0x200404, 0xcc3322, 0x702c0a, 0xffeedd, 0x003800, 0x96c963,
	0x225511, 0xcceebb, 0x047010, 0x000000, 0xffffbb, 0xeeeeff, 0x77bbee, 0x4477aa,
	0x113388
};
static QVector<uint32_t> cols_c = {
	0, 0x0000C0, 0x66b9CC, 0xFFFFFF, 0, 0x004488, 0x5c91bf, 0xFFFFFF,
	0x502210, 0xd75c14, 0, 0xFF8040, 0, 0xFFD700, 0xFFDD99, 0, 0xCCCCEE
};
static QVector<uint32_t> cols_d = {
	0xce371c, 0xee8833, 0xf6b87b, 0x000000, 0xfa6a3e, 0xce250b, 0x301d1a, 0xefefef,
	0xfbb6c2, 0x9a4b5c, 0x000000, 0x6a99ff, 0x0e233a, 0xffffff, 0xbdc6f0, 0x1c1d6a,
	0x69afca, 0x252e4a, 0x0c1cd3, 0xffe0b9, 0x361f15, 0x1f1a7a, 0xa1d2fc, 0x696bee,
	0x1c1270, 0x000000, 0xffffbb, 0xe3e3f3, 0x77bbee, 0x040842, 0xa9b3fd, 0x000000
};
static QVector<uint32_t> cols_e = {
	0x5e69b5, 0xffffff, 0xeecf04, 0xe68300,
	0x200404, 0x7a1d14, 0xffffff, 0x2930b2,
	0x200404, 0xfe962e, 0xffffbb, 0xffffff,
	0xa5b4ee, 0x1c1e8d, 0x000000, 0x6b120a,
	0x000020, 0xffffff,
};

static QVector<uint32_t> cols_f = {
	0x00001e, 0x62a233, 0x88dd40, 0x62a22c,
	0x000000, 0x2f6f23, 0x79da40, 0x2f6f23,
	0xfefdfe, 0x13260a, 0x496522, 0xaaee55,
	0x7bae3c, 0x000000, 0xa33d0d, 0xd75c14,
	0x000000, 0xfea050, 0xbb7438, 0x000000,
	0xfec678, 0x000000, 0xfedd99, 0xfefdfe,
	0x000000, 0xfed9c3, 0xfea050, 0xbb6b31,
	0x000000, 0x9d410c, 0xfefdfe, 0xd75c14,
	0x9d410c
};

static QVector<uint32_t> cols_g = {
	0x000000, 0x200a05, 0x3b160c, 0x4d1e12,
	0x9a8e8c, 0xe2dfdf, 0xffffff, 0x000000,
	0x361309, 0x5f2614, 0x77321c, 0xad918d,
	0xe7e0df, 0xffffff, 0x000000, 0x5c2a05,
	0x994b0c, 0xb96012, 0xd3a18c, 0xf1e4df,
	0xffffff, 0x000000, 0x793c00, 0xc56900,
	0xe68300, 0xefb38b, 0xfae8df, 0xffffff,
	0x000000, 0x83501d, 0xd38836, 0xf4a646,
	0xf8c797, 0xfceee2, 0xffffff, 0x000000,
	0x895e2b, 0xdc9d4e, 0xfbbd63, 0xfdd6a3,
	0xfef2e4, 0xffffff, 0x000000, 0x895e2b,
	0xdc9d4e, 0xfbbd63, 0xfdd6a3, 0xfef2e4,
	0xffffff, 0x000000, 0x83501d, 0xd38836,
	0xf4a646, 0xf8c797, 0xfceee2, 0xffffff,
	0x000000, 0x793c00, 0xc56900, 0xe68300,
	0xefb38b, 0xfae8df, 0xffffff, 0x000000,
	0x5c2a05, 0x994b0c, 0xb96012, 0xd3a18c,
	0xf1e4df, 0xffffff, 0x000000, 0x361309,
	0x5f2614, 0x77321c, 0xad918d, 0xe7e0df,
	0xffffff, 0x000000, 0x461107, 0x792310,
	0x952f17, 0xbd918d, 0xebe0df, 0xffffff,
	0x000000, 0x681c0b, 0xab3419, 0xcc4422,
	0xdf978e, 0xf5e1e0, 0xffffff, 0x000000,
	0x681c0b, 0xab3419, 0xcc4422, 0xdf978e,
	0xf5e1e0, 0xffffff, 0x000000, 0x4d1308,
	0x822612, 0xa03219, 0xc4918d, 0xede0df,
	0xffffff, 0x000000, 0x461107, 0x792310,
	0x952f17, 0xbd918d, 0xebe0df, 0xffffff,
	0x000000, 0x230502, 0x400c07, 0x53120b,
	0x9c8c8b, 0xe3dfdf, 0xffffff, 0x000000,
	0x230502, 0x400c07, 0x53120b, 0x9c8c8b,
	0xe3dfdf, 0xffffff, 0x000000, 0x4d1308,
	0x822612, 0xa03219, 0xc4918d, 0xede0df,
	0xffffff
};

const QVector<uint32_t> &palette_from_index (int cb_idx, const QVector<uint32_t> &custom)
{
	return (cb_idx == 0 ? cols_bw
		: cb_idx == 1 ? cols_bcw
		: cb_idx == 2 ? cols_bcwc
		: cb_idx == 3 ? cols_a
		: cb_idx == 4 ? cols_b
		: cb_idx == 5 ? cols_c
		: cb_idx == 6 ? cols_d
		: cb_idx == 7 ? cols_e
		: cb_idx == 8 ? cols_f
		: cb_idx == 9 ? cols_g
		: custom);
}

static int power_from_fp (frac_params &fp)
{
	int power = fp.power;
	if (fp.fm == formula::sqtwice_a || fp.fm == formula::sqtwice_b)
		power += 2;
	return power;
}

#if 0
static int col_sum (uint32_t col)
{
	return (col & 255) + ((col >> 8) & 255) + (col >> 16);
}
#endif

static inline double niter_transfer (double niter, int type)
{
	switch (type) {
	default:
		return niter;
	case 1:
		return sqrt (niter);
	case 2:
		return cbrt (niter);
	case 3:
		return pow (niter, 1/4.);
	case 4:
		return log (niter);
	case 5:
		return log (sqrt (niter));
	case 6:
		return log (log (niter));
	}
}

static inline QRgb color_from_niter (const QVector<uint32_t> &palette, double niter, double sub_val, int type, double steps,
				     int slider)
{
	niter = niter_transfer (niter, type) - sub_val;
	niter *= interpolation_factor;
	int x = niter / steps;
	double y = niter - x * steps;
	x += 0.5 * slider * interpolation_factor;
	size_t size = palette.size ();
	// return palette[x % size];
	double m1 = y / steps;
	uint32_t col1 = palette[(x + 1) % size];
	uint32_t col2 = palette[x % size];
	uint32_t primary = color_merge (col1, col2, m1);
	return primary | 0xFF000000;
}

static inline QRgb color_from_period (uint32_t p)
{
	uint32_t str = ((p & 2) << 5) | (p << 7);
	str &= 0xC0;
	str = 0xFF - str;
	uint32_t col = ((p & 1) | ((p & 2) << 7) | ((p & 4) << 14)) * str;
	p >>= 3;
	uint32_t col2 = ((p & 1) | ((p & 2) << 7) | ((p & 4) << 14)) * 0x80;
	col ^= col2;
	return col;
}

static inline uint32_t modify_color (uint32_t col, double v)
{
	QColor c = QColor::fromRgb (col);

	int cr = c.red ();
	int cg = c.green ();
	int cb = c.blue ();

	// Lighten dark colors, darken light colors.  Try to distinguish between
	// the two with something better than just the HSV value.
	int sum = (cr * 13762 + cg * 47186 + cb * 4588) / 65536;
	v -= (1 - sum / 255.) / 2;
	if (v < 0) {
		cr = 255 - (255 - cr) * (v + 1);
		cg = 255 - (255 - cg) * (v + 1);
		cb = 255 - (255 - cb) * (v + 1);
	} else {
		cr *= 1 - v;
		cg *= 1 - v;
		cb *= 1 - v;
	}
	col = QColor::fromRgb (cr, cg, cb).rgb ();
	return col;
}

inline std::pair<double, int> iter_value_at (frac_desc *fd, int idx, int power, render_params::smooth_t smooth)
{
	double v = fd->pic_result[idx];
	if (v == 0)
		return { 0, 0 };
	double re = fd->pic_zprev[idx * 2 * fd->n_prev];
	double im = fd->pic_zprev[idx * 2 * fd->n_prev + 1];
	double re2 = re * re;
	double im2 = im * im;
	double radius = fd->dem ? 10 : 100;
	int attractor = re2 + im2 < 16 ? 1 : 0;
	if (smooth == render_params::smooth_t::makin && attractor != 0) {
		// Formula from fractalforums, by David Makin
		double log_tolerance = log (1.0 / (1 << 28)) / 2;
		double prev1re = fd->pic_zprev[idx * 2 * fd->n_prev + 2];
		double prev1im = fd->pic_zprev[idx * 2 * fd->n_prev + 3];
		double prev2re = fd->pic_zprev[idx * 2 * fd->n_prev + 4];
		double prev2im = fd->pic_zprev[idx * 2 * fd->n_prev + 5];
		double rdiff0 = prev1re - re;
		double idiff0 = prev1im - im;
		double diff0 = rdiff0 * rdiff0 + idiff0 * idiff0;
		double rdiff1 = prev2re - prev1re;
		double idiff1 = prev2im - prev1im;
		double diff1 = rdiff1 * rdiff1 + idiff1 * idiff1;
		double log0 = log (diff0);
		double log1 = log (diff1);
		v += (log (-log_tolerance) - log (-log0/2)) / log (log0 / log1);
		return { v, attractor };
	}
	// Some attractor other than infinity.
	if (attractor != 0)
		return { v, 0 };
	double correction = log (0.5 * log ((double)re2 + im2) / log (radius)) / log (power);
	fd->cmin = std::min (fd->cmin, correction);
	fd->cmax = std::max (fd->cmax, correction);
	return { v + 5 - correction, 0 };
}

static inline double compute_sac (double iter_val, double *prev_vals, int n_prev, int count, double density, double radius, double power, int fade_amount)
{
	double re = prev_vals[0];
	double im = prev_vals[1];
	double re2 = re * re;
	double im2 = im * im;
	double firstval = 0.5 * sin (density * atan2 (im, re));
	double lastval = 0;
	double sum = 0;
	for (int i = 1; i < count; i++) {
		double lastre = prev_vals[i * 2];
		double lastim = prev_vals[i * 2 + 1];
		double val = 0.5 * sin (density * atan2 (lastim, lastre));
		lastval = val;
		sum += val;
	}
	int avg2_count = count == n_prev ? count - 1 : count;
	if (count < n_prev)
		lastval = 0;
	double avg1 = count > 1 ? sum / (count - 1) : 0;
	double avg2 = avg2_count > 0 ? (sum + firstval - lastval) / avg2_count : 0;
	double mixfactor = log (0.5 * log (re2 + im2) / log (radius)) / log (power);
	double result1 = avg1 * mixfactor + avg2 * (1 - mixfactor);
	if (fade_amount != 0) {
		iter_val /= fade_amount;
		result1 *= iter_val / (1 + iter_val);
	}
	return result1 + 0.5;
}

static inline double compute_tia (double *prev_vals, double cre, double cim, int count,
				  double tia_power, double radius, double power)
{
	double re = prev_vals[0];
	double im = prev_vals[1];
	double re2 = re * re;
	double im2 = im * im;
	double mixfactor = log (0.5 * log (re2 + im2) / log (radius)) / log (power);
	double cmag = sqrt (cre * cre + cim * cim);
	if (cmag == 0)
		count = 0;
	double sum = 0;
	double firstval = 0;
	double curre = re;
	double curim = im;
	for (int i = 1; i < count; i++) {
		double lastre = prev_vals[i * 2];
		double lastim = prev_vals[i * 2 + 1];
		double mag = sqrt (re2 + im2);
		// We need |zprev^n|, which is (for normal formulas) equal to
		// |z-c|.
		double zlre = curre - cre;
		double zlim = curim - cim;
		double zlmag = sqrt (zlre * zlre + zlim * zlim);
		double lowbound = fabs (zlmag - cmag);
		double mod = mag - lowbound;
		mod /= zlmag + cmag - lowbound;
		mod = pow (mod, tia_power);
		if (i == 1)
			firstval = mod;
		else
			sum += mod;
		curre = lastre;
		curim = lastim;
		re2 = lastre * lastre;
		im2 = lastim * lastim;
	}
	double avg1 = count > 2 ? sum / (count - 2) : 0;
	double avg2 = count > 1 ? (sum + firstval) / (count - 1) : 0;
	return avg1 * mixfactor + avg2 * (1 - mixfactor);
}

static inline uint32_t color_for_sac_common (uint32_t col, double mod, const frac_desc *fd, const render_params &rp)
{
	mod = std::clamp (mod, 0.0, 1.0);

	if (rp.sac_contrast) {
		double prev_width = fd->max_stripeval - fd->min_stripeval;
		mod = std::clamp ((mod - fd->min_stripeval) / prev_width, 0.0, 1.0);
	}
	if (rp.angle_colour)
		col = modify_color (col, mod);
	else {
		QColor base = QColor::fromRgb (rp.sac_tint);
		int cr = base.red ();
		int cg = base.green ();
		int cb = base.blue ();
		cr *= mod;
		cg *= mod;
		cb *= mod;
		col = QColor::fromRgb (cr, cg, cb, floor (mod * 255)).rgba ();
	}
	return col;
}

class runner : public QRunnable
{
	QSemaphore *completion_sem;
	std::atomic<bool> *success, *abort_render;
	const render_params &rp;
	int w, y0, y0e;
	frac_desc *fd;
	double minimum;
	QRgb *data;
	double dstep;
	double colour_sub_val;
	int power;
	int fade_amount;
public:
	runner (QSemaphore *sem, std::atomic<bool> *succ_in, std::atomic<bool> *abrt, const render_params &rp_in,
		int w_in, int y0_in, int y0e_in, frac_desc *fd_in, double min_in, QRgb *data_in)
		: completion_sem (sem), success (succ_in), abort_render (abrt),
		  rp (rp_in), w (w_in), y0 (y0_in), y0e (y0e_in), fd (fd_in), minimum (min_in), data (data_in)
	{
		setAutoDelete (true);
		colour_sub_val = rp.sub ? niter_transfer (rp.sub_val, rp.mod_type) : 0;
		fade_amount = rp.sac_fade ? rp.sac_fade_amount : 0;
	}

	void compute_color (size_t idx, int &r, int &g, int &b, int &alpha, int &outcolor)
	{
		int n_prev = fd->n_prev;
		double re = fd->pic_zprev[idx * 2 * n_prev];
		double im = fd->pic_zprev[idx * 2 * n_prev + 1];

		auto [ v, attractor ] = iter_value_at (fd, idx, power, rp.smooth);
		double v1 = v;
		if (v != 0) {
			outcolor++;

			if (rp.sub)
				v -= minimum - rp.sub_val;
			uint32_t col;
			if (rp.oc_atom)
				col = color_from_period (fd->pic_period[idx]);
			else
				col = color_from_niter (rp.palette, v, colour_sub_val, rp.mod_type, rp.steps, rp.col_off + attractor * rp.basin_col_off);

			if (rp.sac && attractor == 0) {
				double radius = fd->dem ? 10 : 100;
				int thisnp = std::min ((uint32_t)n_prev, fd->pic_result[idx]);
				double mod = compute_sac (v - rp.sub_val, fd->pic_zprev + idx * 2 * n_prev, n_prev, thisnp, rp.sac_factor, radius, power, fade_amount);
				col = color_for_sac_common (col, mod, fd, rp);
			} else if (rp.tia && attractor == 0) {
				double radius = fd->dem ? 10 : 100;
				int thisnp = std::min ((uint32_t)n_prev, fd->pic_result[idx]);
				double cre = fd->pic_t[idx * 2];
				double cim = fd->pic_t[idx * 2 + 1];
				double mod = compute_tia (fd->pic_zprev + idx * 2 * n_prev, cre, cim,
							  thisnp, rp.tia_power, radius, power);
				col = color_for_sac_common (col, mod, fd, rp);
			}
			else if (rp.angle == 1) {
				double re2 = re * re;
				double im2 = im * im;
				double size = sqrt (re2 + im2);
				double init_angle = atan2 (im, re);
				init_angle += M_PI;
#if 0
				double compare_angle = M_PI * v1 * 2;
				compare_angle -= M_PI * 2 * floor (compare_angle / (M_PI * 2));
				init_angle = init_angle - compare_angle + M_PI * 2;
#endif
#if 1
				double ang = init_angle / 2;
				double v = fabs (cos (ang));
#else
				double ang = init_angle - M_PI * 2 * floor (init_angle / (M_PI * 2));
				ang += M_PI / 2;
				ang = fmod (ang, M_PI * 2);
				double v = fabs (ang / (M_PI * 2) - 0.5);
#endif
				col = modify_color (col, v);
			}
			else if (rp.angle == 2) {
				if (im < 0)
					col = rp.bin_a;
				else if (!rp.angle_colour)
					col = rp.bin_b;
			}
			double dem_shade = 1;
			if (rp.dem || rp.dem_shade) {
				double re2 = re * re;
				double im2 = im * im;
				double rezder = fd->pic_zder[idx * 2];
				double imzder = fd->pic_zder[idx * 2 + 1];
				double rezder2 = rezder * rezder;
				double imzder2 = imzder * imzder;
				if (rp.dem_shade) {
					double reu = re * rezder + im * imzder;
					double imu = im * rezder - re * imzder;
					double udiv = rezder2 + imzder2;
					reu /= udiv;
					imu /= udiv;
					double uabs = sqrt (reu * reu + imu * imu);
					if (uabs != 0) {
						reu /= uabs;
						imu /= uabs;
					}
					double light_height = 1.5;
					double light_re = -M_SQRT1_2;
					double light_im = -M_SQRT1_2;
					double t = light_re * reu + light_im * imu + light_height;
					t = t / (1 + light_height);
					if (t < 0)
						t = 0;
					dem_shade = 1 - (1 - t) * rp.dem_strength;
				}
				if (rp.dem) {
					double divisor = sqrt (rezder2 + imzder2);
					double dist = divisor == 0 ? rp.dem_param + 1 : log (re2 + im2) * sqrt (re2 + im2) / divisor;
					double dem_dist = 0;
					if (dist > rp.dem_param)
						dem_dist = 1;
					else {
						dem_dist = dist / rp.dem_param;
						// The result by itself seems too steep to be useful for rendering
						// shades of grey. Could make this configurable, but taking cbrt
						// gives images that I find pleasing.
						dem_dist = cbrt (dem_dist);
						dem_dist = 1 - (1 - dem_dist) * rp.dem_strength;
						uint32_t v = color_merge (rp.dem_start, rp.dem_stop, dem_dist);
						r += (v >> 16) & 255;
						g += (v >> 8) & 255;
						b += v & 255;
						return;
					}
				}
				if (!rp.dem_colour && !rp.sac && !rp.tia && !rp.angle) {
					alpha += 255;
					r += dem_shade * 255;
					g += dem_shade * 255;
					b += dem_shade * 255;
					return;
				}
			}
			double dem_factor = dem_shade;
			alpha += ((col >> 24) & 0xFF) * dem_factor;
			r += ((col >> 16) & 0xFF) * dem_factor;
			g += ((col >> 8) & 0xFF) * dem_factor;
			b += (col & 0xFF) * dem_factor;
		} else if (rp.ic_atom) {
			uint32_t p = fd->pic_period[idx];
			if (p != 0) {
				outcolor++;
				uint32_t col = color_from_period (p);
#if 0
				double val = fd->pic_doubles[idx];
				val = 1 / (1 + 1000 * val);
#else
				double val = 1;
#endif
				r += ((col >> 16) & 0xFF) * val;
				g += ((col >> 8) & 0xFF) * val;
				b += (col & 0xFF) * val;
			}
		}
	}

	void run () override
	{
		bool any_found = false;
		uint32_t in_color1 = rp.incol & 0xFF;
		int sample_steps = fd->samples;
		power = power_from_fp (*fd);
		int ss2 = sample_steps * sample_steps;
		dstep = to_double (fd->step);
		if (dstep == 0)
			dstep = DBL_EPSILON;
		for (int y = y0; y < y0e; y++) {
			if (abort_render->load ()) {
				break;
			}
			for (int x = 0; x < w; x++) {
				int r = 0, g = 0, b = 0, alpha = 0;
				int valid = 0;
				int outcolor = 0;
				unsigned int cx1 = x * sample_steps;
				unsigned int cy1 = y * sample_steps;
				unsigned int idx = cy1 * fd->pixel_width + cx1;
				if (!fd->pic_pixels_done.test_bit (idx) /* && fd->pic_result[idx] == 0 */) {
					unsigned int mask = 0;
					mask = ~mask;
					mask *= sample_steps;
					for (int i = 0; i < 1; i++) {
						mask <<= 1;
						cy1 &= mask;
						cx1 &= mask;
						if (fd->pic_pixels_done.test_bit (cy1 * fd->pixel_width + cx1)) {
							valid++;
							compute_color (cy1 * fd->pixel_width + cx1, r, g, b, alpha, outcolor);
							break;
						}
					}
				} else {
					for (int y1 = 0; y1 < sample_steps; y1++) {
						for (int x1 = 0; x1 < sample_steps; x1++) {
							unsigned int cx = cx1 + x1;
							unsigned int cy = cy1 + y1;
							if (!fd->pic_pixels_done.test_bit (cy * fd->pixel_width + cx))
								continue;
							valid++;
							compute_color (cy * fd->pixel_width + cx, r, g, b, alpha, outcolor);
						}
					}
				}
				if (outcolor > 0) {
					any_found = true;
					int div = outcolor;
					if (valid == ss2) {
						alpha += 0xFF * (ss2 - outcolor);
						r += in_color1 * (ss2 - outcolor);
						g += in_color1 * (ss2 - outcolor);
						b += in_color1 * (ss2 - outcolor);
						div = ss2;
					}
					*data++ = ((alpha / div) << 24) + ((r / div) << 16) + ((g / div) << 8) + (b / div);
				} else if (valid) {
					*data++ = rp.incol;
				} else
					// todo: use previous scaled data
					*data++ = 0;
			}
		}
		if (any_found)
			success->store (true);
		completion_sem->release ();
	}
};

/* The view is null if we are being called for a batch render.  In
   that case we use values computed earlier in live view rather than
   recomputing them, since we only ever get to see a slice of the
   entire render.  */
void Renderer::do_render (const render_params &rp, int w, int h, int yoff, frac_desc *fd, QGraphicsView *view, int gen)
{
	QMutexLocker lock (&mutex);

	if (fd->generation != gen)
		return;

	fd->pic_pixels_done = fd->pixels_done;

	int power = power_from_fp (*fd);

	// printf ("dstep %f\n", dstep);
	double minimum = m_minimum;
	if (rp.sub && gen != m_min_gen) {
		minimum = 0;
		for (int i = 0; i < fd->n_pixels; i++) {
			if (!fd->pic_pixels_done.test_bit (i))
				continue;
			auto [ v, attractor ] = iter_value_at (fd, i, power, rp.smooth);
			v = floor (v);
			if (v > 0 && (minimum == 0 || v < minimum))
				minimum = v;
		}
	}
	bool recompute_sactia = false;
	if (view != nullptr) {
		recompute_sactia = ((rp.sac && (m_sac_density != rp.sac_factor || m_sac_fade != rp.sac_fade || (rp.sac_fade && m_fade_amount != rp.sac_fade_amount)))
				    || (rp.tia && (m_tia_power != rp.tia_power)));
		m_tia_power = rp.tia ? rp.tia_power : 0;
		m_sac_density = rp.sac ? rp.sac_factor : 0;
		m_sac_fade = rp.sac_fade;
		m_fade_amount = rp.sac_fade ? rp.sac_fade_amount : 0;
	}
	if ((rp.sac || rp.tia) && (gen != m_min_gen || recompute_sactia)) {
		/* If we are in batch mode, this value should have been precomputed and the generation
		   set up appropriately.  */
		assert (view != nullptr);
		/* ??? This is still suboptimal.  If we computed min/max first, we could have a smaller
		   number of bins - 256 should be enough.  But another loop through all the pixels would
		   be expensive.  */
		constexpr int n_bins = 50000;
		std::array<int, n_bins> bins {};
		double radius = fd->dem ? 10 : 100;
		int count = 0;
		// fd->min_stripeval = 1;
		// fd->max_stripeval = 0;

		for (int i = 0; i < fd->n_pixels; i += 13) {
			if (!fd->pic_pixels_done.test_bit (i))
				continue;
			int r = fd->pic_result[i];
			if (r == 0)
				continue;
			r = std::min (fd->n_prev, r);
			count++;
			double v;
			if (rp.sac) {
				auto [ iterval, attractor ] = iter_value_at (fd, i, power, rp.smooth);
				v = compute_sac (iterval - minimum, fd->pic_zprev + i * 2 * fd->n_prev, fd->n_prev, r, rp.sac_factor, radius, power, m_fade_amount);
			} else {
				double cre = fd->pic_t[i * 2];
				double cim = fd->pic_t[i * 2 + 1];
				v = compute_tia (fd->pic_zprev + i * 2 * fd->n_prev, cre, cim,
						 r, rp.tia_power, radius, power);
			}
			v = std::clamp (v, 0.0, 1.0);
			// fd->min_stripeval = std::min (v, fd->min_stripeval);
			// fd->max_stripeval = std::max (v, fd->max_stripeval);
			bins[trunc (v * (n_bins - 1))]++;
		}
		if (1) {
			double margin1 = count / 400;
			double margin2 = margin1 * 399;
			int total = 0;
			for (int i = 0; i < n_bins; i++) {
				int newt = total + bins[i];
				if (total < margin1 && newt >= margin1)
					fd->min_stripeval = (double)i / (n_bins - 1);
				if (total < margin2 && newt >= margin2)
					fd->max_stripeval = (double)(i + 1) / (n_bins - 1);
				total = newt;
			}
		}
	}
	if (fd->n_completed == fd->n_pixels) {
		m_minimum = minimum;
		m_min_gen = gen;
	}

	QSemaphore completion_sem (0);

	QRgb *data = (QRgb *)result_image.bits ();
	data += w * yoff;
	int tc = std::max (1, m_pool.maxThreadCount ());
	int lines_per_thread = std::max (20, (h + tc - 1) / tc);
	int n_started = 0;
	std::atomic<bool> any_found = false;
	for (int y0 = 0; y0 < h; y0 += lines_per_thread) {
		int y0e = y0 + lines_per_thread;
		if (y0e > h)
			y0e = h;
		m_pool.start (new runner (&completion_sem, &any_found, &abort_render,
					  rp, w, y0, y0e, fd, minimum, data + w * y0));
#ifdef DEBUG_SINGLETHREAD
		completion_sem.acquire (1);
#else
		n_started++;
#endif
	}
	completion_sem.acquire (n_started);
	if (!any_found || abort_render)
		return;

	if (view != nullptr)
		emit signal_render_complete (view, fd, result_image, minimum);
}

void Renderer::slot_render (frac_desc *fd, QGraphicsView *view, int generation)
{
	queue_mutex.lock ();
	assert (queued);
	queued = false;
	abort_render = false;
	render_params this_p = next_rp;
	// Shouldn't happen, just making sure...
	if (fd->pic_t == nullptr)
		this_p.tia = false;
	if (fd->n_prev == 1)
		this_p.sac = this_p.tia = false;
	if (fd->n_prev < 4)
		this_p.smooth = render_params::smooth_t::std;
	if (!fd->incolor)
		this_p.ic_atom = this_p.oc_atom = false;

	int w = render_width;
	int h = render_height;
	queue_mutex.unlock ();
	result_image = QImage (w, h, QImage::Format_RGB32);
	do_render (this_p, w, h, 0, fd, view, generation);
}

