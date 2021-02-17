#ifndef FRACTAL_H
#define FRACTAL_H

#include <cuda.h>
#include "bitarray.h"
#include "fpvec.h"
#include "formulas.h"

// Describe the location of a fractal image.
struct frac_params
{
	bool julia = false;

	formula fm;

	vpvec center_x, center_y;
	vpvec width;
	vpvec matrix[2][2];

	vpvec param_p, param_q, critpoint;

	int nwords = 0;
	int power = 2;

	uint32_t maxiter = 100;
	// Adjusted during computation
	uint32_t maxiter_found = 0;
	// Adjusted during rendering
	double min_stripeval = 0;
	double max_stripeval = 1;

	// Used for the initial position
	int bounds_w = 0;
	int bounds_h = 0;

	uint32_t hybrid_code = 0;
	int hybrid_len = 0;

	int n_prev = 1;

	void resize (int max_nwords)
	{
		center_x.resize (max_nwords);
		center_y.resize (max_nwords);

		matrix[0][0].resize (max_nwords);
		matrix[0][1].resize (max_nwords);
		matrix[1][0].resize (max_nwords);
		matrix[1][1].resize (max_nwords);
		width.resize (max_nwords);
		param_p.resize (max_nwords * 2);
		param_q.resize (max_nwords * 2);
		critpoint.resize (max_nwords * 2);
	}
};

// Describe other parameters necessary to compute the image.
struct frac_desc : public frac_params
{
	int n_pixels = 0;
	int n_threads = 0;
	int pixel_width, pixel_height;
	int pixel_step;
	int samples = 1;
	bool dem = false;
	int nvals_allocated = 0;

	bit_array pixels_done, pixels_started, pic_pixels_done;
	int n_completed = 0;
	uint32_t start_idx;

	vpvec step;
	/* Updated by the GUI and used to compute the actual rotation matrix.  */
	double rotation_angle = 0;

	double *pic_t;
	double *pic_zprev;
	double *pic_zder;
	uint32_t *pic_result;

	double *pic_iter_value;

	uint32_t *host_cplxvals, *host_coords, *host_result, *host_zpidx;
	double *host_zprev;
	CUdeviceptr cu_ar_origin = 0;
	CUdeviceptr cu_ar_cplxvals = 0;
	CUdeviceptr cu_ar_zprev = 0;
	CUdeviceptr cu_ar_zpidx = 0;
	CUdeviceptr cu_ar_result = 0;
	CUdeviceptr cu_ar_coords = 0;
	CUdeviceptr cu_ar_step = 0;

	double cmin, cmax;

	int generation;

	/* For rendering in several passes during batch export.  */
	int xoff = 0;
	int yoff = 0;
	int full_height;

	void set (const frac_params &newfp)
	{
		*(frac_params *)this = newfp;
	}

private:
	frac_desc &operator=(const frac_desc &) = default;
};
#endif
