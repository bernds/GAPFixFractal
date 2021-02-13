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

	// Used for the initial position
	int bounds_w = 0;
	int bounds_h = 0;

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

	bit_array pixels_done, pixels_started, pic_pixels_done;
	int n_completed = 0;
	uint32_t start_idx;

	vpvec step;
	/* Updated by the GUI and used to compute the actual rotation matrix.  */
	int rotation_angle = 0;

	int32_t *pic_z;
	int32_t *pic_zder;
	uint32_t *pic_z2;
	uint32_t *pic_result;

	double *pic_iter_value;

	uint32_t *host_t, *host_z, *host_zder, *host_z2, *host_coords, *host_result;
	CUdeviceptr cu_ar_origin = 0;
	CUdeviceptr cu_ar_z = 0;
	CUdeviceptr cu_ar_z2 = 0;
	CUdeviceptr cu_ar_zder = 0;
	CUdeviceptr cu_ar_tmp = 0;
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
