#ifndef RENDER_PARAMS_H
#define RENDER_PARAMS_H

struct render_params
{
	QVector<uint32_t> palette;
	enum class smooth_t { std, makin };
	smooth_t smooth;
	uint32_t incol;
	bool oc_atom, ic_atom;
	int mod_type;
	double steps;
	int col_off, basin_col_off;
	int sub_val;
	int angle;
	bool sub;
	bool sac, tia, sac_contrast, sac_fade;
	double sac_factor, tia_power;
	int sac_fade_amount;
	bool dem_colour, angle_colour;
	bool dem;
	bool dem_shade;
	uint32_t dem_start, dem_stop;
	uint32_t bin_a, bin_b;
	uint32_t sac_tint;
	double dem_param, dem_strength;
	// Used for stored parameters, holds either the aspect set in the GUI,
	// or, if that is disabled, the image dimensions.
	double aspect;
	// Minimum niter value, stored for batch rendering (which is done
	// piecewise and therefor cannot compute it)
	double minimum;
};

#endif
