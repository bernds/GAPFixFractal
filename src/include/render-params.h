#ifndef RENDER_PARAMS_H
#define RENDER_PARAMS_H

struct render_params
{
	QVector<uint32_t> palette;
	uint32_t incol;
	int mod_type;
	double steps;
	int slider;
	int sub_val;
	int angle;
	bool sub;
	bool sac, tia, sac_contrast;
	double sac_factor, tia_power;
	bool dem_colour, angle_colour;
	bool dem;
	bool dem_shade;
	uint32_t dem_start, dem_stop;
	uint32_t bin_a, bin_b;
	double dem_param, dem_strength;
	// Used for stored parameters, holds either the aspect set in the GUI,
	// or, if that is disabled, the image dimensions.
	double aspect;
	// Minimum niter value, stored for batch rendering (which is done
	// piecewise and therefor cannot compute it)
	double minimum;
};

#endif
