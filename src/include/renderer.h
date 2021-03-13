#ifndef RENDERER_H
#define RENDERER_H

#include <QMutex>
#include <QImage>
#include <QThreadPool>

#include "fpvec.h"
#include "fractal.h"
#include "render-params.h"

/* We take a color table and interpolate this many times.  When rendering the image, a simple linear
   interpolation is used for whatever small differences are left.  */
constexpr int interpolation_factor = 256;

class QGraphicsView;

class Renderer : public QObject
{
	Q_OBJECT
	QThreadPool m_pool;
	double m_minimum = 0;
	int m_min_gen = -1;
	double m_tia_power = 0;
	double m_sac_density = 0;
public:

	// One big mutex around the drawing function
	QMutex mutex;

	// Communication with Mainwindow
	QMutex queue_mutex;
	bool queued = false;
	render_params next_rp;
	int render_width, render_height;
	QImage result_image;

	void do_render (const render_params &rp, int w, int h, int yoff, frac_desc *, QGraphicsView *, int);
	void slot_render (frac_desc *, QGraphicsView *, int);
	void set_minimum (double m, int gen) { m_minimum = m; m_min_gen = gen; }
signals:
	void signal_render_complete (QGraphicsView *, frac_desc *, QImage, double);
};

extern const QVector<uint32_t> &palette_from_index (int cb_idx, const QVector<uint32_t> &custom);

#endif
