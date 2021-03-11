#include <cstdint>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cfloat>

#include <utility>

#include <QColorDialog>
#include <QProgressDialog>
#include <QVector>
#include <QFileDialog>
#include <QSettings>
#include <QThread>
#include <QCommandLineParser>
#include <QMessageBox>
#include <QTimer>
#include <QMouseEvent>
#include <QDataStream>
#include <QRunnable>
#include <QThreadPool>
#include <QElapsedTimer>

#include "gpuhandler.h"
#include "formulas.h"
#include "hybriddialog.h"
#include "gradeditor.h"
#include "colors.h"

#include "rotationdialog.h"
#include "batchrender.h"
#include "settings.h"
#include "mainwindow.h"
#include "ui_mainwindow.h"

#define PACKAGE "GAPFixFractal"

const formula formula_table[] = {
	formula::standard, formula::lambda, formula::spider, formula::tricorn,
	formula::ship, formula::mix, formula::sqtwice_a, formula::sqtwice_b,
	formula::celtic, formula::magnet_a, formula::facing, formula::facing_b,
	formula::testing
};

constexpr int default_power = 2;

/* We take a color table and interpolate this many times.  When rendering the image, a simple linear
   interpolation is used for whatever small differences are left.  */
constexpr int interpolation_factor = 256;

// Old versions of Qt which we want to support don't have the recommended
// range constructors, while the new version warns about the old style.
// Shut it up.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

QDataStream &operator<< (QDataStream &s, const frac_params &fp)
{
	s << (qint32)2;
	s << fp.julia;
	s << (qint32)fp.fm;
	s << QVector<uint32_t>::fromStdVector (fp.center_x);
	s << QVector<uint32_t>::fromStdVector (fp.center_y);
	s << QVector<uint32_t>::fromStdVector (fp.width);
	s << QVector<uint32_t>::fromStdVector (fp.param_p);
	s << QVector<uint32_t>::fromStdVector (fp.param_q);
	s << QVector<uint32_t>::fromStdVector (fp.critpoint);
	s << fp.nwords << fp.power << fp.maxiter;
	// We stream another formula::standard in case we ever want to support
	// things like hybrids of tricon/ship.
	s << fp.hybrid_code << fp.hybrid_len << (qint32)formula::standard;
	return s;
}

QDataStream &operator>> (QDataStream &s, frac_params &fp)
{
	qint32 version;
	qint32 fm;
	s >> version;
	s >> fp.julia;
	s >> fm;
	fp.fm = (formula)fm;
	QVector<uint32_t> vcx, vcy, vwidth, vparam, vparamq, vcrit;
	s >> vcx >> vcy >> vwidth >> vparam >> vparamq >> vcrit;
	fp.center_x = vcx.toStdVector ();
	fp.center_y = vcy.toStdVector ();
	fp.width = vwidth.toStdVector ();
	fp.param_p = vparam.toStdVector ();
	fp.param_q = vparamq.toStdVector ();
	fp.critpoint = vcrit.toStdVector ();
	s >> fp.nwords >> fp.power >> fp.maxiter;
	if (version >= 2) {
		qint32 dummy;
		s >> fp.hybrid_code >> fp.hybrid_len >> dummy;
	}

	set (fp.matrix[0][0], 1);
	set (fp.matrix[0][1], 0);
	set (fp.matrix[1][0], 0);
	set (fp.matrix[1][1], 1);

	return s;
}
#pragma GCC diagnostic pop

static int power_from_fp (frac_params &fp)
{
	int power = fp.power;
	if (fp.fm == formula::sqtwice_a || fp.fm == formula::sqtwice_b)
		power += 2;
	return power;
}

void MainWindow::layout_stored_params ()
{
	QGraphicsView *view = ui->storedView;
	if (m_stored.size () == 0) {
		ui->action_BatchRender->setEnabled (false);
		return;
	}
	if (!view->isVisible ())
		return;

	int w = view->width ();
	int h = view->height ();
	int extra_height = 0;

	int n_elts = m_stored.size ();

	int min_size = 220;
	int min_w = min_size;
	int min_h = min_size;
	min_h += extra_height;

	int sb_size = qApp->style()->pixelMetric (QStyle::PM_ScrollBarExtent);

	int smaller = w < h ? w : h;
	int smaller_min = w < h ? min_w : min_h;
	int larger = w < h ? h : w;
	int larger_min = w < h ? min_h : min_w;

	int new_width = 0;
	int new_height = 0;

	/* Two attempts: try to fit everything without scrollbars first, then reduce the
	   available area to account for the scroll bar.  */
	for (int assume_sb = 0; assume_sb < 2; assume_sb++) {
		/* The number of minimally-sized previews we can fit into the longer side.  */
		int n_larger_max = std::max (1, larger / larger_min);
		/* The number of minimally-sized previews we can fit on the smaller side.  */
		int n_smaller_max = std::max (1, smaller / smaller_min);
		/* The number of previews that need to fit on the smaller side to make sure that the
		   above number is enough to fit all on screen.  */
		int n_small_1 = std::max (1, (n_elts + n_larger_max - 1) / n_larger_max);
		/* How many previews to actually fit on the smaller side: the minimum of the two
		   previously calculated numbers (so ideally the number we want, but not exceeding
		   the maximum).  */
		int n_small = std::min (n_small_1, n_smaller_max);
		int n_large = (n_elts + n_small - 1) / n_small;
		/* Compute picture sizes for filling the area.  */
		int ps1 = smaller / n_small;
		int ps2 = larger / n_large;
		// printf ("ps1 %d [%d %d] ps2 %d [%d %d] %% %d\n", ps1, smaller, n_small, ps2, larger, n_large, sb_size);
		if (w < h)
			ps2 -= extra_height;
		else
			ps1 -= extra_height;
		int min_ps = ps1 < min_size ? ps1 : std::min (ps1, ps2);
		if (min_ps >= min_size || (ps1 < min_size && (ps2 > ps1 || assume_sb))) {
			/* Everything fits (or can't fit, if the window is too narrow).
			   But try to find an even better layout.  */
			while (n_small < n_elts) {
				int n_small_new = n_small + 1;
				int n_large_new = (n_elts + n_small_new - 1) / n_small_new;
				int ps1b = smaller / n_small_new;
				int ps2b = larger / n_large_new;
				if (w < h)
					ps2b -= extra_height;
				else
					ps1b -= extra_height;
				int t = std::min (ps1b, ps2b);
				if (t < min_ps)
					break;
				min_ps = t;
				n_small++;
			}
			new_width = min_ps;
			new_height = min_ps + extra_height;
			break;
		} else {
			/* Use scrollbars and use the maximum width of the smaller side.  */
			if (assume_sb == 0) {
				smaller -= sb_size;
				continue;
			}
			// ps1 = std::max (ps1, min_size);
			new_width = ps1;
			new_height = ps1 + extra_height;
			break;
		}
	}

	int row = 0;
	int col = 0;
	int max_row = 0;
	int max_col = 0;
	for (auto &p: m_stored) {
		p.x = col;
		p.y = row;
		max_row = std::max (max_row, row);
		max_col = std::max (max_col, col);
		if (w < h) {
			col += new_width;
			if (col + new_width > w) {
				col = 0;
				row += new_height;
			}
		} else {
			row += new_height;
			if (row + new_height > h) {
				row = 0;
				col += new_width;
			}
		}
	}
	view->setSceneRect (0, 0, max_col + new_width, max_row + new_height);

	m_stored_canvas.clear ();
	for (size_t i = 0; i < m_stored.size (); i++) {
		auto &p = m_stored[i];
		auto callback = [this] (QGraphicsSceneMouseEvent *) -> void {  };
		auto menu_callback = [this, &p, i] (QGraphicsSceneContextMenuEvent *e) -> bool
		{
			QMenu menu;
			menu.addAction (tr ("Go here"), [this, &p] () { restore_params (p.params.fp); });
			menu.addAction (tr ("Forget"),
					[this, i] () {
						m_stored_canvas.clear ();
						m_stored.erase (m_stored.begin () + i);
						layout_stored_params ();
					});
			menu.exec (e->screenPos ());

			return true;
		};
		constexpr int margin = 8;
		p.pixmap = new ClickablePixmap (QPixmap::fromImage (p.thumbnail.scaled (QSize { new_width - margin, new_height - margin })),
						callback, menu_callback);
		m_stored_canvas.addItem (p.pixmap);
		int szdiff_x = margin;
		int szdiff_y = margin;
		p.pixmap->setPos (p.x + szdiff_x / 2, p.y + szdiff_y / 2);
	}
}

GPU_handler *gpu_handler;

void MainWindow::start_threads ()
{
	m_gpu_thread = new QThread;
	m_gpu_thread->start (QThread::LowPriority);
	gpu_handler = new GPU_handler;

	gpu_handler->moveToThread (m_gpu_thread);
	connect (m_gpu_thread, &QThread::finished, gpu_handler, &QObject::deleteLater);

	connect (gpu_handler, &GPU_handler::signal_new_data, this, &MainWindow::slot_new_data);
	connect (gpu_handler, &GPU_handler::signal_kernel_complete, this, &MainWindow::slot_kernel_complete);
	connect (this, &MainWindow::signal_init_cuda, gpu_handler, &GPU_handler::slot_init_cuda);
	connect (this, &MainWindow::signal_alloc_mem, gpu_handler, &GPU_handler::slot_alloc_mem);
	connect (this, &MainWindow::signal_invalidate, gpu_handler, &GPU_handler::slot_invalidate);
	connect (this, &MainWindow::signal_start_kernel, gpu_handler, &GPU_handler::slot_start_kernel);
	connect (this, &MainWindow::signal_compile_kernel, gpu_handler, &GPU_handler::slot_compile_kernel);
	QString errstr;
	emit signal_init_cuda (&errstr);
	gpu_handler->done_sem.acquire ();
	if (!errstr.isEmpty ()) {
		QMessageBox::critical (this, PACKAGE, errstr);
		close ();
	}

	m_render_thread = new QThread;
	m_render_thread->start ();
	m_renderer = new Renderer;
	m_renderer->moveToThread (m_render_thread);
	connect (m_render_thread, &QThread::finished, m_renderer, &QObject::deleteLater);
	connect (m_renderer, &Renderer::signal_render_complete, this, &MainWindow::slot_render_complete);
	connect (this, &MainWindow::signal_render, m_renderer, &Renderer::slot_render);

	m_preview_thread = new QThread;
	m_preview_thread->start ();
	m_preview_renderer = new Renderer;
	m_preview_renderer->moveToThread (m_preview_thread);
	connect (m_preview_thread, &QThread::finished, m_preview_renderer, &QObject::deleteLater);
	connect (m_preview_renderer, &Renderer::signal_render_complete, this, &MainWindow::slot_render_complete);
	connect (this, &MainWindow::signal_render_preview, m_preview_renderer, &Renderer::slot_render);
}

void MainWindow::build_points (frac_desc &fd, int w, int h)
{
	vpvec step = div1 (fd.width, std::min (w, h) * fd.samples);
	fd.step = step;
}

static double sin_deg (int angle)
{
	switch (angle) {
	case 30: case 150:
		return 0.5;
	case 45: case 135:
		return M_SQRT1_2;
	case 90:
		return 1;
	case 270:
		return -1;
	case 225: case 315:
		return -M_SQRT1_2;
	case 210: case 330:
		return -0.5;
	}
	return sin (angle * M_PI / 180);
}

static double cos_deg (int angle)
{
	switch (angle) {
	case 60: case 300:
		return 0.5;
	case 45: case 315:
		return M_SQRT1_2;
	case 90: case 270:
		return 0;
	case 135: case 225:
		return -M_SQRT1_2;
	case 120: case 240:
		return -0.5;
	}
	return cos (angle * M_PI / 180);
}

void MainWindow::set_rotation (frac_desc &fd, double angle)
{
	angle = fmod (angle, 360);
	fd.rotation_angle = angle;
	double shear = shear_slider_value ();
	double scale = scale_slider_value ();
	int shidx = ui->shearComboBox->currentIndex ();
	int scidx = ui->scaleComboBox->currentIndex ();
	double c = angle == floor (angle) ? cos_deg (angle) : cos (angle * M_PI / 180);
	double s = angle == floor (angle) ? sin_deg (angle) : sin (angle * M_PI / 180);
	double xshear = shidx == 1 ? shear : 0;
	double yshear = shidx == 2 ? shear : 0;
	double xscale = scidx == 1 ? scale : 1;
	double yscale = scidx == 2 ? scale : 1;

	set (fd.matrix[0][0], (c - s * yshear) * xscale);
	set (fd.matrix[0][1], (c * xshear - s) * yscale);
	set (fd.matrix[1][0], (s + c * yshear) * xscale);
	set (fd.matrix[1][1], (s * xshear + c) * yscale);
}

void MainWindow::inc_rotation (frac_desc &fd, int angle)
{
	set_rotation (fd, fd.rotation_angle + angle);
}

void MainWindow::enter_rotation (bool)
{
	auto &fd = current_fd ();
	RotationDialog dlg (this, fd.rotation_angle);
	connect (&dlg, &RotationDialog::apply_rotation,
		 [this, &fd] (double v) { printf ("rotate %f\n", v);set_rotation (fd, v); update_settings (false); });
	dlg.exec ();
}

double MainWindow::shear_slider_value ()
{
	double maxv = ui->shearSlider->maximum ();
	double v = ui->shearSlider->value ();
	v /= maxv;
	/* Arrange for finer control near the center.  */
	v = 1 - cbrt (1 - v);
	/* Scale to the actual maximum we want to support.  */
	v *= 20;
	return v;
}

double MainWindow::scale_slider_value ()
{
	double maxv = ui->scaleSlider->maximum ();
	double v = ui->scaleSlider->value ();
	/* Force into range 0..1.  */
	v = (v - 1) / (maxv - 1);
	/* Arrange for finer control near the center.  */
	v = v * v * v;
	/* Scale to the actual maximum we want to support.  */
	v = 1 + v * 20;
	return v;
}

void MainWindow::discard_fd_data (frac_desc &fd)
{
	fd.generation++;
	// Do this always even if there is nothing to free.
	// Calling invalidate and Waiting on done_sem also ensures that we are not within
	// slot_start_kernel.
	emit signal_invalidate (&fd);
	gpu_handler->done_sem.acquire ();

	if (fd.n_pixels == 0)
		return;

	delete[] fd.host_cplxvals;
	delete[] fd.host_zprev;
	delete[] fd.host_zpidx;
	delete[] fd.host_result;
	delete[] fd.host_coords;

	delete[] fd.pic_zprev;
	delete[] fd.pic_t;
	delete[] fd.pic_result;
	delete[] fd.pic_iter_value;
	fd.n_pixels = 0;
}

void MainWindow::compute_fractal (frac_desc &fd, int nwords, int n_prev, int w, int h, int full_h,
				  int maxiter, int ss, bool isdem, bool preview, bool batch)
{
	QMutexLocker render_lock (&m_renderer->mutex);
	QMutexLocker preview_lock (&m_preview_renderer->mutex);

	ss = 1 << ss;

	// We need to ensure nthreads is always big enough to handle the initial_setup phase.
	int nthreads = w * h;
	int npixels = w * h * ss * ss;
	int nvals = n_formula_cplx_vals (fd.fm, isdem);
	if (fd.nvals_allocated != nvals || fd.samples != ss
	    || fd.n_pixels != npixels
	    || fd.n_threads != nthreads
	    || fd.nwords != nwords
	    || fd.n_prev != n_prev
	    || fd.dem != isdem)
	{
		discard_fd_data (fd);

		fd.nvals_allocated = nvals;
		fd.dem = isdem;
		fd.pixel_width = w * ss;
		fd.pixel_height = h * ss;
		fd.n_pixels = npixels;
		fd.n_threads = nthreads;
		fd.nwords = nwords;
		fd.n_prev = n_prev;
		fd.samples = ss;

		int nvals = n_formula_cplx_vals (fd.fm, isdem);
		fd.host_cplxvals = new uint32_t[nwords * 2 * nvals * nthreads];
		fd.host_zprev = new double[n_prev * 2 * nthreads];
		fd.host_zpidx = new uint32_t[nthreads];
		fd.host_coords = new uint32_t[nthreads];
		fd.host_result = new uint32_t[nthreads];

		fd.pic_t = nullptr;
		if (n_prev > 1)
			fd.pic_t = new double[2 * npixels];
		fd.pic_zprev = new double[2 * n_prev * npixels];
		fd.pic_zder = nullptr;
		if (fd.dem)
			fd.pic_zder = new double[2 * npixels];
		fd.pic_result = new uint32_t[npixels];
		fd.pixels_done = bit_array (npixels);
		fd.pixels_started = bit_array (npixels);
		fd.pic_iter_value = nullptr;
	}
	fd.full_height = full_h * ss;
	fd.cmin = 10000;
	fd.cmax = -10000;
	if (!batch) {
		fd.min_stripeval = 0;
		fd.max_stripeval = 1;
	}
	fd.pixel_step = preview || batch ? 1 : ss * 2;
	fd.n_completed = 0;
	fd.start_idx = 0;
	fd.pixels_done.clear ();
	fd.pixels_started.clear ();
	memset (fd.pic_result, 0, npixels * sizeof (uint32_t));
	if (fd.pic_iter_value != nullptr)
		memset (fd.pic_iter_value, 0, npixels * sizeof (double));
	if (fd.dem)
		memset (fd.pic_zder, 0, 2 * npixels * sizeof (double));
	fd.maxiter = preview ? iter_steps : maxiter;
	fd.maxiter_found = 0;
	build_points (fd, w, full_h);

	QString errstr;
	emit signal_alloc_mem (&fd, max_nwords, nwords, w, h, &errstr);
	gpu_handler->done_sem.acquire ();
	if (!errstr.isEmpty ()) {
		QMessageBox::critical (this, PACKAGE, tr ("Could not allocate CUDA memory: ") + errstr);
		return;
	}

	m_generation++;
	fd.generation = m_generation;
	if (m_working)
		abort ();
	m_working = !batch;
	gpu_handler->processing_data = false;

	emit signal_start_kernel (&fd, m_generation, max_nwords, iter_steps, batch);
}

void MainWindow::enable_sac_or_tia ()
{
	bool sac = ui->action_SAC->isChecked ();
	bool tia = ui->action_TIA->isChecked ();
	bool on = sac || tia;
	ui->stripesWidget->setEnabled (sac);
	if (!on) {
		ui->superBox->setEnabled (true);
		n_prev_requested = 1;
		update_views ();
		return;
	}
	QSettings settings;
	n_prev_requested = 1 << settings.value ("coloring/nprev").toInt ();
	bool safety = settings.value ("coloring/nosuper-sac").toBool ();
	ui->superBox->setEnabled (!safety);
	if (safety)
		ui->sampleSpinBox->setValue (0);
	update_settings (false);
}

/* Called whenever an angle group menu item other than SAC is chosen.  */
void MainWindow::slot_disable_sac (bool on)
{
	if (n_prev_requested != 1) {
		n_prev_requested = 1;
		update_settings (false);
	} else
		update_views ();
}

void MainWindow::closeEvent (QCloseEvent *event)
{
	QSettings settings;
	settings.setValue ("mainwin/geometry", saveGeometry ());
	settings.setValue ("mainwin/windowState", saveState ());

	m_generation++;
	m_paused = true;
	abort_computation ();
	emit signal_invalidate (&m_fd_mandel);
	gpu_handler->done_sem.acquire ();
	emit signal_invalidate (&m_fd_julia);
	gpu_handler->done_sem.acquire ();

	QMainWindow::closeEvent (event);
}

void MainWindow::restore_geometry ()
{
	QSettings settings;
	restoreGeometry (settings.value("mainwin/geometry").toByteArray());
	restoreState (settings.value("mainwin/windowState").toByteArray());
}

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

#if 0
static int col_sum (uint32_t col)
{
	return (col & 255) + ((col >> 8) & 255) + (col >> 16);
}
#endif
static inline QRgb color_from_niter (const QVector<uint32_t> &palette, double niter, int type, double steps,
				     int slider)
{
	if (type == 1) {
		niter = sqrt (niter);
	} else if (type == 2) {
		niter = cbrt (niter);
	} else if (type == 3) {
		niter = pow (niter, 1/4.);
	} else if (type == 4) {
		niter = log (niter);
	} else if (type == 5) {
		niter = log (sqrt (niter));
	} else if (type == 6) {
		niter = log (log (niter));
	}
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

static inline uint32_t modify_color (uint32_t col, double v)
{
	QColor c = QColor::fromRgb (col);
	int cl = c.value ();
	// Lighten dark colors, darken light colors
	v -= (1 - cl / 255.) / 2;
#if 0
	/* Both HSV and HSL introduce ugly artifacts.  */
	int ch = c.hsvHue ();
	int cs = c.hsvSaturation ();
	if (v < 0) {
		cl = 255 - (255 - cl) * (v + 1);
	} else {
		cl *= 1 - v;
	}
	col = QColor::fromHsv (ch, cs, cl).rgb ();
#else
	int cr = c.red ();
	int cg = c.green ();
	int cb = c.blue ();
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
#endif
	return col;
}

inline std::pair<double, int> iter_value_at (frac_desc *fd, int idx, int power)
{
	double v = fd->pic_result[idx];
	if (v == 0)
		return { 0, 0 };
	double re2 = fd->pic_zprev[idx * 2 * fd->n_prev];
	double im2 = fd->pic_zprev[idx * 2 * fd->n_prev + 1];
	re2 *= re2;
	im2 *= im2;
	// Some attractor other than infinity.
	if (re2 + im2 < 16)
		return { v, 1 };
	double radius = fd->dem ? 10 : 100;
	double correction = log (0.5 * log ((double)re2 + im2) / log (radius)) / log (power);
	fd->cmin = std::min (fd->cmin, correction);
	fd->cmax = std::max (fd->cmax, correction);
	return { v + 5 - correction, 0 };
}

class runner : public QRunnable
{
	QSemaphore *completion_sem;
	std::atomic<bool> *success;
	const render_params &rp;
	int w, y0, y0e;
	frac_desc *fd;
	double minimum;
	QRgb *data;
	double dstep;
	int power;
	double found_stripe_min = 1.0;
	double found_stripe_max = 0.0;
	double *new_stripe_min, *new_stripe_max;
	QMutex *update_mutex;
public:
	runner (QSemaphore *sem, std::atomic<bool> *succ_in, const render_params &rp_in,
		int w_in, int y0_in, int y0e_in, frac_desc *fd_in, double min_in, QRgb *data_in,
		double *nsmin, double *nsmax, QMutex *umutex)
		: completion_sem (sem), success (succ_in),
		  rp (rp_in), w (w_in), y0 (y0_in), y0e (y0e_in), fd (fd_in), minimum (min_in), data (data_in),
		  new_stripe_min (nsmin), new_stripe_max (nsmax), update_mutex (umutex)
	{
		setAutoDelete (true);
	}

	void compute_color (size_t idx, int &r, int &g, int &b, int &alpha, int &outcolor)
	{
		int n_prev = fd->n_prev;
		double re = fd->pic_zprev[idx * 2 * n_prev];
		double im = fd->pic_zprev[idx * 2 * n_prev + 1];

		auto [ v, attractor ] = iter_value_at (fd, idx, power);
		double v1 = v;
		if (v != 0) {
			outcolor++;

			if (rp.sub)
				v -= minimum - rp.sub_val;
			uint32_t col = color_from_niter (rp.palette, v, rp.mod_type, rp.steps, rp.slider);

			if (rp.sac && attractor == 0) {
				double density = rp.sac_factor;
				double re2 = re * re;
				double im2 = im * im;
				int thisnp = std::min ((uint32_t)n_prev, fd->pic_result[idx]);
				double firstval = 0.5 * sin (density * atan2 (im, re)) + 0.5;
				double sum = 0;
				for (int i = 1; i < thisnp; i++) {
					double lastre = fd->pic_zprev[idx * 2 * n_prev + i * 2];
					double lastim = fd->pic_zprev[idx * 2 * n_prev + i * 2 + 1];
					double val = 0.5 * sin (density * atan2 (lastim, lastre)) + 0.5;
					sum += val;
				}
				double avg1 = thisnp > 1 ? sum / (thisnp - 1) : 0;
				double avg2 = (sum + firstval) / thisnp;
				double radius = fd->dem ? 10 : 100;
				double mixfactor = log (0.5 * log (re2 + im2) / log (radius)) / log (power);
				double mod = avg1 * mixfactor + avg2 * (1 - mixfactor);
				found_stripe_min = std::min (mod, found_stripe_min);
				found_stripe_max = std::max (mod, found_stripe_max);
				double prev_width = fd->max_stripeval - fd->min_stripeval;
				if (rp.sac_contrast)
					mod = std::max (0.0, std::min (1.0, (mod - fd->min_stripeval) / prev_width));
				if (rp.angle_colour)
					col = modify_color (col, mod);
				else {
					col = 0x01010101 * floor (mod * 255);
				}
			} else if (rp.tia && fd->pic_t != nullptr && attractor == 0) {
				double cre = fd->pic_t[idx * 2];
				double cim = fd->pic_t[idx * 2 + 1];
				double re2 = re * re;
				double im2 = im * im;
				double radius = fd->dem ? 10 : 100;
				double mixfactor = log (0.5 * log (re2 + im2) / log (radius)) / log (power);
				double cmag = sqrt (cre * cre + cim * cim);
				int thisnp = std::min ((uint32_t)n_prev, fd->pic_result[idx]);
				if (cmag == 0)
					thisnp = 0;
				double sum = 0;
				double firstval = 0;
				for (int i = 1; i < thisnp; i++) {
					double lastre = fd->pic_zprev[idx * 2 * n_prev + i * 2];
					double lastim = fd->pic_zprev[idx * 2 * n_prev + i * 2 + 1];
					double mag = sqrt (re2 + im2);
                                        double zlmag = pow (lastre * lastre + lastim * lastim, power / 2.0);
                                        mag -= zlmag - cmag;
					mag /= 2 * cmag;
					if (i == 1)
						firstval = mag;
					else
						sum += mag;
                                        re2 = lastre * lastre;
                                        im2 = lastim * lastim;
				}
				double avg1 = thisnp > 2 ? sum / (thisnp - 2) : 0;
				double avg2 = thisnp > 1 ? (sum + firstval) / (thisnp - 1) : 0;
				double mod = avg1 * mixfactor + avg2 * (1 - mixfactor);
				if (mod < 0 || mod > 1)
					// ??? This isn't supposed to happen, but let's not let it mess up
					// the min/max values if it does happen.
					col = rp.angle_colour ? col : 0;
				else {
					found_stripe_min = std::min (mod, found_stripe_min);
					found_stripe_max = std::max (mod, found_stripe_max);
					double prev_width = fd->max_stripeval - fd->min_stripeval;
					if (rp.sac_contrast)
						mod = std::max (0.0, std::min (1.0, (mod - fd->min_stripeval) / prev_width));
					if (rp.angle_colour)
						col = modify_color (col, mod);
					else {
						col = 0x01010101 * floor (mod * 255);
					}
				}
			}
			else if (rp.angle == 1) {
				double re2 = re * re;
				double im2 = im * im;
				double size = sqrt (re2 + im2);
				double init_angle = atan2 (re, im);
				init_angle += M_PI * 2;
#if 0
				double compare_angle = M_PI * v1 * 2;
				compare_angle -= M_PI * 2 * floor (compare_angle / (M_PI * 2));
				init_angle = init_angle - compare_angle + M_PI * 2;
#endif
				double ang = init_angle - M_PI * 2 * floor (init_angle / (M_PI * 2));
#if 0
				double v = fabs (sin (ang)) / 2;
#else
				ang += M_PI / 2;
				ang = fmod (ang, M_PI * 2);
				double v = fabs (ang / (M_PI * 2) - 0.5);
#endif
				col = modify_color (col, v);
			}
			else if (rp.angle == 2) {
				if (im < 0)
					col = 0;
				else if (!rp.angle_colour)
					col = 0xFFFFFFFF;
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
		update_mutex->lock ();
		*new_stripe_min = std::min (*new_stripe_min, found_stripe_min);
		*new_stripe_max = std::max (*new_stripe_max, found_stripe_max);
		update_mutex->unlock ();
		// Must be outside the lock, otherwise the lock could go away from under us.
		completion_sem->release ();
	}
};

/* The view is null if we are being called for a batch render.  */
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
			auto [ v, attractor ] = iter_value_at (fd, i, power);
			v = floor (v);
			if (v > 0 && (minimum == 0 || v < minimum))
				minimum = v;
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
	double new_stripe_min = 1;
	double new_stripe_max = 0;
	QMutex update_mutex;
	for (int y0 = 0; y0 < h; y0 += lines_per_thread) {
		int y0e = y0 + lines_per_thread;
		if (y0e > h)
			y0e = h;
		m_pool.start (new runner (&completion_sem, &any_found, rp, w, y0, y0e, fd, minimum, data + w * y0,
					  &new_stripe_min, &new_stripe_max, &update_mutex));
		n_started++;
	}
	completion_sem.acquire (n_started);
	if (!any_found)
		return;

	if (view != nullptr) {
		fd->min_stripeval = new_stripe_min;
		fd->max_stripeval = new_stripe_max;
		emit signal_render_complete (view, fd, result_image, minimum);
	}
}

void Renderer::slot_render (frac_desc *fd, QGraphicsView *view, int generation)
{
	queue_mutex.lock ();
	assert (queued);
	queued = false;
	render_params this_p = next_rp;
	int w = render_width;
	int h = render_height;
	queue_mutex.unlock ();
	result_image = QImage (w, h, QImage::Format_RGB32);
	do_render (this_p, w, h, 0, fd, view, generation);
}

void MainWindow::slot_render_complete (QGraphicsView *view, frac_desc *fd, QImage img, double minimum)
{
	if (fd->julia) {
		m_img_julia = img;
		m_min_julia = minimum;
	} else {
		m_img_mandel = img;
		m_min_mandel = minimum;
	}
	QPixmap pm = QPixmap::fromImage (img);
	QGraphicsScene *canvas = view->scene ();
	canvas->clear ();
	canvas->addItem (new QGraphicsPixmapItem (pm));

//	ui->fractalView->setBackgroundBrush (img);
//	ui->fractalView->setCacheMode(QGraphicsView::CacheBackground);
}

void MainWindow::set_render_params (render_params &p)
{
	p.palette = m_palette;
	p.incol = ui->action_ICWhite->isChecked () ? 0xFFFFFF : 0;
	p.mod_type = ui->modifyComboBox->currentIndex ();
	int steps_spin = ui->widthSpinBox->value ();
	p.steps = (pow (steps_spin + 1, 1.3) - 2) / 2;
	p.angle = !!ui->action_AngleSmooth->isChecked () + 2 * !!ui->action_AngleBin->isChecked ();
	p.tia = ui->action_TIA->isChecked ();
	p.sac = ui->action_SAC->isChecked ();
	p.sac_factor = ui->stripesSpinBox->value ();
	p.sac_contrast = ui->action_Contrast->isChecked ();
	p.sub = !ui->action_ShiftNone->isChecked ();
	p.sub_val = ui->action_Shift10->isChecked () ? 10 : ui->action_Shift100->isChecked () ? 100 : 0;
	p.slider = ui->colStepSlider->value ();
	p.dem = ui->action_DEMDist->isChecked () || ui->action_DEMBoth->isChecked ();
	p.dem_shade = ui->action_DEMShading->isChecked () || ui->action_DEMBoth->isChecked ();
	p.dem_colour = ui->action_DEMColour->isChecked ();
	p.angle_colour = ui->action_AngleColour->isChecked ();
	p.dem_param = ui->demParamSpinBox->value () * (1 << ui->sampleSpinBox->value ());
	p.dem_strength = ui->demStrengthSpinBox->value ();
	p.dem_start = m_dem_start;
	p.dem_stop = m_dem_stop;
	p.aspect = chosen_aspect ();
}

void MainWindow::update_display (QGraphicsView *view)
{
	if (!view->isVisible ())
		return;

	frac_desc &fd = view == ui->previewView ? m_fd_julia : current_fd ();
	if (fd.n_completed == 0)
		return;

	Renderer *r = view == ui->previewView ? m_preview_renderer : m_renderer;
	{
		QMutexLocker lock (&r->queue_mutex);
		set_render_params (r->next_rp);
		r->render_width = view->width ();
		r->render_height = view->height ();
		if (r->queued)
			return;

		r->queued = true;
	}
	if (r == m_renderer)
		emit signal_render (&fd, view, m_generation);
	else
		emit signal_render_preview (&fd, view, m_generation);
}

void MainWindow::autoprec (frac_desc &fd)
{
	if (!ui->action_AutoPrec->isChecked ())
		return;

	bool_changer (m_inhibit_updates, true);

	bool dem = !ui->action_DEMOff->isChecked ();
	int w = ui->fractalView->width ();
	int h = ui->fractalView->height ();
	vpvec step = div1 (fd.width, std::max (w, h));
	int i = max_nwords;
	while (i-- > 0)
		if (step[i] != 0)
			break;
	if (i > 0 && step[i] < (dem ? 128 : 32))
		i--;
	int required = max_nwords - i;
	if ((fd.fm == formula::facing || fd.fm == formula::facing_b) && required < max_nwords)
		required++;
	if (required > m_nwords) {
		m_nwords = required;
		ui->precSpinBox->setValue (m_nwords);
		m_recompile = true;
	}
}

void MainWindow::render_fractal ()
{
	if (!isVisible ())
		return;

	frac_desc &fd = current_fd ();

	int w = ui->fractalView->width ();
	int h = ui->fractalView->height ();
	m_canvas.setSceneRect (0, 0, w, h);
	printf ("render_fractal size %d %d\n", w, h);

	bool isdem = !ui->action_DEMOff->isChecked ();
	compute_fractal (fd, m_nwords, n_prev_requested, w, h, h, default_maxiter, ui->sampleSpinBox->value (), isdem, false);
}

void MainWindow::render_preview ()
{
	if (!ui->previewView->isVisible ())
		return;

	int w = ui->previewView->width ();
	int h = ui->previewView->height ();

	if (w == 0 || h == 0)
		return;

	m_preview_canvas.setSceneRect (0, 0, w, h);
	printf ("render_preview size %d %d\n", w, h);

	bool isdem = !ui->action_DEMOff->isChecked ();
	compute_fractal (m_fd_julia, m_nwords, n_prev_requested, w, h, h, default_maxiter, 0, isdem, true);
}

void MainWindow::slot_new_data (frac_desc *fd, int generation, bool success)
{
	if (!success) {
		m_working = false;
		QMessageBox::warning (this, PACKAGE, tr ("Failed to run CUDA kernel."));
		return;
	}

	if (generation != m_generation)
		return;

	double percent = (double)(fd->n_completed + 1) * 100 / (fd->n_pixels + 1);
	setWindowTitle (PACKAGE " - " + QString::number (percent, 'f', 2) + "% complete");

	if (fd->n_completed == 0)
		return;

	{
		QMutexLocker lock (&gpu_handler->data_mutex);

		gpu_handler->processing_data = true;
		gpu_handler->data_available = false;
	}

	if (fd != &current_fd ()) {
		m_preview_uptodate = true;
		update_display (ui->previewView);
	} else
		update_display (ui->fractalView);

	{
		QMutexLocker lock (&gpu_handler->data_mutex);
		bool more_data = gpu_handler->data_available && !m_recompile && !m_reinit_render;
		if (!more_data) {
			gpu_handler->processing_data = false;
			return;
		}
	}

	QMetaObject::invokeMethod (this,
				   [=] () { slot_new_data (fd, generation, true); },
				   Qt::QueuedConnection);
}

void MainWindow::do_compile ()
{
	QMutexLocker render_lock (&m_renderer->mutex);
	QMutexLocker preview_lock (&m_preview_renderer->mutex);

	QString errstr;
	emit signal_invalidate (&m_fd_mandel);
	gpu_handler->done_sem.acquire ();
	emit signal_invalidate (&m_fd_julia);
	gpu_handler->done_sem.acquire ();

	int power = ui->powerSpinBox->value ();
	auto it = std::find (std::begin (formula_table), std::end (formula_table), m_formula);
	int fidx = it - std::begin (formula_table);
	emit signal_compile_kernel (fidx, power, m_nwords, max_nwords, &errstr);

	m_power = power;
	gpu_handler->done_sem.acquire ();
	if (!errstr.isEmpty ()) {
		QMessageBox::critical (this, PACKAGE, errstr);
		close ();
		return;
	}
	m_fd_mandel.power = power;
	m_fd_julia.power = power;
	m_fd_mandel.fm = m_formula;
	m_fd_julia.fm = m_formula;

	m_recompile = false;
	m_preview_uptodate = false;
	m_reinit_render = true;
}

void MainWindow::slot_kernel_complete ()
{
	m_working = false;

	if (m_paused)
		return;

	if (m_recompile)
		do_compile ();

	if (!current_fd ().julia && !m_preview_uptodate && ui->previewView->isVisible ()) {
		render_preview ();
	} else if (m_reinit_render) {
		m_reinit_render = false;
		render_fractal ();
	} else {
		frac_desc &fd = current_fd ();

		if (fd.n_completed < fd.n_pixels) {
			m_working = true;
			emit signal_start_kernel (&fd, m_generation, max_nwords, iter_steps, false);
		}
	}

}

bool MainWindow::abort_computation ()
{
	QMutexLocker lock (&gpu_handler->data_mutex);
	if (m_working) {
		// printf ("forcing abort\n");
		gpu_handler->abort_computation = true;
		return true;
	}
	return false;
}

void MainWindow::restart_computation ()
{
	if (abort_computation ())
		/* slot_kernel complete will be signalled from the GPU handler thread.  */
		return;

	slot_kernel_complete ();
}

void MainWindow::preview_wheel_event (QWheelEvent *e)
{
	if (e->angleDelta().y() == 0)
		return;

	frac_desc &fd = m_fd_julia;

	int factor = ui->zoomSpinBox->value () * 10;
	if (e->angleDelta ().y () < 0)
		fd.width = div1 (mul1 (fd.width, factor), 10);
	else
		fd.width = div1 (mul1 (fd.width, 10), factor);
	fd.bounds_w = fd.bounds_h = 0;

	m_preview_uptodate = false;
	restart_computation ();
}

void MainWindow::fractal_wheel_event (QWheelEvent *e)
{
	if (e->angleDelta().y() == 0)
		return;
	if (e->angleDelta().y() < 0)
		zoom_out ();
	else
		zoom_in ();
}

/* Try to adjust clicks so that they land on "integer" points, where we extend
   the definition of "integer" slightly to all multiples of 1/4.
   The idea is to allow the picture (or the Julia parameter)to be centered on
   potentially interesting values - like the real line for symmetry.  */

static vpvec aim_assist_1 (vpvec v, vpvec step)
{
	size_t sz = v.size ();
	vpvec range = mul1 (step, 16);
	vpvec l = sub (v, range);
	vpvec h = add (v, range);
	if (l[sz - 1] == h[sz - 1] && (l[sz - 2] & 0xC0000000) == (h[sz - 2] & 0xC0000000))
		return v;
	for (int i = 0; i < sz - 2; i++)
		h[i] = 0;
	h[sz - 2] &= 0xC0000000;
	return h;
}

void MainWindow::fractal_mouse_event (QMouseEvent *e)
{
	if (e->type () != QEvent::MouseButtonPress)
		return;

	auto pos = e->pos ();
	auto scene_pos = ui->fractalView->mapToScene (pos);
	frac_desc &fd = current_fd ();

	bool shf = (e->modifiers () & Qt::ShiftModifier) != 0;
	bool param = (e->modifiers () & Qt::ControlModifier) != 0;
	// printf ("click %d %d\n", (int)p.x (), (int)p.y ());

	int w = ui->fractalView->width ();
	int h = ui->fractalView->height ();

	int px = scene_pos.x ();
	int py = scene_pos.y ();
	double pxn = px - w / 2;
	double pyn = py - h / 2;
	int real_px = to_double (fd.matrix[0][0]) * pxn + to_double (fd.matrix[0][1]) * pyn;
	int real_py = to_double (fd.matrix[1][0]) * pxn + to_double (fd.matrix[1][1]) * pyn;
	auto pxstep = mul1 (fd.step, fd.samples * abs (real_px));
	auto pystep = mul1 (fd.step, fd.samples * abs (real_py));
	vpvec a = real_px < 0 ? sub (fd.center_x, pxstep) : add (fd.center_x, pxstep);
	vpvec b = real_py < 0 ? sub (fd.center_y, pystep) : add (fd.center_y, pystep);
	if (ui->action_AimAssist->isChecked ()) {
		a = aim_assist_1 (a, fd.step);
		b = aim_assist_1 (b, fd.step);
	}
	if (param) {
		memcpy (&m_fd_julia.param_p[0], &a[0], max_nwords * sizeof (uint32_t));
		memcpy (&m_fd_julia.param_p[max_nwords], &b[0], max_nwords * sizeof (uint32_t));

		if (ui->typeComboBox->currentIndex () == 1) {
			m_reinit_render = true;
			restart_computation ();
		} else if (ui->previewView->isVisible ()) {
			m_preview_uptodate = false;
			restart_computation ();
		}
	} else {
		if (!shf) {
			fd.center_x = a;
			fd.center_y = b;
		}
		abort_computation ();
		fd.bounds_w = fd.bounds_h = 0;
		fd.width = div1 (mul1 (fd.width, 2), 5);
		autoprec (fd);
		m_reinit_render = true;
		restart_computation ();
	}
}

/* Reset is true if we are changing the power, in which case we reset the coords.  */
void MainWindow::update_settings (bool reset)
{
	if (m_inhibit_updates)
		return;

	m_nwords = ui->precSpinBox->value ();
	int power = ui->powerSpinBox->value ();
	if (m_power != power || m_fd_julia.nwords != m_nwords || m_fd_mandel.nwords != m_nwords
	    || m_fd_julia.fm != m_formula || m_fd_mandel.fm != m_formula)
		m_recompile = true;
	m_reinit_render = true;
	if (reset) {
		reset_coords (m_fd_mandel);
		reset_coords (m_fd_julia);
	} else
		autoprec (current_fd ());
	restart_computation ();
}

void MainWindow::update_views (int)
{
	// No need to update in response to user inputs if we're going into slot_new_data again soon.
	if (!gpu_handler->processing_data) {
		update_display (ui->fractalView);
	}
	if (m_preview_uptodate)
		update_display (ui->previewView);
}

frac_desc &MainWindow::current_fd ()
{
	bool julia = ui->typeComboBox->currentIndex () == 1;
	return julia ? m_fd_julia : m_fd_mandel;
}

void MainWindow::reset_coords (frac_desc &fd)
{
	vpvec zero (max_nwords, 0);

	fd.width = zero;
	fd.bounds_w = fd.bounds_h = 4;
	fd.width[max_nwords - 1] = 4;

	fd.center_x = zero;
	fd.center_y = zero;

	int power = ui->powerSpinBox->value ();
	if (m_formula == formula::standard) {
		if (power == 2 && !fd.julia) {
			fd.center_x[max_nwords - 2] = 0x40000000;
			fd.center_x[max_nwords - 1] = 0xffffffff;
			fd.bounds_h = 3;
		} else if (power == 3 && !fd.julia) {
			fd.bounds_h = 3;
			fd.bounds_w = 2;
		} else if (power > 3) {
			fd.bounds_h = 3;
			fd.bounds_w = 3;
		}
	} else if (m_formula == formula::lambda) {
		if (!fd.julia) {
			if (power == 3)
				fd.bounds_h = fd.bounds_w = 2;
			else if (power > 3)
				fd.bounds_h = fd.bounds_w = 1;
		} else {
			if (power == 3) {
				fd.bounds_w = 6;
				fd.bounds_h = 4;
			} else
				fd.bounds_h = fd.bounds_w = 4;
		}
	} else if (m_formula == formula::mix && !fd.julia) {
		fd.center_x[max_nwords - 2] = 0xe0000000;
		fd.center_x[max_nwords - 1] = 0xffffffff;
		fd.width[max_nwords - 1] = 2;
		fd.bounds_w = 2;
		fd.bounds_h = 1;
	} else if (m_formula == formula::magnet_a) {
		if (!fd.julia)
			fd.center_x[max_nwords - 1] = 1;
		if (fd.julia)
			fd.bounds_h = fd.bounds_w = 10;
		else {
			fd.bounds_w = 5;
			fd.bounds_h = 6;
		}
		fd.width[max_nwords - 1] = fd.bounds_h;
	} else if (m_formula == formula::facing) {
		if (fd.julia)
			fd.bounds_h = fd.bounds_w = 3;
		else if (power == 2)
			fd.bounds_h = fd.bounds_w = 2;
		else
			fd.bounds_h = fd.bounds_w = 1;
		fd.width[max_nwords - 1] = fd.bounds_w;
	}
	if (m_formula == formula::standard && !fd.julia && power == 4) {
		fd.center_x[max_nwords - 2] = 0xd0000000;
		fd.center_x[max_nwords - 1] = 0xffffffff;
	} else if (m_formula == formula::lambda && !fd.julia && power == 4) {
		fd.center_x[max_nwords - 2] = 0xe8000000;
		fd.center_x[max_nwords - 1] = 0xffffffff;
	} else if (m_formula == formula::tricorn && !fd.julia) {
		fd.center_x[max_nwords - 2] = 0xb0000000;
		fd.center_x[max_nwords - 1] = 0xffffffff;
	}
	adjust_width_for_bounds (fd);

	bool_changer (m_inhibit_updates, true);
	ui->scaleComboBox->setCurrentIndex (0);
	ui->shearComboBox->setCurrentIndex (0);
	ui->shearSlider->setValue (0);
	ui->scaleSlider->setValue (1);
	set_rotation (fd, 0);
}

void MainWindow::set_q (int qr, int qi)
{
	vpvec val = vpvec (max_nwords * 2);
	val[max_nwords - 1] = qr;
	val[max_nwords * 2 - 1] = qi;
	m_fd_mandel.param_q = val;
	m_fd_julia.param_q = val;
}

void MainWindow::init_formula (formula f)
{
	vpvec cplx_zero = vpvec (max_nwords * 2);
	m_fd_mandel.critpoint = cplx_zero;
	m_fd_mandel.param_p = cplx_zero;
	m_fd_mandel.param_q = cplx_zero;
	m_fd_julia.critpoint = cplx_zero;
	m_fd_julia.param_p = cplx_zero;
	m_fd_julia.param_q = cplx_zero;

	m_fd_mandel.hybrid_len = 0;
	m_fd_julia.hybrid_len = 0;
	if (f == formula::lambda) {
		vpvec one = cplx_zero;
		one[max_nwords - 1] = 1;
		m_fd_mandel.critpoint = one;
	}
	if (f == formula::mix || f == formula::facing || f == formula::facing_b) {
		vpvec one = cplx_zero;
		one[max_nwords - 1] = 1;
		m_fd_mandel.critpoint = one;
	}
	if (f == formula::mix)
		set_q (2, 0);
}

void MainWindow::enable_interface_for_settings ()
{
	QSettings settings;
	bool largemem = settings.value ("largemem").toBool ();
	ui->action_SAC->setEnabled (largemem);
	ui->action_TIA->setEnabled (largemem);
	if (!largemem && (ui->action_SAC->isChecked () || ui->action_TIA->isChecked ()))
		ui->action_AngleNone->setChecked (true);
}

void MainWindow::enable_interface_for_formula (formula f)
{
	bool_changer (m_inhibit_updates, true);

	QAction *fa = (f == formula::tricorn ? ui->action_FormulaTricorn
		       : f == formula::ship ? ui->action_FormulaShip
		       : f == formula::celtic ? ui->action_FormulaCeltic
		       : f == formula::lambda ? ui->action_FormulaLambda
		       : f == formula::spider ? ui->action_FormulaSpider
		       : f == formula::mix ? ui->action_FormulaMix
		       : f == formula::sqtwice_a ? ui->action_FormulaSqTwiceA
		       : f == formula::sqtwice_b ? ui->action_FormulaSqTwiceB
		       : f == formula::testing ? ui->action_FormulaTest
		       : f == formula::magnet_a ? ui->action_FormulaMagnetA
		       : f == formula::facing ? ui->action_FormulaFacing
		       : f == formula::facing_b ? ui->action_FormulaFacingB
		       : ui->action_FormulaStandard);
	fa->setChecked (true);

	bool dem_allowed = formula_supports_dem (f);
	ui->menuDEM->setEnabled (dem_allowed);
	ui->action_DEMOff->setEnabled (dem_allowed);
	ui->action_DEMDist->setEnabled (dem_allowed);
	ui->action_DEMShading->setEnabled (dem_allowed);
	ui->action_DEMBoth->setEnabled (dem_allowed);
	ui->DEMGroup->setEnabled (dem_allowed);
	if (!dem_allowed)
		ui->action_DEMOff->setChecked (true);
	ui->powerSpinBox->setEnabled (f == formula::standard || f== formula::lambda || f == formula::tricorn
				      || f == formula::ship || f == formula::sqtwice_a || f == formula::sqtwice_b
				      || f == formula::celtic || f == formula::facing || f == formula::facing_b
				      || f == formula::testing);
	ui->menuHybrid->setEnabled (formula_supports_hybrid (f));
	ui->action_HybridOff->setChecked (true);

	ui->action_q_1->setEnabled (f == formula::mix);
	ui->action_q_m1->setEnabled (f == formula::mix);
	ui->action_q_2->setEnabled (f == formula::mix);
}

void MainWindow::formula_chosen (formula f, int power)
{
	if (m_inhibit_updates || (m_formula == f && power == m_power))
		return;

	bool old_inhibit_updates = m_inhibit_updates;
	m_inhibit_updates = true;
	ui->powerSpinBox->setValue (power);
	if (m_formula != f) {
		m_formula = f;
		init_formula (f);
	}
	enable_interface_for_formula (f);

	m_inhibit_updates = old_inhibit_updates;
	update_settings (true);
}

void MainWindow::update_dem_settings (QAction *)
{
	if (m_inhibit_updates)
		return;
	frac_desc &fd = current_fd ();
	bool isdem = !ui->action_DEMOff->isChecked ();
	if (isdem != fd.dem)
		update_settings (false);
	else
		update_views ();
}

void MainWindow::update_dem_color_buttons ()
{
	QPixmap p (16, 16);
	p.fill (QColor::fromRgb (m_dem_start));
	QIcon i1 (p);
	ui->DEMStartButton->setIcon (i1);
	QPixmap p2 (16, 16);
	p2.fill (QColor::fromRgb (m_dem_stop));
	QIcon i2 (p2);
	ui->DEMEndButton->setIcon (i2);
}

void MainWindow::choose_dem_color (int col)
{
	QColor old = QColor::fromRgb (col == 0 ? m_dem_start : m_dem_stop);
	QColor n = QColorDialog::getColor (old, this, tr ("Choose a color for the DEM gradient"));
	if (col == 0)
		m_dem_start = n.rgb ();
	else
		m_dem_stop = n.rgb ();
	update_dem_color_buttons ();
	update_views ();
}

void MainWindow::do_reset (bool)
{
	reset_coords (current_fd ());
	set_rotation (current_fd (), 0);
	m_reinit_render = true;
	restart_computation ();
}

void MainWindow::do_pause (bool on)
{
	m_paused = on;
	if (on)
		abort_computation ();
	else
		restart_computation ();
}

void MainWindow::do_wind_down (bool)
{
	frac_desc &fd = current_fd ();
	QMutexLocker lock (&gpu_handler->data_mutex);
	if (fd.maxiter_found > 1000)
		fd.maxiter = fd.maxiter_found;
}

void MainWindow::zoom_out (bool)
{
	int factor = ui->zoomSpinBox->value () * 10;
	frac_desc &fd = current_fd ();
	fd.bounds_w = fd.bounds_h = 0;

	fd.width = div1 (mul1 (fd.width, factor), 10);
	m_reinit_render = true;
	restart_computation ();
}

void MainWindow::zoom_in (bool)
{
	int factor = ui->zoomSpinBox->value () * 10;
	frac_desc &fd = current_fd ();
	fd.bounds_w = fd.bounds_h = 0;

	fd.width = div1 (mul1 (fd.width, 10), factor);
	m_reinit_render = true;
	restart_computation ();
}

void MainWindow::update_fractal_type (int t)
{
	ui->previewDock->setVisible (t == 2);
	ui->storePreviewButton->setEnabled (t == 2);
	update_settings (false);
}

double MainWindow::chosen_aspect ()
{
	if (!ui->aspectBox->isChecked ())
		return 0;
	int idx = ui->aspectComboBox->currentIndex ();
	return (idx == 0 ? 1
		: idx == 1 ? 4 / 3.0
		: idx == 2 ? 3 / 4.0
		: idx == 3 ? 16 / 9.0
		: idx == 4 ? 16 / 10.0
		: idx == 5 ? M_SQRT2
		: idx == 6 ? M_SQRT1_2
		: 1);
}

void MainWindow::update_aspect ()
{
       double aspect = chosen_aspect ();
       if (aspect == 0)
               ui->fractalAspectWidget->set_aspect (1, false);
       else
               ui->fractalAspectWidget->set_aspect (aspect, true);
}

void MainWindow::slot_save_as (bool)
{
	QSettings settings;
	QString ipath = settings.value ("paths/images").toString ();
	QFileDialog dlg (this, tr ("Save image file"), ipath, "PNG (*.png)");
	dlg.setAcceptMode (QFileDialog::AcceptSave);
	dlg.setDefaultSuffix (".png");
	if (!dlg.exec ()) {
		m_paused = false;
		restart_computation ();
		return;
	}

	QStringList flist = dlg.selectedFiles ();
	if (flist.isEmpty ()) {
		m_paused = false;
		restart_computation ();
		return;
	}
	QString filename = flist[0];

	QImage img = current_fd ().julia ? m_img_julia : m_img_mandel;
	if (!img.save (filename, "png"))
		QMessageBox::warning (this, PACKAGE, tr ("Failed to save image!"));
}

void MainWindow::slot_save_params ()
{
	bool_changer (m_paused, true);
	abort_computation ();

	QSettings settings;
	QString ppath = settings.value ("paths/params").toString ();
	QFileDialog dlg (this, tr ("Save parameters"), ppath, "GFF fractal params (*.fparm)");
	dlg.setAcceptMode (QFileDialog::AcceptSave);
	dlg.setDefaultSuffix (".fparm");
	if (!dlg.exec ()) {
		return;
	}

	QStringList flist = dlg.selectedFiles ();
	if (flist.isEmpty ()) {
		return;
	}
	QString filename = flist[0];
	auto fd =  current_fd ();
	QFile f (filename);
	if (!f.open (QIODevice::WriteOnly)) {
                QMessageBox::warning (this, PACKAGE, QObject::tr ("Cannot open parameter file for saving."));
                return;
	}
	QDataStream s (&f);
	s << fd;
	if (!f.flush ()) {
                QMessageBox::warning (this, PACKAGE, QObject::tr ("Error while saving to file."));
                return;
	}
	f.close ();
}

void MainWindow::slot_load_params ()
{
	bool_changer (m_paused, true);
	abort_computation ();

	QSettings settings;
	QString ppath = settings.value ("paths/params").toString ();
	QFileDialog dlg (this, tr ("Load parameters"), ppath, "GFF fractal params (*.fparm)");
	dlg.setAcceptMode (QFileDialog::AcceptOpen);
	dlg.setFileMode (QFileDialog::ExistingFile);
	dlg.setDefaultSuffix (".fparm");
	if (!dlg.exec ()) {
		return;
	}

	QStringList flist = dlg.selectedFiles ();
	if (flist.isEmpty ()) {
		return;
	}
	QString filename = flist[0];

	if (filename.isEmpty ()) {
		return;
	}

	QFile f (filename);
	f.open (QIODevice::ReadOnly);
	QDataStream s (&f);
	frac_params newfd;
	newfd.resize (max_nwords);
	s >> newfd;
	f.close ();

	restore_params (newfd);

	m_recompile = true;
}

void MainWindow::slot_save_palette ()
{
	bool_changer (m_paused, true);
	abort_computation ();

	QSettings settings;
	QString cpath = settings.value ("paths/palettes").toString ();
	QFileDialog dlg (this, tr ("Save parameters"), cpath, "GFF fractal palette (*.fpal)");
	dlg.setAcceptMode (QFileDialog::AcceptSave);
	dlg.setDefaultSuffix (".fpal");
	if (!dlg.exec ()) {
		return;
	}

	QStringList flist = dlg.selectedFiles ();
	if (flist.isEmpty ()) {
		return;
	}
	QString filename = flist[0];
	auto fd =  current_fd ();
	QFile f (filename);
	if (!f.open (QIODevice::WriteOnly)) {
                QMessageBox::warning (this, PACKAGE, QObject::tr ("Cannot open palette file for saving."));
                return;
	}
	QDataStream s (&f);
	s << m_custom_palette;
	if (!f.flush ()) {
                QMessageBox::warning (this, PACKAGE, QObject::tr ("Error while saving to file."));
                return;
	}
	f.close ();
}

void MainWindow::slot_load_palette ()
{
	bool_changer (m_paused, true);
	abort_computation ();

	QSettings settings;
	QString cpath = settings.value ("paths/palettes").toString ();
	QFileDialog dlg (this, tr ("Load parameters"), cpath, "GFF fractal palette (*.fpal)");
	dlg.setAcceptMode (QFileDialog::AcceptOpen);
	dlg.setFileMode (QFileDialog::ExistingFile);
	dlg.setDefaultSuffix (".fpal");
	if (!dlg.exec ()) {
		return;
	}

	QStringList flist = dlg.selectedFiles ();
	if (flist.isEmpty ()) {
		return;
	}
	QString filename = flist[0];

	if (filename.isEmpty ()) {
		return;
	}

	QFile f (filename);
	f.open (QIODevice::ReadOnly);
	QDataStream s (&f);
	QVector<uint32_t> vec;
	s >> m_custom_palette;
	f.close ();

	if (ui->gradComboBox->currentIndex () == ui->gradComboBox->count () - 1)
		update_palette ();
}

void MainWindow::restore_params (const frac_params &p)
{
	abort_computation ();

	{
		QMutexLocker render_lock (&m_renderer->mutex);
		QMutexLocker preview_lock (&m_preview_renderer->mutex);
		discard_fd_data (m_fd_mandel);
		discard_fd_data (m_fd_julia);
	}
	m_inhibit_updates = true;

	m_fd_mandel.param_q = p.param_q;
	m_fd_julia.param_q = p.param_q;

	m_fd_mandel.critpoint = p.critpoint;
	m_fd_julia.critpoint = p.critpoint;

	m_fd_mandel.fm = p.fm;
	m_fd_julia.fm = p.fm;

	ui->precSpinBox->setValue (p.nwords);
	ui->powerSpinBox->setValue (p.power);
	ui->typeComboBox->setCurrentIndex (p.julia ? 1 : 0);

	m_power = std::max (2, std::min (8, p.power));
	m_nwords = std::max (2, std::min (16, p.nwords));
	ui->powerSpinBox->setValue (m_power);
	ui->precSpinBox->setValue (m_nwords);
	if (p.julia) {
		m_fd_julia.set (p);
		ui->typeComboBox->setCurrentIndex (1);
	} else {
		m_fd_mandel.set (p);
		ui->typeComboBox->setCurrentIndex (0);
	}
	m_formula = p.fm;
	enable_interface_for_formula (m_formula);

	reset_coords (m_fd_mandel);
	reset_coords (m_fd_julia);

	auto &fd = p.julia ? m_fd_julia : m_fd_mandel;
	fd.center_x = p.center_x;
	fd.center_y = p.center_y;
	fd.width = p.width;
	fd.maxiter = p.maxiter;
	fd.param_p = p.param_p;
	fd.bounds_w = p.bounds_w;
	fd.bounds_h = p.bounds_h;

	m_inhibit_updates = false;

	m_recompile = true;
	restart_computation ();
}

void MainWindow::store_params (bool preview)
{
	constexpr int tn_sz = 400;
	stored_params newsp;
	newsp.fp = preview ? m_fd_julia : current_fd ();
	set_render_params (newsp.rp);
	newsp.rp.minimum = newsp.fp.julia ? m_min_julia : m_min_mandel;
	if (newsp.rp.aspect == 0)
		newsp.rp.aspect = (double)ui->fractalView->width () / ui->fractalView->height ();
	QImage img = newsp.fp.julia ? m_img_julia : m_img_mandel;
	QImage thumbnail = img.scaled (QSize (tn_sz, tn_sz), Qt::KeepAspectRatioByExpanding);
	int w = thumbnail.width ();
	int h = thumbnail.height ();
	if (w > tn_sz) {
		thumbnail = thumbnail.copy (QRect ((w - tn_sz) / 2, 0, tn_sz, tn_sz));
	} else if (h > tn_sz) {
	        thumbnail = thumbnail.copy (QRect (0, (h - tn_sz) / 2, tn_sz, tn_sz));
	}
	m_stored.emplace_back (newsp, thumbnail);
	if (m_stored.size () == 1) {
		ui->action_BatchRender->setEnabled (true);
		ui->storedDock->show ();
	}
	layout_stored_params ();
}

void MainWindow::adjust_width_for_bounds (frac_desc &fd)
{
	if (fd.bounds_w != 0 && fd.bounds_h != 0) {
		QGraphicsView *v = &fd == &current_fd () ? ui->fractalView : ui->previewView;
		double w = v->width ();
		double h = v->height ();
		double bnd;
		if (w < h) {
			bnd = fd.bounds_w;
			if (bnd * h / w < fd.bounds_h)
				bnd = fd.bounds_h / h * w;
		} else {
			bnd = fd.bounds_h;
			if (bnd * w / h < fd.bounds_w)
				bnd = fd.bounds_w / w * h;
		}
		set (fd.width, bnd);
	}
}

void MainWindow::perform_resizes ()
{
	if (m_resized_fractal) {
		adjust_width_for_bounds (current_fd ());
		m_reinit_render = true;
	}
	if (m_resized_preview) {
		adjust_width_for_bounds (m_fd_julia);
		m_preview_uptodate = false;
	}
	restart_computation ();
}

const QVector<uint32_t> &MainWindow::palette_from_index (int cb_idx)
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
		: m_custom_palette);
}

void MainWindow::update_palette ()
{
	int cb_idx = ui->gradComboBox->currentIndex ();
	auto &pal = palette_from_index (cb_idx);

	/* Switching to custom palette before any has been edited -> initialize with current.  */
	if (pal.size () == 0)
		m_custom_palette = palette_from_index (m_last_pal_idx);

	m_last_pal_idx = cb_idx;

	ui->action_SavePalette->setEnabled (m_custom_palette.size () > 0);
	bool narrowb = ui->action_NarrowB->isChecked ();
	bool narroww = ui->action_NarrowW->isChecked ();
	int nfactor = ui->action_NFactor4->isChecked () ? 4 : 2;
	int hue_shift = ui->hueSlider->value ();
	bool trad = ui->action_RotateTrad->isChecked ();
	m_palette = interpolate_colors (pal, interpolation_factor, hue_shift, trad, narrowb, narroww, nfactor);
	if (ui->structureGroup->isChecked ()) {
		int strtype = ui->action_StructLight->isChecked () ? 0 : ui->action_StructDark->isChecked () ? 1 : 2;
		int str = ui->structureSlider->value ();
		for (int i = 0; i < m_palette.size (); i++) {
			double y2 = i;
			y2 /= interpolation_factor;
			y2 *= str;
			double m1b = y2 - floor (y2);

			m1b = (m1b - 0.5) * 2;
			double m1bs = 1 - m1b * m1b * m1b * m1b / 4;
			int primary = m_palette[i];
			int r = (primary >> 16) & 255;
			int g = (primary >> 8) & 255;
			int b = (primary & 255);
			if ((strtype == 2 && m1b < 0) || strtype == 1) {
				r = r * m1bs;
				g = g * m1bs;
				b = b * m1bs;
			} else {
				r = 255 - (255 - r) * m1bs;
				g = 255 - (255 - g) * m1bs;
				b = 255 - (255 - b) * m1bs;
			}
			m_palette[i] = (r << 16) + (g << 8) + b;
		}
	}
	update_views ();
}

void MainWindow::gradient_edit (bool)
{
	bool old_paused = m_paused;
	m_paused = true;
	abort_computation ();
	QVector<uint32_t> saved_pal = m_palette;
	QVector<uint32_t> new_pal = m_palette;

	const QVector<uint32_t> &pal = palette_from_index (ui->gradComboBox->currentIndex ());
	bool narrowb = ui->action_NarrowB->isChecked ();
	bool narroww = ui->action_NarrowW->isChecked ();
	int nfactor = ui->action_NFactor4->isChecked () ? 4 : 2;
	int hue_shift = ui->hueSlider->value ();
	bool trad = ui->action_RotateTrad->isChecked ();
	GradEditor edit (this, pal);
	connect (&edit, &GradEditor::colors_changed,
		 [this, &new_pal, narrowb, hue_shift, trad, narroww, nfactor] (const QVector<uint32_t> &newcols)
		 {
			 new_pal = newcols;
			 m_palette = interpolate_colors (newcols, interpolation_factor, hue_shift, trad,
							 narrowb, narroww, nfactor);
			 update_views ();
		 });
	if (edit.exec ()) {
		m_custom_palette = new_pal;
		ui->action_SavePalette->setEnabled (m_custom_palette.size () > 0);
		ui->gradComboBox->setCurrentIndex (ui->gradComboBox->count () - 1);
	} else {
		m_palette = saved_pal;
		update_views ();
	}
	m_paused = old_paused;
	restart_computation ();
}

void MainWindow::slot_batchrender (bool)
{
	bool old_paused = m_paused;
	m_paused = true;
	abort_computation ();

	BatchRenderDialog dlg (this);
	if (!dlg.exec ())
		return;

	QSettings settings;

	QProgressDialog pdlg (tr ("Rendering images..."), tr ("Abort"), 0, 100, this);
	pdlg.setWindowModality (Qt::WindowModal);
	pdlg.setMinimumDuration (0);
	pdlg.show ();

	QString pattern = dlg.get_file_template ();
	int count = 1;
	int samples = dlg.get_samples ();
	int prev_maxiter_factor = dlg.get_prev_maxiter ();
	bool preserve = dlg.get_preserve_aspect ();
	double progress = 0;
	double pro_step = 100. / m_stored.size ();
	printf ("pro_step %f\n", pro_step);
	for (auto st: m_stored) {
		QString v_str = QString::number (count++);
		while (v_str.length () < 4)
			v_str = "0" + v_str;
		QString filename = pattern;
		filename.replace (QRegExp ("%n"), v_str);
		QString label = tr ("Working on file: ") + filename;
		pdlg.setLabelText (label);

		frac_params &fp = st.params.fp;
		render_params &rp = st.params.rp;
		int w = dlg.get_width ();
		int h = dlg.get_height ();

		double dlg_aspect = (double)w / h;
		if (dlg_aspect != rp.aspect && (dlg_aspect > 1) != (rp.aspect > 1)) {
			std::swap (w, h);
			dlg_aspect = 1 / dlg_aspect;
		}
		if (preserve) {
			if (dlg_aspect > rp.aspect)
				w = rp.aspect * h;
			else
				h = w / rp.aspect;
		}

		QString errstr;
		auto it = std::find (std::begin (formula_table), std::end (formula_table), fp.fm);
		int fidx = it - std::begin (formula_table);
		int n_prev = 1;
		if (rp.sac) {
			n_prev = 1 << settings.value ("coloring/nprev").toInt ();
		}
		emit signal_compile_kernel (fidx, fp.power, fp.nwords, max_nwords, &errstr);
		gpu_handler->done_sem.acquire ();
		if (!errstr.isEmpty ()) {
			QMessageBox::critical (this, PACKAGE, errstr);
			close ();
			return;
		}
		QFile f (filename);
		if (f.exists () && !dlg.get_overwrite ()) {
			QMessageBox::StandardButton choice;
			choice = QMessageBox::warning (this, tr ("File exists"),
						       tr ("A filename matching the pattern and current number already exists.  Overwrite?"),
						       QMessageBox::Yes | QMessageBox::No, QMessageBox::Yes);
			if (choice == QMessageBox::No)
				break;
		}

		frac_desc temp_fd;
		temp_fd.set (fp);
		temp_fd.generation = 0;
		int spixels = 1 << samples;
		double d = (double)w * spixels * spixels * n_prev;
		int batch_size = std::max (20.0, 1000000.0 / d);
		int steps = (h + batch_size - 1) / batch_size;
		batch_size = h / steps;

		Renderer renderer;
		renderer.queued = true;
		renderer.next_rp = rp;
		renderer.render_width = w;
		bool alpha = settings.value ("render/alpha").toBool ();
		renderer.result_image = QImage (w, h, alpha ? QImage::Format_ARGB32 : QImage::Format_RGB32);
		int maxiter = default_maxiter;
		int dlg_maxiter = dlg.get_maxiter ();
		int prev_maxiter = fp.maxiter_found;
		if (dlg_maxiter != 0)
			maxiter = std::min (maxiter, default_maxiter);
		if (prev_maxiter_factor != 0)
			maxiter = std::min (maxiter, (int)(prev_maxiter / 100.0 * prev_maxiter_factor));
		printf ("maxiter is %d (%d %d %d)\n", maxiter, default_maxiter, dlg_maxiter, prev_maxiter);
		for (int y0 = 0; y0 < h; y0 += batch_size) {
			pdlg.setValue (progress + pro_step * ((double)y0 / h));

			temp_fd.yoff = y0 * (1 << samples);
			int this_h = std::min (h - y0, batch_size);
			compute_fractal (temp_fd, temp_fd.nwords, n_prev, w, this_h, h, maxiter,
					 samples, rp.dem || rp.dem_shade, false, true);
			gpu_handler->done_sem.acquire ();
			if (pdlg.wasCanceled ())
				break;
			renderer.render_height = this_h;
			renderer.set_minimum (rp.minimum, temp_fd.generation);
			renderer.do_render (rp, w, this_h, y0, &temp_fd, nullptr, temp_fd.generation);
			if (pdlg.wasCanceled ())
				break;
		}

		discard_fd_data (temp_fd);

		if (pdlg.wasCanceled ())
			break;

		if (!renderer.result_image.save (&f)) {
			QMessageBox::warning (this, tr ("Error while saving"),
					      tr ("One of the files could not be saved.\nPlease verify the filename pattern is correct."));
			break;
		}
		f.close ();
		progress += pro_step;
	}

	m_paused = old_paused;
	m_recompile = true;
	restart_computation ();
}

void MainWindow::help_about ()
{
	QString txt = "<p>" PACKAGE "</p>";
	txt = "<p>Copyright \u00a9 2021\nBernd Schmidt &lt;bernds_cb1@t-online.de&gt;</p>";
	txt += "<p>This is a fractal image generator using arbitrary precision arithmetic on a GPU. It is still experimental. It is distributed in the hope that you will find it fun to use already, but expect it to be rough around the edges.</p>";

	QMessageBox mb (this);
	mb.setWindowTitle (PACKAGE);
	mb.setTextFormat (Qt::RichText);
	mb.setText (txt);
	mb.exec ();
}

void MainWindow::choose_hybrid (bool on)
{
	frac_desc &fd = current_fd ();
	if (!on) {
		ui->action_HybridOn->setText (tr ("On"));
		fd.hybrid_code = 0;
		fd.hybrid_len = 0;
		update_settings (true);
		return;
	}
	if (m_inhibit_updates)
		return;

	HybridDialog dlg (this);
	if (!dlg.exec ()) {
		ui->action_HybridOff->setChecked (true);
	}
	QString code = dlg.get_code ();
	int len = code.length ();
	if (len < 2 || len > 32) {
		QMessageBox::warning (this, tr ("Invalid hybrid code"),
				      tr ("The code should have a length between 2 and 32."));
		ui->action_HybridOn->setText (tr ("On"));
		ui->action_HybridOn->setText ("On");
		ui->action_HybridOff->setChecked (true);
		return;
	}
	int c0 = code.count ('0');
	int c1 = code.count ('1');
	if (c0 + c1 != len) {
		QMessageBox::warning (this, tr ("Invalid hybrid code"),
				      tr ("The code should contain a pattern made of only 0s and 1s."));
		ui->action_HybridOn->setText (tr ("On"));
		ui->action_HybridOff->setChecked (true);
		return;
	}
	uint32_t bin_code = 0;
	for (int i = 0; i < len; i++) {
		bin_code *= 2;
		if (code[i] == '1')
			bin_code++;
	}
	ui->action_HybridOn->setText (tr ("On: ") + code);
	m_fd_mandel.hybrid_code = bin_code;
	m_fd_mandel.hybrid_len = len;
	m_fd_julia.hybrid_code = bin_code;
	m_fd_julia.hybrid_len = len;
	update_settings (true);
}

MainWindow::MainWindow ()
	: ui (new Ui::MainWindow)
{
	ui->setupUi (this);
	ui->fractalView->setScene (&m_canvas);
	ui->previewView->setScene (&m_preview_canvas);
	ui->storedView->setScene (&m_stored_canvas);

	restore_geometry ();

	start_threads ();

	init_formula (formula::standard);
	m_fd_mandel.fm = m_fd_julia.fm = m_formula;
	m_fd_julia.julia = true;

	m_fd_mandel.resize (max_nwords);
	m_fd_julia.resize (max_nwords);

	set_rotation (m_fd_mandel, 0);
	set_rotation (m_fd_julia, 0);

#ifndef TESTING
	ui->action_FormulaTest->setVisible (false);
#endif

	QString errstr;
	m_power = default_power;
	auto it = std::find (std::begin (formula_table), std::end (formula_table), m_formula);
	int fidx = it - std::begin (formula_table);
	emit signal_compile_kernel (fidx, default_power, m_nwords, max_nwords, &errstr);
	gpu_handler->done_sem.acquire ();
	if (!errstr.isEmpty ()) {
		QMessageBox::critical (this, PACKAGE, errstr);
		close ();
		return;
	}

	ui->precSpinBox->setValue (m_nwords);
	ui->precSpinBox->setMinimum (2);
	ui->precSpinBox->setMaximum (max_nwords);
	ui->powerSpinBox->setValue (default_power);
	ui->powerSpinBox->setMinimum (2);
	ui->powerSpinBox->setMaximum (7);
	ui->widthSpinBox->setValue (1);
	ui->widthSpinBox->setMinimum (1);
	ui->widthSpinBox->setMaximum (20);

	ui->zoomSpinBox->setMinimum (1.1);
	ui->zoomSpinBox->setMaximum (5);
	ui->sampleSpinBox->setValue (1);
	ui->action_NarrowB->setChecked (true);
	ui->action_NarrowW->setChecked (true);
	ui->action_RotateHSV->setChecked (true);
	ui->action_DEMOff->setChecked (true);
	ui->action_DEMColour->setChecked (false);
	ui->action_AngleColour->setChecked (false);
	ui->action_AngleNone->setChecked (true);
	ui->action_AimAssist->setChecked (true);

	reset_coords (m_fd_mandel);
	reset_coords (m_fd_julia);

	ui->storePreviewButton->setEnabled (ui->typeComboBox->currentIndex () == 2);
	ui->previewDock->setVisible (ui->typeComboBox->currentIndex () == 2);
	ui->storedDock->hide ();
	ui->extraDock->hide ();
	ui->action_BatchRender->setEnabled (false);
	ui->menuHybrid->setEnabled (false);
	ui->stripesWidget->setEnabled (false);
	render_fractal ();

	m_resize_timer.setSingleShot (true);
	connect (&m_resize_timer, &QTimer::timeout, this, &MainWindow::perform_resizes);

	connect (ui->fractalView, &SizeGraphicsView::resized, [this] () { m_resized_fractal = true; m_resize_timer.start (500); });
	connect (ui->previewView, &SizeGraphicsView::resized, [this] () { m_resized_preview = true; m_resize_timer.start (500); });
	connect (ui->storedView, &SizeGraphicsView::resized, this, &MainWindow::layout_stored_params);
	connect (ui->storedDock->toggleViewAction (), &QAction::toggled,
		 [this] (bool on) { ui->storedDock->setVisible (on); if (on) layout_stored_params (); });

	ui->previewAspectWidget->set_child (ui->previewView);
	ui->fractalAspectWidget->set_child (ui->fractalView);

	update_aspect ();

	constexpr int default_palette = 3;
	m_last_pal_idx = default_palette;
	ui->gradComboBox->setCurrentIndex (default_palette);
	update_palette ();

	m_struct_group = new QActionGroup (this);
	m_struct_group->addAction (ui->action_StructLight);
	m_struct_group->addAction (ui->action_StructDark);
	m_struct_group->addAction (ui->action_StructBoth);

	m_power_group = new QActionGroup (this);
	m_power_group->addAction (ui->action_FD2);
	m_power_group->addAction (ui->action_FD3);
	m_power_group->addAction (ui->action_FD4);
	m_power_group->addAction (ui->action_FD5);
	m_power_group->addAction (ui->action_FD6);
	m_power_group->addAction (ui->action_FD7);

	m_narrow_group = new QActionGroup (this);
	m_narrow_group->addAction (ui->action_NFactor2);
	m_narrow_group->addAction (ui->action_NFactor4);

	m_rotate_group = new QActionGroup (this);
	m_rotate_group->addAction (ui->action_RotateTrad);
	m_rotate_group->addAction (ui->action_RotateHSV);

	// The idea behind these is that if we shift to 0 and apply something like
	// cbrt afterwards, we end up with much high colour variation on the outside.
	// Shifting to a higher minimum reduces that effect.
	m_sub_group = new QActionGroup (this);
	m_sub_group->addAction (ui->action_ShiftNone);
	m_sub_group->addAction (ui->action_Shift0);
	m_sub_group->addAction (ui->action_Shift10);
	m_sub_group->addAction (ui->action_Shift100);

	m_dem_group = new QActionGroup (this);
	m_dem_group->addAction (ui->action_DEMOff);
	m_dem_group->addAction (ui->action_DEMDist);
	m_dem_group->addAction (ui->action_DEMShading);
	m_dem_group->addAction (ui->action_DEMBoth);

	m_hybrid_group = new QActionGroup (this);
	m_hybrid_group->addAction (ui->action_HybridOff);
	m_hybrid_group->addAction (ui->action_HybridOn);

	m_incolor_group = new QActionGroup (this);
	m_incolor_group->addAction (ui->action_ICBlack);
	m_incolor_group->addAction (ui->action_ICWhite);

	m_angles_group = new QActionGroup (this);
	m_angles_group->addAction (ui->action_AngleNone);
	m_angles_group->addAction (ui->action_AngleSmooth);
	m_angles_group->addAction (ui->action_AngleBin);
	m_angles_group->addAction (ui->action_SAC);
	m_angles_group->addAction (ui->action_TIA);

	m_formula_group = new QActionGroup (this);
	m_formula_group->addAction (ui->action_FormulaStandard);
	m_formula_group->addAction (ui->action_FormulaTricorn);
	m_formula_group->addAction (ui->action_FormulaShip);
	m_formula_group->addAction (ui->action_FormulaCeltic);
	m_formula_group->addAction (ui->action_FormulaLambda);
	m_formula_group->addAction (ui->action_FormulaSpider);
	m_formula_group->addAction (ui->action_FormulaMix);
	m_formula_group->addAction (ui->action_FormulaSqTwiceA);
	m_formula_group->addAction (ui->action_FormulaSqTwiceB);
	m_formula_group->addAction (ui->action_FormulaMagnetA);
	m_formula_group->addAction (ui->action_FormulaFacing);
	m_formula_group->addAction (ui->action_FormulaFacingB);
	m_formula_group->addAction (ui->action_FormulaTest);

	ui->action_Shift10->setChecked (true);
	ui->action_NFactor4->setChecked (true);
	ui->action_StructDark->setChecked (true);
	ui->action_ICBlack->setChecked (true);
	ui->action_Contrast->setChecked (true);
	ui->action_FD2->setChecked (true);
	enable_interface_for_formula (m_formula);
	enable_interface_for_settings ();

	ui->menu_View->insertAction (nullptr, ui->storedDock->toggleViewAction ());
	ui->menu_View->insertAction (nullptr, ui->extraDock->toggleViewAction ());
	ui->extraDock->toggleViewAction ()->setShortcut (Qt::CTRL + Qt::Key_X);

	connect (ui->fractalView, &SizeGraphicsView::mouse_event, this, &MainWindow::fractal_mouse_event);
	connect (ui->fractalView, &SizeGraphicsView::wheel_event, this, &MainWindow::fractal_wheel_event);
	connect (ui->previewView, &SizeGraphicsView::wheel_event, this, &MainWindow::preview_wheel_event);

	void (QComboBox::*cic) (int) = &QComboBox::currentIndexChanged;
	void (QSpinBox::*changed) (int) = &QSpinBox::valueChanged;
	void (QDoubleSpinBox::*dchanged) (double) = &QDoubleSpinBox::valueChanged;
	connect (ui->precSpinBox, changed,
		 // The spinbox value can be changed programmatically, so we check before calling update_settings.
		 [this] (int v) { if (v != m_nwords) update_settings (false); });
	connect (ui->powerSpinBox, changed, [this] (int) { update_settings (true); });
	connect (ui->sampleSpinBox, changed, [this] (int) { update_settings (false); });
	connect (ui->typeComboBox, cic, this, &MainWindow::update_fractal_type);
	connect (ui->widthSpinBox, changed, this, &MainWindow::update_views);
	connect (ui->stripesSpinBox, dchanged, this, &MainWindow::update_views);
	connect (ui->modifyComboBox, cic, this, &MainWindow::update_views);
	connect (ui->gradComboBox, cic,  [this] (int) { update_palette (); });
	connect (ui->colStepSlider, &QSlider::valueChanged, this, &MainWindow::update_views);
	connect (m_sub_group, &QActionGroup::triggered, [this] (QAction *) { update_views (); });
	connect (m_dem_group, &QActionGroup::triggered, this, &MainWindow::update_dem_settings);

	update_dem_color_buttons ();
	connect (ui->DEMStartButton, &QToolButton::clicked, [this] (bool) { choose_dem_color (0); });
	connect (ui->DEMEndButton, &QToolButton::clicked, [this] (bool) { choose_dem_color (1); });

	connect (ui->structureGroup, &QGroupBox::toggled, [this] (bool) { update_palette (); });
	connect (ui->structureSlider, &QSlider::valueChanged, [this] (bool) { update_palette (); });
	connect (ui->hueSlider, &QSlider::valueChanged, [this] (bool) { update_palette (); });
	connect (ui->action_NarrowB, &QAction::toggled, [this] (bool) { update_palette (); });
	connect (ui->action_NarrowW, &QAction::toggled, [this] (bool) { update_palette (); });
	connect (ui->action_NFactor2, &QAction::toggled, [this] (bool) { update_palette (); });
	connect (ui->action_NFactor4, &QAction::toggled, [this] (bool) { update_palette (); });
	connect (ui->action_RotateTrad, &QAction::toggled, [this] (bool) { update_palette (); });
	connect (ui->action_RotateHSV, &QAction::toggled, [this] (bool) { update_palette (); });
	connect (ui->action_StructLight, &QAction::toggled, [this] (bool) { update_palette (); });
	connect (ui->action_StructDark, &QAction::toggled, [this] (bool) { update_palette (); });
	connect (ui->action_StructBoth, &QAction::toggled, [this] (bool) { update_palette (); });
	connect (ui->action_AngleNone, &QAction::toggled, this, &MainWindow::slot_disable_sac);
	connect (ui->action_AngleSmooth, &QAction::toggled, this, &MainWindow::slot_disable_sac);
	connect (ui->action_AngleBin, &QAction::toggled, this, &MainWindow::slot_disable_sac);
	connect (ui->action_SAC, &QAction::toggled, [this] (bool) { enable_sac_or_tia (); });
	connect (ui->action_TIA, &QAction::toggled, [this] (bool) { enable_sac_or_tia (); });

	connect (ui->action_ICBlack, &QAction::toggled, [this] (bool) { update_views (); });
	connect (ui->action_ICWhite, &QAction::toggled, [this] (bool) { update_views (); });

	connect (ui->action_DEMColour, &QAction::toggled,
		 [this] (bool) { if (!ui->action_DEMOff->isChecked ()) update_views (); });
	connect (ui->action_AngleColour, &QAction::toggled,
		 [this] (bool) { if (!ui->action_AngleNone->isChecked ()) update_views (); });

	connect (ui->action_Reset, &QAction::triggered, this, &MainWindow::do_reset);

	connect (ui->pauseButton, &QPushButton::toggled, this, &MainWindow::do_pause);
	connect (ui->windDownButton, &QPushButton::clicked, this, &MainWindow::do_wind_down);
	connect (ui->zinButton, &QPushButton::clicked, this, &MainWindow::zoom_in);
	connect (ui->zoutButton, &QPushButton::clicked, this, &MainWindow::zoom_out);
	connect (ui->storeButton, &QPushButton::clicked, [this] (bool) { store_params (false); });
	connect (ui->storePreviewButton, &QPushButton::clicked, [this] (bool) { store_params (true); });

	connect (ui->demParamSpinBox, dchanged, [this] (bool) { update_views (); });
	connect (ui->demStrengthSpinBox, dchanged, [this] (bool) { update_views (); });
	connect (ui->aspectBox, &QGroupBox::toggled, [this] (bool) { update_aspect (); });
	connect (ui->aspectComboBox, cic, [this] (int) { update_aspect (); });

	connect (ui->action_Mandelbrot, &QAction::triggered,
		 [this] (bool) { ui->typeComboBox->setCurrentIndex (0); });
	connect (ui->action_Julia, &QAction::triggered,
		 [this] (bool) { ui->typeComboBox->setCurrentIndex (1); });
	connect (ui->action_MJpreview, &QAction::triggered,
		 [this] (bool) { ui->typeComboBox->setCurrentIndex (2); });

	connect (ui->action_IncPrec, &QAction::triggered,
		 [this] (bool) { ui->precSpinBox->stepBy (1); });
	connect (ui->action_DecPrec, &QAction::triggered,
		 [this] (bool) { ui->precSpinBox->stepBy (-1); });

	connect (ui->action_Rotate0, &QAction::triggered,
		 [this] (bool) { set_rotation (current_fd (), 0); update_settings (false); });
	connect (ui->action_RotateEnter, &QAction::triggered, this, &MainWindow::enter_rotation);
	connect (ui->action_Rotate5, &QAction::triggered,
		 [this] (bool) { inc_rotation (current_fd (), 5); update_settings (false); });
	connect (ui->action_Rotate30, &QAction::triggered,
		 [this] (bool) { inc_rotation (current_fd (), 30); update_settings (false); });
	connect (ui->action_Rotate45, &QAction::triggered,
		 [this] (bool) { inc_rotation (current_fd (), 45); update_settings (false); });
	connect (ui->action_Rotate90, &QAction::triggered,
		 [this] (bool) { inc_rotation (current_fd (), 90); update_settings (false); });
	connect (ui->shearComboBox, cic,
		 [this] (bool) { inc_rotation (current_fd (), 0); update_settings (false); });
	connect (ui->shearSlider, &QSlider::valueChanged,
		 [this] (bool) { inc_rotation (current_fd (), 0); update_settings (false); });
	connect (ui->scaleComboBox, cic,
		 [this] (bool) { inc_rotation (current_fd (), 0); update_settings (false); });
	connect (ui->scaleSlider, &QSlider::valueChanged,
		 [this] (bool) { inc_rotation (current_fd (), 0); update_settings (false); });
	connect (ui->action_q_1, &QAction::triggered, [this] (bool) { set_q (1, 0); update_settings (true); });
	connect (ui->action_q_m1, &QAction::triggered, [this] (bool) { set_q (-1, 0); update_settings (true); });
	connect (ui->action_q_2, &QAction::triggered, [this] (bool) { set_q (2, 0); update_settings (true); });

	connect (ui->action_FD2, &QAction::triggered, [this] (bool) { formula_chosen (m_formula, 2); });
	connect (ui->action_FD3, &QAction::triggered, [this] (bool) { formula_chosen (m_formula, 3); });
	connect (ui->action_FD4, &QAction::triggered, [this] (bool) { formula_chosen (m_formula, 4); });
	connect (ui->action_FD5, &QAction::triggered, [this] (bool) { formula_chosen (m_formula, 5); });
	connect (ui->action_FD6, &QAction::triggered, [this] (bool) { formula_chosen (m_formula, 6); });
	connect (ui->action_FD7, &QAction::triggered, [this] (bool) { formula_chosen (m_formula, 6); });
	connect (ui->action_FormulaStandard, &QAction::triggered,
		 [this] (bool) { formula_chosen (formula::standard, 2); });
	connect (ui->action_FormulaLambda, &QAction::triggered,
		 [this] (bool) { formula_chosen (formula::lambda, 3); });
	connect (ui->action_FormulaSpider, &QAction::triggered,
		 [this] (bool) { formula_chosen (formula::spider, 2); });
	connect (ui->action_FormulaTricorn, &QAction::triggered,
		 [this] (bool) { formula_chosen (formula::tricorn, 2); });
	connect (ui->action_FormulaShip, &QAction::triggered,
		 [this] (bool) { formula_chosen (formula::ship, 2); });
	connect (ui->action_FormulaCeltic, &QAction::triggered,
		 [this] (bool) { formula_chosen (formula::celtic, 2); });
	connect (ui->action_FormulaMix, &QAction::triggered,
		 [this] (bool) { formula_chosen (formula::mix, 3); });
	connect (ui->action_FormulaSqTwiceA, &QAction::triggered,
		 [this] (bool) { formula_chosen (formula::sqtwice_a, 2); });
	connect (ui->action_FormulaSqTwiceB, &QAction::triggered,
		 [this] (bool) { formula_chosen (formula::sqtwice_b, 2); });
	connect (ui->action_FormulaMagnetA, &QAction::triggered,
		 [this] (bool) { formula_chosen (formula::magnet_a, 2); });
	connect (ui->action_FormulaFacing, &QAction::triggered,
		 [this] (bool) { formula_chosen (formula::facing, 2); });
	connect (ui->action_FormulaFacingB, &QAction::triggered,
		 [this] (bool) { formula_chosen (formula::facing_b, 2); });
	connect (ui->action_FormulaTest, &QAction::triggered,
		 [this] (bool) { formula_chosen (formula::testing, 2); });

	connect (ui->action_HybridOn, &QAction::triggered, [this] (bool) { choose_hybrid (true); });
	connect (ui->action_HybridOff, &QAction::triggered, [this] (bool) { choose_hybrid (false); });

	ui->action_SavePalette->setEnabled (m_custom_palette.size () > 0);
	connect (ui->action_SaveImageAs, &QAction::triggered, this, &MainWindow::slot_save_as);
	connect (ui->action_SaveParams, &QAction::triggered, [this] (bool) { slot_save_params (); restart_computation (); });
	connect (ui->action_LoadParams, &QAction::triggered, [this] (bool) { slot_load_params (); restart_computation (); });
	connect (ui->action_SavePalette, &QAction::triggered, [this] (bool) { slot_save_palette (); restart_computation (); });
	connect (ui->action_LoadPalette, &QAction::triggered, [this] (bool) { slot_load_palette (); restart_computation (); });
	connect (ui->action_GradEditor, &QAction::triggered, this, &MainWindow::gradient_edit);
	connect (ui->action_BatchRender, &QAction::triggered, this, &MainWindow::slot_batchrender);

	connect (ui->action_Prefs, &QAction::triggered,
		 [this] (bool) {
			 PrefsDialog dlg (this);
			 if (dlg.exec ()) {
				 abort_computation ();
				 m_recompile = true;
				 m_inhibit_updates = true;
				 enable_interface_for_settings ();
				 enable_sac_or_tia ();
				 m_inhibit_updates = false;
				 restart_computation ();
			 }
		 });
	connect (ui->action_About, &QAction::triggered, [=] (bool) { help_about (); });
	connect (ui->action_AboutQt, &QAction::triggered, [=] (bool) { QMessageBox::aboutQt (this); });

	addActions ({ ui->action_Mandelbrot, ui->action_Julia, ui->action_MJpreview });
	addActions ({ ui->action_IncPrec, ui->action_DecPrec });
	addActions ({ ui->action_FD2, ui->action_FD3, ui->action_FD4, ui->action_FD5 });
	addActions ({ ui->action_FD6, ui->action_FD7  });
	addActions ({ ui->action_AngleSmooth  });
	addActions ({ ui->action_SaveImageAs });
	addActions ({ ui->action_GradEditor });
}

MainWindow::~MainWindow ()
{
	delete m_narrow_group;
	delete m_struct_group;
	delete m_formula_group;
	delete m_power_group;
	delete m_rotate_group;
	delete m_sub_group;
	delete m_dem_group;
	delete m_angles_group;
	delete m_hybrid_group;
	delete ui;
}

int main (int argc, char **argv)
{
	QApplication::setAttribute (Qt::AA_EnableHighDpiScaling);
	QApplication myapp (argc, argv);

	myapp.setOrganizationName ("bernds");
	myapp.setApplicationName (PACKAGE);
#if 0
	QCommandLineParser cmdp;

	cmdp.addHelpOption ();

	cmdp.process (myapp);
#endif

	QSettings settings;
	if (!settings.contains ("coloring/nosuper-sac"))
		settings.setValue ("coloring/nosuper-sac", true);
	bool shown = settings.contains ("helpshown");
	if (!shown) {
		QMessageBox::information (nullptr, PACKAGE,
					  QObject::tr ("<p>Welcome to " PACKAGE "!</p>")
					  + QObject::tr ("<p>This program is stll in development, but hopefully already fun to use.</p><p>Left-click to zoom. Use Ctrl-click to select a parameter for Julia sets. Use shortcuts M, J and P to switch between Mandelbrot, Julia or Mandelbrot with preview modes.</p><p>This dialog will not be shown again.</p>"));
		settings.setValue ("helpshown", true);
	}

	auto w = new MainWindow ();
	w->show ();
	auto retval = myapp.exec ();

	return retval;
}
