#include <QMainWindow>
#include <QMutex>
#include <QGraphicsScene>
#include <QImage>
#include <QTimer>
#include <QThreadPool>

#include "fpvec.h"
#include "fractal.h"

// RAII wrapper around temporarily setting m_inhibit_updates in MainWindow
class bool_changer
{
	bool &m_var;
	bool m_old;

public:
	bool_changer (bool &var, bool set) : m_var (var), m_old (var)
	{
		m_var = set;
	}
	~bool_changer ()
	{
		m_var = m_old;
	}
};

namespace Ui
{
	class MainWindow;
};

class ClickablePixmap;
class QActionGroup;

struct render_params
{
	QVector<uint32_t> palette;
	int incol;
	int mod_type;
	double steps;
	int slider;
	int sub_val;
	int angle;
	bool sub;
	bool dem_colour, angle_colour;
	bool dem;
	bool dem_shade;
	uint32_t dem_start, dem_stop;
	double dem_param, dem_strength;
	// Used for stored parameters, holds either the aspect set in the GUI,
	// or, if that is disabled, the image dimensions.
	double aspect;
	// Minimum niter value, stored for batch rendering (which is done
	// piecewise and therefor cannot compute it)
	double minimum;
};

struct stored_params
{
	frac_params fp;
	render_params rp;
};

struct stored_preview {
	ClickablePixmap *pixmap {};
	stored_params params;
	QImage thumbnail;
	int x = 0, y = 0;

	stored_preview (const stored_params &p, QImage tn)
		: params (p), thumbnail (tn)
	{
	}
	stored_preview (const stored_preview &) = default;
	stored_preview (stored_preview &&) = default;
	stored_preview &operator= (const stored_preview &) = default;
	stored_preview &operator= (stored_preview &&) = default;
};

class Renderer : public QObject
{
	Q_OBJECT
	QThreadPool m_pool;
	double m_minimum = 0;
	int m_min_gen = -1;

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

class MainWindow: public QMainWindow
{
	Q_OBJECT

	Ui::MainWindow *ui;

	QThread *m_gpu_thread {}, *m_render_thread {}, *m_preview_thread {};

	QTimer m_resize_timer;
	bool m_resized_fractal = false, m_resized_preview = false;

	static constexpr int max_nwords = 16;
	static constexpr int iter_steps = 2000;
	static constexpr int default_maxiter = 500000;

	formula m_formula = formula::standard;

	frac_desc m_fd_mandel;
	frac_desc m_fd_julia;

	/* Color gradient for DEM.  */
	uint32_t m_dem_start = 0xFFFFFF;
	uint32_t m_dem_stop = 0;
	QImage m_img_mandel, m_img_julia;
	// Minimum niter values, computed during rendering
	double m_min_mandel = 0, m_min_julia = 0;

	vector<stored_preview> m_stored;

	int m_generation = 0;
	int m_nwords = 2;
	int m_power;

	bool m_new_data_queued = false;
	bool m_preview_uptodate = false;

	bool m_reinit_render = false;
	bool m_recompile = false;

	bool m_working = false;
	bool m_inhibit_updates = false;
	bool m_paused = false;

	Renderer *m_renderer;
	Renderer *m_preview_renderer;

	QGraphicsScene m_canvas;
	QGraphicsScene m_preview_canvas;
	QGraphicsScene m_stored_canvas;

	QActionGroup *m_formula_group {};
	QActionGroup *m_sub_group {};
	QActionGroup *m_dem_group {};
	QActionGroup *m_hybrid_group {};
	QActionGroup *m_power_group {};
	QActionGroup *m_struct_group {};
	QActionGroup *m_narrow_group {};
	QActionGroup *m_rotate_group {};
	QActionGroup *m_angles_group {};

	int m_last_pal_idx = 0;
	QVector<uint32_t> m_palette;
	QVector<uint32_t> m_custom_palette;

	void restore_geometry ();
	void start_threads ();
	void do_compile ();

	void perform_resizes ();

	void inc_rotation (frac_desc &, int);
	void set_rotation (frac_desc &, double);
	void enter_rotation (bool);
	double shear_slider_value ();
	double scale_slider_value ();
	void build_points (frac_desc &, int w, int h);
	void compute_fractal (frac_desc &, int nwords, int w, int h, int full_h,
			      int maxiter, int ss, bool dem, bool preview, bool batch = false);
	void render_fractal ();
	void render_preview ();
	bool abort_computation ();
	void restart_computation ();

	void autoprec (frac_desc &);
	void update_settings (bool);
	void update_dem_settings (QAction *);

	// double iter_value_at (frac_desc &, int);
	void update_display (QGraphicsView *);
	void precompute_iter_value (frac_desc *);
	double chosen_aspect ();
	void update_aspect ();
	void update_views (int = 0);
	const QVector<uint32_t> &palette_from_index (int);
	void update_palette ();
	void update_fractal_type (int = 0);
	void adjust_width_for_bounds (frac_desc &);
	void reset_coords (frac_desc &);
	void set_q (int, int);
	frac_desc &current_fd ();
	void discard_fd_data (frac_desc &);

	void do_pause (bool);
	void do_wind_down (bool);
	void do_reset (bool);
	void zoom_in (bool = false);
	void zoom_out (bool = false);

	void choose_dem_color (int);
	void update_dem_color_buttons ();

	void set_render_params (render_params &);
        void slot_new_data (frac_desc *, int, bool);
	void slot_kernel_complete ();
	void slot_render_complete (QGraphicsView *view, frac_desc *, QImage, double);

	void layout_stored_params ();
	void store_params (bool);
	void restore_params (const frac_params &);

	void slot_save_as (bool);
	void slot_save_params ();
	void slot_load_params ();
	void slot_save_palette ();
	void slot_load_palette ();
	void slot_batchrender (bool);

	void gradient_edit (bool);

	void fractal_mouse_event (QMouseEvent *);
	void fractal_wheel_event (QWheelEvent *);
	void preview_wheel_event (QWheelEvent *);

	void enable_interface_for_formula (formula);
	void init_formula (formula);
	void formula_chosen (formula, int);
	void choose_hybrid (bool);

	void help_about ();
protected:
	void closeEvent(QCloseEvent *event) override;

public:
	MainWindow ();
	~MainWindow ();

signals:
	void signal_init_cuda (QString *err);
	void signal_start_kernel (frac_desc *, int generation, int, int, bool);
	void signal_alloc_mem (frac_desc *, int, int, int, int, QString *err);
	void signal_invalidate (frac_desc *);
	void signal_compile_kernel (int, int, int, int, QString *err);

	void signal_render (frac_desc *, QGraphicsView *, int);
	void signal_render_preview (frac_desc *, QGraphicsView *, int);
};
