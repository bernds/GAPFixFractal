#include <QMainWindow>
#include <QMutex>
#include <QGraphicsScene>
#include <QImage>
#include <QTimer>
#include <QThreadPool>

#include "fpvec.h"
#include "fractal.h"

namespace Ui
{
	class MainWindow;
};

class ClickablePixmap;
class QActionGroup;

struct stored_preview {
	ClickablePixmap *pixmap {};
	frac_params params;
	QImage thumbnail;
	int x = 0, y = 0;

	stored_preview (const frac_params &p, QImage tn)
		: params (p), thumbnail (tn)
	{
	}
	stored_preview (const stored_preview &) = default;
	stored_preview (stored_preview &&) = default;
	stored_preview &operator= (const stored_preview &) = default;
	stored_preview &operator= (stored_preview &&) = default;
};

struct render_params
{
	QVector<uint32_t> palette;
	int incol;
	int mod_type;
	double steps;
	int slider;
	bool sub;
	bool angle;
	double dem_param;
};

class Renderer : public QObject
{
	Q_OBJECT
	QThreadPool m_pool;
	double m_minimum = 0;
	int m_min_gen = -1;

	void do_render (const render_params &rp, int w, int h, frac_desc *, QGraphicsView *, int);

public:

	// One big mutex around the drawing function
	QMutex mutex;

	// Communication with Mainwindow
	QMutex queue_mutex;
	bool queued = false;
	render_params next_rp;
	int render_width, render_height;

	void slot_render (frac_desc *, QGraphicsView *, int);

signals:
	void signal_render_complete (QGraphicsView *, frac_desc *, QImage);
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
	static constexpr int maxiter = 500000;

	formula m_formula = formula::standard;

	frac_desc m_fd_mandel;
	frac_desc m_fd_julia;

	QImage m_img_mandel, m_img_julia;

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
	QActionGroup *m_power_group {};
	QActionGroup *m_struct_group {};
	QActionGroup *m_narrow_group {};

	int m_last_pal_idx = 0;
	QVector<uint32_t> m_palette;
	QVector<uint32_t> m_custom_palette;

	void restore_geometry ();
	void start_threads ();
	void do_compile ();

	void perform_resizes ();

	void build_points (frac_desc &, int w, int h);
	void compute_fractal (frac_desc &, int nwords, int w, int h, int ss, bool preview);
	void render_fractal ();
	void render_preview ();
	bool abort_computation ();
	void restart_computation ();

	void autoprec (frac_desc &);
	void update_settings (bool);
	// double iter_value_at (frac_desc &, int);
	void update_display (QGraphicsView *);
	void precompute_iter_value (frac_desc *);
	void update_aspect ();
	void update_views (int = 0);
	const QVector<uint32_t> &palette_from_index (int);
	void update_palette ();
	void update_fractal_type (int = 0);
	void adjust_width_for_bounds (frac_desc &);
	void reset_coords (frac_desc &);
	void set_q (int, int);
	frac_desc &current_fd ();
	void do_pause (bool);
	void do_reset (bool);
	void zoom_in (bool = false);
	void zoom_out (bool = false);

	void set_render_params (render_params &);
        void slot_new_data (frac_desc *, int, bool);
	void slot_kernel_complete ();
	void slot_render_complete (QGraphicsView *view, frac_desc *, QImage);

	void layout_stored_params ();
	void store_params (bool);
	void restore_params (const frac_params &);

	void slot_save_as (bool);
	void slot_save_params (bool);
	void slot_load_params (bool);
	void slot_save_palette (bool);
	void slot_load_palette (bool);

	void gradient_edit (bool);

	void fractal_mouse_event (QMouseEvent *);
	void fractal_wheel_event (QWheelEvent *);
	void preview_wheel_event (QWheelEvent *);

	void init_formula (formula);
	void formula_chosen (formula, int);

	void help_about ();
protected:
	void closeEvent(QCloseEvent *event) override;

public:
	MainWindow ();
	~MainWindow ();

signals:
	void signal_init_cuda (QString *err);
	void signal_start_kernel (frac_desc *, int generation, int, int);
	void signal_alloc_mem (frac_desc *, int, int, int, int, QString *err);
	void signal_invalidate (frac_desc *);
	void signal_compile_kernel (int, int, int, int, QString *err);

	void signal_render (frac_desc *, QGraphicsView *, int);
	void signal_render_preview (frac_desc *, QGraphicsView *, int);
};
