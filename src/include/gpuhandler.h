#include <cuda.h>

#include <QObject>
#include <QSemaphore>
#include <QMutex>

#include "fpvec.h"
#include "fractal.h"
#include "bitarray.h"

class GPU_handler : public QObject
{
        Q_OBJECT

	CUmodule m_module;
	CUdevice m_device;
	CUcontext m_context;
	bool m_have_module = false;
	CUfunction m_mandel, m_julia;
	CUfunction m_mandel_dem, m_julia_dem;
	CUfunction m_mandel_hybrid, m_julia_hybrid;

	void free_cuda_data (frac_desc *);
	int initial_setup (frac_desc *);
	int continue_setup (frac_desc *);
	int batch_setup (frac_desc *);

	bit_array compute_ss_pixels (int, int, const bit_array &, const bit_array &);

public:
	/* Used to communicate with the owner (MainWindow).  */
	QSemaphore done_sem;
	QMutex data_mutex;
	bool data_available = false;
	bool processing_data = false;
	bool abort_computation = false;

public slots:
	void slot_init_cuda (QString *);
	void slot_alloc_mem (frac_desc *, int max_nwords, int nwords, int w, int h, QString *);
	void slot_invalidate (frac_desc *fd)
	{
		free_cuda_data (fd);
		done_sem.release ();
	}
	void slot_start_kernel (frac_desc *, int generation, int max_nwords, int steps, bool batch);
	void slot_compile_kernel (int, int power, int off, int nwords, int max_nwords, bool incolor, QString *);
signals:
        void signal_new_data (frac_desc *, int generation, bool);
        void signal_kernel_complete ();
};

