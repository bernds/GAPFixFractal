#include <cstdint>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cfloat>

#include <utility>

#include <QFileDialog>
#include <QThread>
#include <QTimer>
#include <QElapsedTimer>

#include "gpuhandler.h"
#include "formulas.h"
#include "genkernel.h"

inline void tryCuda (CUresult err)
{
	if (err != CUDA_SUCCESS) {
		const char *strerr;

		cuGetErrorName (err, &strerr);
		throw strerr;
	}
}
void GPU_handler::slot_init_cuda (QString *errstr)
{
	int dcount = 0;
	CUresult err = cuInit(0);
	int major = 0, minor = 0;

	if (err == CUDA_SUCCESS)
		err = cuDeviceGetCount (&dcount);

	if (dcount == 0 || err != CUDA_SUCCESS) {
		*errstr = tr ("No CUDA devices found");
		done_sem.release ();
		return;
	}
	err = cuDeviceGet (&m_device, 0);
	if (err != CUDA_SUCCESS)
		*errstr = tr ("Could not obtain CUDA device.");
	else {
		err = cuDeviceGetAttribute (&major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, m_device);
		err = cuDeviceGetAttribute (&minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, m_device);

		err = cuCtxCreate (&m_context, 0, m_device);
		if (err != CUDA_SUCCESS)
			*errstr = tr ("Could not initialize CUDA.");
	}

	done_sem.release ();
}

static uint32_t encode_coord (int x, int y, int w, int h)
{
	uint16_t x1 = x;
	uint16_t y1 = y;
	x1 -= w / 2;
	y1 -= h / 2;
	return (uint32_t)y1 * 65536 + x1;
}

int GPU_handler::initial_setup (frac_desc *fd)
{
	int w = fd->pixel_width;
	int h = fd->pixel_height;
	int pstep = fd->pixel_step;
	int idx = 0;
	for (int y = 0; y < h; y += pstep)
		for (int x = 0; x < w; x += pstep) {
			fd->pixels_started.set_bit (y * w + x);
			fd->host_coords[idx++] = encode_coord (x, y, w, h);
		}
	return idx;
}

int GPU_handler::batch_setup (frac_desc *fd)
{
	int idx = fd->start_idx;
	if (idx == fd->n_threads)
		return idx;
	int w = fd->pixel_width;
	int h = fd->pixel_height;
	int pix_idx = fd->pixels_started.ffz (0);
	int y0 = pix_idx / w;
	int x0 = pix_idx % w;
	for (int y = y0; y < h; y++) {
		for (int x = x0; x < w; x++) {
			fd->pixels_started.set_bit (y * w + x);
			fd->host_coords[idx] = encode_coord (x, y + fd->yoff, w, fd->full_height);
			if (++idx == fd->n_threads)
				return idx;
		}
		x0 = 0;
	}
	return idx;
}

bit_array GPU_handler::compute_ss_pixels (int pstep, int w, const bit_array &done, const bit_array &started)
{
	bit_array ss = done;
	ss.ior (done, -pstep);
	ss.ior (done, pstep);
	ss.ior (done, pstep * w);
	ss.ior (done, -pstep * w);
	ss.andnot (started);
	return ss;
}

int GPU_handler::continue_setup (frac_desc *fd)
{
	int w = fd->pixel_width;
	int h = fd->pixel_height;

	int idx = fd->start_idx;
	if (idx == fd->n_threads)
		return idx;

	int ss = fd->samples;
	int pstep = fd->pixel_step;
	if (pstep > 1 && pstep >= ss) {
		pstep /= 2;
		fd->pixel_step = pstep;
	}
	// The idea here is that we wait for all the regular pixels to start computing,
	// and get into supersampling after that. Priority is given to neighbours of
	// already completed pixels, to quickly improve portions that will end up with
	// a colour.
	// The outer loops are intended to distribute the supersampling effort
	// over the entire picture, rather than doing it in bands.
	if (pstep < ss && fd->pixels_started.popcnt () >= (uint32_t)w * h / ss / ss) {
		bit_array ss_pixels = compute_ss_pixels (pstep, w, fd->pixels_done, fd->pixels_started);
		if (pstep > 1 && ss_pixels.popcnt () < fd->n_threads - fd->start_idx) {
			pstep /= 2;
			fd->pixel_step = pstep;
			ss_pixels = compute_ss_pixels (pstep, w, fd->pixels_done, fd->pixels_started);
		}
		for (int y0 = 0; y0 < ss; y0 += pstep)
			for (int x0 = 0; x0 < ss; x0 += pstep)
				for (int y = y0; y < h; y += ss)
					for (int x = x0; x < w; x += ss) {
						if (!ss_pixels.test_bit (y * w + x))
							continue;
						fd->pixels_started.set_bit (y * w + x);
						fd->host_coords[idx] = encode_coord (x, y, w, h);
						if (++idx == fd->n_threads)
							return idx;
					}
	}

	for (int y = 0; y < h; y += pstep)
		for (int x = 0; x < w; x += pstep) {
			if (fd->pixels_started.test_bit (y * w + x))
				continue;
			fd->pixels_started.set_bit (y * w + x);
			fd->host_coords[idx] = encode_coord (x, y, w, h);
			if (++idx == fd->n_threads)
				return idx;
		}
	return idx;
}

void GPU_handler::slot_start_kernel (frac_desc *fd, int generation, int max_nwords, int steps, bool batch)
{
	{
		QMutexLocker lock (&data_mutex);
		data_available = false;
	}
	if (fd->n_completed == fd->n_pixels)
		abort ();

	bool init_fail = false;
	init_fail |= cuMemcpyHtoD (fd->cu_ar_step, &fd->step[0], 4 * max_nwords) != CUDA_SUCCESS;

	CUdeviceptr originx, originy, param_p, param_q, critpoint, matrix00, matrix01, matrix10, matrix11;
	size_t bytes;
	init_fail |= cuModuleGetGlobal(&originx, &bytes, m_module, "const_origin_x") != CUDA_SUCCESS;
	init_fail |= cuMemcpyHtoD (originx, &fd->center_x[0], 4 * max_nwords) != CUDA_SUCCESS;
	init_fail |= cuModuleGetGlobal(&originy, &bytes, m_module, "const_origin_y") != CUDA_SUCCESS;
	init_fail |= cuMemcpyHtoD (originy, &fd->center_y[0], 4 * max_nwords) != CUDA_SUCCESS;

	init_fail |= cuModuleGetGlobal(&param_p, &bytes, m_module, "const_param_p") != CUDA_SUCCESS;
	init_fail |= cuMemcpyHtoD (param_p, &fd->param_p[0], 4 * max_nwords * 2) != CUDA_SUCCESS;
	init_fail |= cuModuleGetGlobal(&param_q, &bytes, m_module, "const_param_q") != CUDA_SUCCESS;
	init_fail |= cuMemcpyHtoD (param_q, &fd->param_q[0], 4 * max_nwords * 2) != CUDA_SUCCESS;
	init_fail |= cuModuleGetGlobal(&critpoint, &bytes, m_module, "const_critpoint") != CUDA_SUCCESS;
	init_fail |= cuMemcpyHtoD (critpoint, &fd->critpoint[0], 4 * max_nwords * 2) != CUDA_SUCCESS;

	init_fail |= cuModuleGetGlobal(&matrix00, &bytes, m_module, "const_matrix00") != CUDA_SUCCESS;
	init_fail |= cuMemcpyHtoD (matrix00, &fd->matrix[0][0][0], 4 * max_nwords) != CUDA_SUCCESS;
	init_fail |= cuModuleGetGlobal(&matrix01, &bytes, m_module, "const_matrix01") != CUDA_SUCCESS;
	init_fail |= cuMemcpyHtoD (matrix01, &fd->matrix[0][1][0], 4 * max_nwords) != CUDA_SUCCESS;
	init_fail |= cuModuleGetGlobal(&matrix10, &bytes, m_module, "const_matrix10") != CUDA_SUCCESS;
	init_fail |= cuMemcpyHtoD (matrix10, &fd->matrix[1][0][0], 4 * max_nwords) != CUDA_SUCCESS;
	init_fail |= cuModuleGetGlobal(&matrix11, &bytes, m_module, "const_matrix11") != CUDA_SUCCESS;
	init_fail |= cuMemcpyHtoD (matrix11, &fd->matrix[1][1][0], 4 * max_nwords) != CUDA_SUCCESS;

	if (init_fail)
		printf ("init fail\n");

	double iter_scale_factor = 1;

	for (;;) {
		uint32_t count = steps * iter_scale_factor;
		uint32_t maxidx = fd->n_threads;
		int last_n_completed = fd->n_completed;

		if (count > fd->maxiter)
			count = fd->maxiter;
		if (fd->hybrid_len > 0) {
			count = (count + fd->hybrid_len - 1) / fd->hybrid_len * fd->hybrid_len;
		}
		else if (iter_scale_factor != 1)
			printf ("scaling niter from %d to %d\n", steps, count);

		if (batch)
			maxidx = batch_setup (fd);
		else if (fd->n_completed == 0 && fd->start_idx == 0)
			// First-time setup
			// This is an initial pass over fewer pixels to produce some image data quickly
			maxidx = initial_setup (fd);
		else if (fd->n_completed == 0)
			// We've been through at least one round, but no results.  Continue with current set.
			maxidx = fd->start_idx;
		else
			maxidx = continue_setup (fd);

		tryCuda (cuMemcpyHtoD (fd->cu_ar_coords, fd->host_coords, 4 * maxidx));

		uint32_t hybrid_mask = 1;
		hybrid_mask <<= fd->hybrid_len - 1;

		void *args[] = {
			&fd->cu_ar_z, &fd->cu_ar_z2,
			&fd->cu_ar_coords, &fd->cu_ar_step, &fd->cu_ar_tmp,
			&maxidx, &fd->cu_ar_result, &count, &fd->start_idx,
			&fd->hybrid_code, &hybrid_mask,
			&fd->cu_ar_zder
		};

		QElapsedTimer timer;
		timer.start();

		int nthreads = 1024;
		for (;;) {
			int blocks = (maxidx + nthreads - 1) / nthreads;
			CUfunction kernel = (fd->julia
					     ? (fd->hybrid_len != 0 ? m_julia_hybrid : fd->dem ? m_julia_dem : m_julia)
					     : (fd->hybrid_len != 0 ? m_mandel_hybrid : fd->dem ? m_mandel_dem : m_mandel));
			auto err = cuLaunchKernel (kernel,
						   blocks, 1, 1,
						   nthreads, 1, 1,
						   0, 0, args, 0);
			if (err == CUDA_SUCCESS) {
				printf ("launch success %d %d\n", blocks, nthreads);
				break;
			}
			nthreads /= 2;
			if (nthreads < 32) {
				fprintf (stderr, "error launching\n");
				emit signal_new_data (fd, generation, false);
				QMutexLocker lock (&data_mutex);
				return;
			}
		}
		bool fail = false;
		fail |= cuMemcpyDtoH (fd->host_z, fd->cu_ar_z, 4 * fd->nwords * 2 * maxidx) != CUDA_SUCCESS;
		if (fd->dem)
			fail |= cuMemcpyDtoH (fd->host_zder, fd->cu_ar_zder, 4 * fd->nwords * 2 * maxidx) != CUDA_SUCCESS;
		fail |= cuMemcpyDtoH (fd->host_t, fd->cu_ar_tmp, 4 * fd->nwords * 2 * maxidx) != CUDA_SUCCESS;
		fail |= cuMemcpyDtoH (fd->host_z2, fd->cu_ar_z2, 4 * fd->nwords * 2 * maxidx) != CUDA_SUCCESS;
		fail |= cuMemcpyDtoH (fd->host_result, fd->cu_ar_result, 4 * maxidx) != CUDA_SUCCESS;

		qint64 ms = std::max ((qint64)5, timer.elapsed ());
		uint32_t j = 0;
		int w = fd->pixel_width;
		int full_h = fd->full_height;
		int z_size = fd->nwords * 2;
		for (uint32_t i = 0; i < maxidx; i++) {
			uint32_t coord = fd->host_coords[i];
			int32_t hcx = (int16_t)(coord & 65535);
			int32_t hcy = (int16_t)(coord >> 16);
			hcx += w / 2;
			hcy += full_h / 2;
			int result = fd->host_result[i];
			int idx = (hcy - fd->yoff) * w + hcx;
			if (result != 0) {
				fd->pic_result[idx] += result;
				fd->maxiter_found = std::max (fd->maxiter_found, fd->pic_result[idx]);
				fd->pic_z[idx * 2] = to_double (&fd->host_z[i * z_size], fd->nwords);
				fd->pic_z[idx * 2 + 1] = to_double (&fd->host_z[i * z_size + fd->nwords], fd->nwords);
				if (fd->dem) {
					fd->pic_zder[idx * 2] = to_double (&fd->host_zder[i * z_size], fd->nwords);
					fd->pic_zder[idx * 2 + 1] = to_double (&fd->host_zder[i * z_size + fd->nwords], fd->nwords);
				}
				fd->pixels_done.set_bit (idx);
				fd->n_completed++;
			} else {
				fd->pic_result[idx] += count;
				if (fd->pic_result[idx] >= fd->maxiter) {
					fd->pic_result[idx] = 0;
					fd->pic_z[idx * 2] = to_double (&fd->host_z[i * z_size], fd->nwords);
					fd->pic_z[idx * 2 + 1] = to_double (&fd->host_z[i * z_size + fd->nwords], fd->nwords);
					if (fd->dem) {
						fd->pic_zder[idx * 2] = to_double (&fd->host_zder[i * z_size], fd->nwords);
						fd->pic_zder[idx * 2 + 1] = to_double (&fd->host_zder[i * z_size + fd->nwords], fd->nwords);
					}
					fd->pixels_done.set_bit (idx);
					fd->n_completed++;
				} else if (i != j) {
					fd->host_coords[j] = coord;
					memcpy (&fd->host_z2[j * z_size], &fd->host_z2[i * z_size], z_size * 4);
					memcpy (&fd->host_z[j * z_size], &fd->host_z[i * z_size], z_size * 4);
					if (fd->dem)
						memcpy (&fd->host_zder[j * z_size], &fd->host_zder[i * z_size], z_size * 4);
					memcpy (&fd->host_t[j * z_size], &fd->host_t[i * z_size], z_size * 4);
					j++;
				} else
					j++;
			}
		}
		printf ("left with %d out of %d, completed %d vs %d in %d ms\n", j, maxidx, fd->n_completed, fd->n_pixels, (int)ms);
		fd->start_idx = j;

		if (j > 0 && j != maxidx) {
			fail |= cuMemcpyHtoD (fd->cu_ar_tmp, fd->host_t, 4 * fd->nwords * 2 * j) != CUDA_SUCCESS;
			fail |= cuMemcpyHtoD (fd->cu_ar_z, fd->host_z, 4 * fd->nwords * 2 * j) != CUDA_SUCCESS;
			fail |= cuMemcpyHtoD (fd->cu_ar_z2, fd->host_z2, 4 * fd->nwords * 2 * j) != CUDA_SUCCESS;
			if (fd->dem)
				fail |= cuMemcpyHtoD (fd->cu_ar_zder, fd->host_zder, 4 * fd->nwords * 2 * j) != CUDA_SUCCESS;
		}

		QMutexLocker lock (&data_mutex);
		if (abort_computation) {
			data_available = abort_computation = false;
			break;
		}
		data_available = true;
		if (!processing_data && !batch) {
			emit signal_new_data (fd, generation, !fail);
		}
		if (fd->n_completed == fd->n_pixels)
			break;
		if (last_n_completed == 0 && fd->n_completed != 0)
			iter_scale_factor = 1;
		else if (ms < 1000) {
			iter_scale_factor = std::min (iter_scale_factor * 1000. / ms, 100.);
		} else if (ms > 2000)
			iter_scale_factor = std::max (1.0, iter_scale_factor / 2);

	}
	if (batch)
		done_sem.release ();
	else
		emit signal_kernel_complete ();
}

/* The caller should invalidate all existing frac_desc structures before calling this.  */
void GPU_handler::slot_compile_kernel (int fidx, int power, int nwords, int max_nwords, QString *errstr)
{
	if (m_have_module) {
		CUresult err;
		err = cuModuleUnload (m_module);
		if (err != CUDA_SUCCESS) {
			*errstr = tr ("Error unloading previous module");
			done_sem.release ();
			return;
		}
	}

	formula f = formula_table[fidx];
	char *ptx = gen_mprec_funcs (f, nwords, max_nwords, power);

	vector<CUjit_option> opts;
	vector<void *>ovals;

	constexpr int logsz = 8192;
	char error_log[logsz], info_log[logsz];
	int myErr = 0;

	opts.push_back (CU_JIT_INFO_LOG_BUFFER);
	ovals.push_back ((void *)info_log);
	opts.push_back (CU_JIT_INFO_LOG_BUFFER_SIZE_BYTES);
	ovals.push_back ((void *)(uintptr_t)logsz);
	opts.push_back (CU_JIT_ERROR_LOG_BUFFER);
	ovals.push_back ((void *)error_log);
	opts.push_back (CU_JIT_ERROR_LOG_BUFFER_SIZE_BYTES);
	ovals.push_back ((void *)(uintptr_t)logsz);
	opts.push_back (CU_JIT_LOG_VERBOSE);
	ovals.push_back ((void *)1);

	CUlinkState lState;
	tryCuda (cuLinkCreate (opts.size (), &opts[0], &ovals[0], &lState));

	myErr = cuLinkAddData (lState, CU_JIT_INPUT_PTX, ptx, strlen(ptx) + 1, 0, 0, 0, 0);

	if (myErr != CUDA_SUCCESS) {
		fprintf (stderr, "Errors:\n%s\n", error_log);
	}

	void *cuOut;
	tryCuda (cuLinkComplete(lState, &cuOut, nullptr));

	// puts (info_log);

	tryCuda(cuModuleLoadData(&m_module, cuOut));
	tryCuda(cuLinkDestroy(lState));

	// Locate the kernel entry poin
	tryCuda (cuModuleGetFunction (&m_mandel, m_module, "iter_mandel"));
	tryCuda (cuModuleGetFunction (&m_julia, m_module, "iter_julia"));
	if (formula_supports_dem (f)) {
		tryCuda (cuModuleGetFunction (&m_mandel_dem, m_module, "iter_mandel_dem"));
		tryCuda (cuModuleGetFunction (&m_julia_dem, m_module, "iter_julia_dem"));
	} else {
		m_mandel_dem = 0;
		m_julia_dem = 0;
	}
	if (formula_supports_hybrid (f)) {
		tryCuda (cuModuleGetFunction (&m_mandel_hybrid, m_module, "iter_mandel_hybrid"));
		tryCuda (cuModuleGetFunction (&m_julia_hybrid, m_module, "iter_julia_hybrid"));
	} else {
		m_mandel_hybrid = 0;
		m_julia_hybrid = 0;
	}

	free (ptx);
	m_have_module = true;
	done_sem.release ();
}

void GPU_handler::slot_alloc_mem (frac_desc *fd, int max_nwords, int nwords, int w, int h, QString *errstr)
{
	try {
		int nthreads = w * h;
		if (fd->cu_ar_origin == 0) {
			tryCuda (cuMemAlloc (&fd->cu_ar_origin, 4 * nwords * 2 * nthreads));
			tryCuda (cuMemAlloc (&fd->cu_ar_z, 4 * nwords * 2 * nthreads));
			if (fd->dem)
				tryCuda (cuMemAlloc (&fd->cu_ar_zder, 4 * nwords * 2 * nthreads));
			tryCuda (cuMemAlloc (&fd->cu_ar_z2, 4 * nwords * 2 * nthreads));
			tryCuda (cuMemAlloc (&fd->cu_ar_step, 4 * max_nwords));
			tryCuda (cuMemAlloc (&fd->cu_ar_tmp, 4 * nwords * 2 * nthreads));
			tryCuda (cuMemAlloc (&fd->cu_ar_result, 4 * nthreads));
			tryCuda (cuMemAlloc (&fd->cu_ar_coords, 4 * nthreads));
		}
	} catch (const char *err) {
		*errstr = err;
		free_cuda_data (fd);
	}
	done_sem.release ();
}

void GPU_handler::free_cuda_data (frac_desc *fd)
{
	if (fd->n_pixels == 0 || fd->cu_ar_origin == 0)
		return;

	cuMemFree (fd->cu_ar_origin);
	cuMemFree (fd->cu_ar_z);
	if (fd->cu_ar_zder != 0)
		cuMemFree (fd->cu_ar_zder);
	cuMemFree (fd->cu_ar_z2);
	cuMemFree (fd->cu_ar_step);
	cuMemFree (fd->cu_ar_coords);
	cuMemFree (fd->cu_ar_tmp);
	cuMemFree (fd->cu_ar_result);
	fd->cu_ar_origin = 0;
	fd->cu_ar_z = 0;
	fd->cu_ar_zder = 0;
	fd->cu_ar_z2 = 0;
	fd->cu_ar_step = 0;
	fd->cu_ar_coords = 0;
	fd->cu_ar_tmp = 0;
	fd->cu_ar_result = 0;
}
