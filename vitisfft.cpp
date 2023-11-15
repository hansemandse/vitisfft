
#include "vitisfft.h"
#include <cassert>
#include "hls_print.h"

// ----------------------------------------------------------------------------
// Top-level FFT function (not meant for calls outside this file!)

void _fft_ip(
	cmplx_data_in_t* xn,
	cmplx_data_out_t* xk,
	config_t* fft_config)
{
#pragma HLS interface mode=ap_fifo port=xn, xk
#pragma HLS dataflow

	// Create the local status buffer
	status_t fft_status;

	// Run the FFT on the given buffers
	hls::fft<config_str>(xn, xk, &fft_status, fft_config);
#ifndef __SYNTHESIS__
	assert(!(fft_status.getOvflo() & 0x1));
#endif
	return;
}

// ----------------------------------------------------------------------------
// Helpers to simplify the interface

void _data_input_mover(
	double* in,
	cmplx_data_in_t* xn,
	const size_t size,
	const bool dir,
	const bool real)
{
	// Move the signal
	MOVE_IN : for (size_t i {}; i < size; ++i)
	{
		if (real && dir) // real forward
		{
			xn[i].real(static_cast<data_in_t>(in[i]));
			xn[i].imag(static_cast<data_in_t>(0.0));
		}
		else // others
		{
			xn[i].real(static_cast<data_in_t>(in[2*i]));
			xn[i].imag(static_cast<data_in_t>(in[2*i+1]));
		}
	}

	// Zero-pad the signal
	ZEROPAD_IN : for (size_t i {size}; i < FFT_LENGTH; ++i)
	{
		xn[i].real(static_cast<data_in_t>(0.0));
		xn[i].imag(static_cast<data_in_t>(0.0));
	}
	return;
}

void _data_output_mover(
	double* out,
	cmplx_data_out_t* xk,
	const unsigned nfft,
	const bool dir,
	const bool real)
{
	// Compute the output size
	const size_t size = 1 << nfft;

	// Move the signal
	MOVE_OUT : for (size_t i {}; i < size; ++i)
	{
		if (real && !dir) // real inverse
		{
			out[i] = static_cast<double>(xk[i].real());
		}
		else // others
		{
			out[2*i]   = static_cast<double>(xk[i].real());
			out[2*i+1] = static_cast<double>(xk[i].imag());
		}
	}

	// Zero-pad the signal
	if (real && !dir)
	{
		ZEROPAD_OUT : for (size_t i {size}; i < 2 * size; ++i) out[i] = 0;
	}
	return;
}

// ----------------------------------------------------------------------------

unsigned log2up(size_t num)
{
	unsigned msbp {}, cnt {};
	LOG2UP : while ((num >> 1) > 0)
	{
		msbp += 1;
		cnt  += (num & 0x1);
		num >>= 1;
	}
	return msbp + ((cnt > 1) ? 1 : 0);
}

void _fft(const bool dir, const bool real, double* signal, const size_t size)
{
	// Compute the log2 FFT size ...
	const unsigned nfft = log2up(size);
#ifndef FFT_USE_FLOATS
	// ... and use it to establish a conservative FFT scaling schedule
	// according to AMD PG109 ps. 29-30, 43-44, and 48. In general, the
	// scaling schedule's entries should sum to nfft.
	unsigned sch {}, stages {};
	switch (FFT_ARCH_OPT)
	{
	case hls::ip_fft::pipelined_streaming_io:
		// When using a radix-4 algorithm, the IP has (nfft+1)/2 stages that
		// each requires at least one and typically two bit shifts (except the
		// last one)
		stages = (nfft+1) / 2;
		PSIO_SCH : for (int i {}; i < stages; ++i)
		{
			if (i == 0) sch |= (nfft & 0x1);
			else if (i == stages-1) sch = (sch << 2) | 0x3;
			else sch = (sch << 2) | 0x2;
		}
		break;
	case hls::ip_fft::radix_2_burst_io:
		// When using a radix-2 algorithm, the IP has nfft stages that each
		// requires at most one bit shift (except the last one)
		stages = nfft;
		R2BIO_SCH : for (int i {}; i < stages-1; ++i)
		{
			if (i == stages/2) sch <<= 2;
			else sch = (sch << 2) | 0x1;
		}
		sch = (sch << 2) | 0x2;
		break;
	default:
		hls::print("Unsupported FFT architecture option %i!\n", FFT_ARCH_OPT);
		assert(0);
	}
#endif

	// Set-up the configuration for this call
	config_t fft_config;
	fft_config.setDir(dir);
	fft_config.setNfft(nfft);
#ifndef FFT_USE_FLOATS
	fft_config.setSch(sch);
#endif

	// Create the right-size containers for input and output to the FFT
	// These must be size (1 << FFT_NFFT_MAX) according to the following
	// https://support.xilinx.com/s/question/0D52E00006iHvpeSAC/sigsegv-with-fft-hls?language=en_US
#ifndef __SYNTHESIS__
	cmplx_data_in_t  *xn = new cmplx_data_in_t [FFT_LENGTH];
	cmplx_data_out_t *xk = new cmplx_data_out_t[FFT_LENGTH];
#else
	cmplx_data_in_t  xn[FFT_LENGTH];
	cmplx_data_out_t xk[FFT_LENGTH];
#endif

	// Cast and move the data into the input container
	_data_input_mover(signal, xn, size, dir, real);

	// Pass the containers to the FFT and retrieve its result
	_fft_ip(xn, xk, &fft_config);

	// Cast and move the data from the output container into the signal container
	_data_output_mover(signal, xk, nfft, dir, real);

#ifndef __SYNTHESIS__
	delete [] xk;
	delete [] xn;
#endif

	return;
}

// ----------------------------------------------------------------------------

// Complex forward FFT (input signal has size complex value pairs)
void vitis_cfft(double* signal, const size_t size)
{
	_fft(true, false, signal, size);
	return;
}

// Complex inverse FFT (input signal has size complex value pairs)
void vitis_cifft(double *signal, const size_t size)
{
	_fft(false, false, signal, size);
	return;
}

// Real forward FFT (input signal has size real values, and the function
// returns the upper half of the resulting frequency spectrum)
void vitis_rfft(double *signal, const size_t size)
{
	_fft(true, true, signal, size);
	return;
}

// Real inverse FFT (input signal has size complex value pairs representing
// the upper half of the frequency spectrum)
void vitis_rifft(double *signal, const size_t size)
{
	_fft(false, true, signal, size);
	return;
}
