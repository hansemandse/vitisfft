
#include "vitisfft.h"
#include <cassert>

// ----------------------------------------------------------------------------
// Helpers to simplify the interface

void _data_in_mover(
	bool dir,
	const size_t size,
	hls::stream<config_t>& fft_config,
	complex<float> in[FFT_LENGTH_MAX],
	hls::stream<cmplx_fft_data_in_t>& xn)
{
	// Set up the configuration
	config_t conf;
	conf.setDir(dir);
#ifndef FFT_USE_FLOATS
	// Establish a conservative FFT scaling schedule according to 
	// AMD PG109 ps. 29-30, 43-44, and 48. In general, the scaling 
	// schedule's entries should sum to nfft.
	conf.setSch(FFT_SCH);
#endif
	fft_config.write(conf);

	// Move the signal
	MOVE_IN : for (size_t i {}; i < size; ++i)
#pragma hls pipeline ii=1
		xn.write(static_cast<cmplx_fft_data_in_t>(in[i]));

	// Zero-pad the signal
	ZEROPAD_IN : for (size_t i {size}; i < FFT_LENGTH_MAX; ++i)
#pragma hls pipeline ii=1
		xn.write(0);
	return;
}

void _data_out_mover(
	hls::stream<status_t>& fft_status,
	hls::stream<cmplx_fft_data_out_t>& xk,
	complex<float> out[FFT_LENGTH_MAX])
{
	status_t status = fft_status.read();
#ifndef __SYNTHESIS__
	// Check the overflow flag
	assert(!(status.getOvflo() & 0x1));
#endif

	// Move the signal
	MOVE_OUT : for (size_t i {}; i < FFT_LENGTH_MAX; ++i)
#pragma hls pipeline ii=1
		out[i] = static_cast<complex<float>>(xk.read());
	return;
}

// ----------------------------------------------------------------------------

void _fft(
	bool dir,
	complex<float> in [FFT_LENGTH_MAX],
	complex<float> out[FFT_LENGTH_MAX],
	const size_t size)
{
#pragma hls dataflow

	// Create the containers for input and output to the FFT
	hls::stream<cmplx_fft_data_in_t>  xn;
	hls::stream<cmplx_fft_data_out_t> xk;
	hls::stream<config_t> fft_config;
	hls::stream<status_t> fft_status;

	_data_in_mover(dir, size, fft_config, in, xn);
	hls::fft<config_str>(xn, xk, fft_status, fft_config);
	_data_out_mover(fft_status, xk, out);
}

// ----------------------------------------------------------------------------

// Complex forward FFT
void cfft(
	complex<float> in[FFT_LENGTH_MAX],
	complex<float> out[FFT_LENGTH_MAX],
	const size_t size)
{
#pragma hls inline off
	_fft(true, in, out, size);
	return;
}

// Complex inverse FFT
void cifft(
	complex<float> in[FFT_LENGTH_MAX],
	complex<float> out[FFT_LENGTH_MAX],
	const size_t size)
{
#pragma hls inline off
	_fft(false, in, out, size);
	return;
}
