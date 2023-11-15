
#ifndef FFT_ALT_H
#define FFT_ALT_H

#include <complex>
#include "ap_fixed.h"
#include "hls_fft.h"

using namespace std;

// ----------------------------------------------------------------------------
// Various type definitions
// Include this definition to use floats instead of 16-bit fixed-point
#define FFT_USE_FLOATS

// NOTE: these can be adjusted as we see fit. I have picked values that
//		 correspond roughly to what seems fitting for resource consumption
//		 and precision
#define FFT_INPUT_WIDTH 		16						// bit-width of inputs to the FFT
#define FFT_OUTPUT_WIDTH		FFT_INPUT_WIDTH			// bit-width of outputs from the FFT
#ifdef FFT_USE_FLOATS
#define FFT_PHASE_FACTOR_WIDTH	24						// bit-width of phase factors in the FFT
#else
#define FFT_PHASE_FACTOR_WIDTH	16
#endif
#define FFT_NFFT_MAX			14						// log2(max FFT size) (2^14 = 16384)
#define	FFT_LENGTH				(1 << FFT_NFFT_MAX)
#define FFT_CONFIG_WIDTH		(2*FFT_NFFT_MAX + 12)	// bit-width of the configuration
#define FFT_ARCH_OPT			hls::ip_fft::radix_2_burst_io

// Doubles are impossible to synthesize, floats are very costly in area.
// Seems like 16-bit fixed-point values with scaling are a suitable alternative.
#ifdef FFT_USE_FLOATS
typedef float data_in_t;
typedef float data_out_t;
#else
typedef ap_fixed<FFT_INPUT_WIDTH, 1> data_in_t;
typedef ap_fixed<FFT_OUTPUT_WIDTH, FFT_OUTPUT_WIDTH - FFT_INPUT_WIDTH + 1> data_out_t;
#endif
typedef hls::x_complex<data_in_t> cmplx_data_in_t;
typedef hls::x_complex<data_out_t> cmplx_data_out_t;

struct config_str : hls::ip_fft::params_t {
	// Outputs ordered naturally as the inputs
	static const unsigned ordering_opt = hls::ip_fft::natural_order;

	// Outputs rounded to avoid large DC bias
	static const unsigned rounding_opt = hls::ip_fft::convergent_rounding;

	// Good balance of throughput, area, and controllability with scaling
	static const unsigned arch_opt = FFT_ARCH_OPT;

	// Depends on the data_in_t: if float, needs to be either 24 or 25,
	// otherwise it can be scaled to, e.g., 16
	static const unsigned phase_factor_width = FFT_PHASE_FACTOR_WIDTH;

	// Also depends on the other parameters passed to the core. Scaling
	// schedule bits take a lot of space. Seems like 2*NFFT_MAX+12 is
	// sufficient for all tested configurations
	static const unsigned config_width = FFT_CONFIG_WIDTH;

	// Maximum log2 FFT size and support for runtime configuration
	static const unsigned max_nfft = FFT_NFFT_MAX;
	static const bool 	  has_nfft = true;
};

typedef hls::ip_fft::config_t<config_str> config_t;
typedef hls::ip_fft::status_t<config_str> status_t;

// ----------------------------------------------------------------------------
/* Provide a full set of FFT-related functions
 *
 * These functions can be summarized as follows:
 * - vitis_cfft computes the complex forward FFT, assuming `signal` contains
 *   `size` complex values (i.e., length 2 * `size`), and returns the full
 *   frequency spectrum of 2 * `2 ^ log2up(size)` complex values
 *
 * - vitis_cifft computes the complex inverse FFT, assuming `signal` contains
 *   a full frequency spectrum of `size` complex values (i.e., length 2 * `size`),
 *   and returns the time domain signal of 2 * `2 ^ log2up(size)` complex values
 *
 * - vitis_rfft computes the real forward FFT, assuming `signal` contains
 *   `size` real values, and returns the full frequency spectrum of
 *   2 * `2 ^ log2up(size)` complex values
 *
 * - vitis_rifft computes the real inverse FFT, assuming `signal` contains
 *   a full frequency spectrum of `size` complex values (i.e., has length
 *   2 * `size`), and returns the time domain signal of `2 ^ log2up(size)`
 *   real values
 *
 * Note that with these follow an assumption of `signal` always being
 * allocated at least 2 * `2 ^ log2up(size)` elements. All FFTs computed for
 * `size` < `2 ^ log2up(size)` are zero-padded.
 */
void vitis_cfft(double*, const size_t);
void vitis_cifft(double*, const size_t);
void vitis_rfft(double*, const size_t);
void vitis_rifft(double*, const size_t);

// Provide a convenient function to find the next power-of-2
unsigned log2up(size_t);

#endif
