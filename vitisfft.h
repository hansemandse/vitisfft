
#ifndef FFT_H
#define FFT_H

#include <complex>
#include "ap_fixed.h"
#include "hls_fft.h"
#include "hls_stream.h"

using namespace std;

// ----------------------------------------------------------------------------
// Various definitions
// Include this definition to use floats instead of 16-bit fixed-point
#define FFT_USE_FLOATS

// Most of these values can be adjusted as needed
#define FFT_INPUT_W 16           // bit-width of inputs to the FFT
#define FFT_OUTPUT_W FFT_INPUT_W // bit-width of outputs from the FFT
#ifdef FFT_USE_FLOATS
#define FFT_PHASE_FACTOR_W 24    // bit-width of phase factors in the FFT
#else
#define FFT_PHASE_FACTOR_W 16

constexpr unsigned gen_sch()
{
	static_assert(
		FFT_ARCH_OPT == hls::ip_fft::pipelined_streaming_io ||
			FFT_ARCH_OPT == hls::ip_fft::radix_2_burst_io,
		"can only instantiate pipelined streaming or radix-2 burst FFT IP");

	unsigned sch{}, stages{};
	switch (FFT_ARCH_OPT)
	{
	case hls::ip_fft::pipelined_streaming_io:
		// When using a radix-4 algorithm, the IP has (nfft+1)/2 stages that
		// each requires at least one and typically two bit shifts (except the
		// last one)
		stages = (FFT_NFFT_MAX + 1) / 2;
		for (int i{}; i < stages; ++i)
		{
			if (i == 0)
				sch |= (FFT_NFFT_MAX & 0x1);
			else if (i == stages - 1)
				sch = (sch << 2) | 0x3;
			else
				sch = (sch << 2) | 0x2;
		}
		break;
	case hls::ip_fft::radix_2_burst_io:
		// When using a radix-2 algorithm, the IP has nfft stages that each
		// requires at most one bit shift (except the last one)
		stages = FFT_NFFT_MAX;
		for (int i{}; i < stages - 1; ++i)
		{
			if (i == stages / 2)
				sch <<= 2;
			else
				sch = (sch << 2) | 0x1;
		}
		sch = (sch << 2) | 0x2;
		break;
	}

	return sch;
}

constexpr unsigned FFT_SCH = gen_sch();
#endif
#define FFT_NFFT_MAX 14 // log2(max FFT size) (2^14 = 16384)
#define FFT_LENGTH_MAX (1 << FFT_NFFT_MAX)
#define FFT_CONFIG_W (2 * FFT_NFFT_MAX + 4) // bit-width of the configuration
#define FFT_ARCH_OPT hls::ip_fft::radix_2_burst_io

// Doubles are impossible to synthesize, floats are very costly in area.
// Seems like 16-bit fixed-point values with scaling are a suitable alternative.
#ifdef FFT_USE_FLOATS
typedef float fft_data_in_t;
typedef float fft_data_out_t;
#else
typedef ap_fixed<FFT_INPUT_W, 1> fft_data_in_t;
typedef ap_fixed<FFT_OUTPUT_W, FFT_OUTPUT_W - FFT_INPUT_W + 1> fft_data_out_t;
#endif
typedef complex<fft_data_in_t> cmplx_fft_data_in_t;
typedef complex<fft_data_out_t> cmplx_fft_data_out_t;

struct config_str : hls::ip_fft::params_t
{
	// Outputs ordered naturally as the inputs
	static const unsigned ordering_opt = hls::ip_fft::natural_order;

	// Outputs rounded to avoid large DC bias
	static const unsigned rounding_opt = hls::ip_fft::convergent_rounding;

	// Good balance of throughput, area, and controllability with scaling
	static const unsigned arch_opt = FFT_ARCH_OPT;

	// Optimize for resources with DSPs
	static const unsigned complex_mult_type = hls::ip_fft::use_mults_resources;

	// Change internal memory structure
	static const bool mem_hybrid = true;
	static const unsigned mem_phase_factors = hls::ip_fft::distributed_ram;
	static const unsigned mem_reorder = hls::ip_fft::distributed_ram;

	// Depends on the data_in_t: if float, needs to be either 24 or 25,
	// otherwise it can be scaled to, e.g., 16
	static const unsigned phase_factor_width = FFT_PHASE_FACTOR_W;

	// Also depends on the other parameters passed to the core. Scaling
	// schedule bits take a lot of space.
	static const unsigned config_width = FFT_CONFIG_W;

	// Maximum log2 FFT size and support for runtime configuration
	static const unsigned max_nfft = FFT_NFFT_MAX;
};

typedef hls::ip_fft::config_t<config_str> config_t;
typedef hls::ip_fft::status_t<config_str> status_t;

// ----------------------------------------------------------------------------
// Provide a set of complex FFT-related functions
void cfft(complex<float>[FFT_LENGTH_MAX], complex<float>[FFT_LENGTH_MAX], const size_t);
void cifft(complex<float>[FFT_LENGTH_MAX], complex<float>[FFT_LENGTH_MAX], const size_t);

#endif
