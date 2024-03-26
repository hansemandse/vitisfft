
#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include "vitisfft.h"

// Enable debug printing of signals and FFT results to files
#define FEXPORT

int main()
{
// ----------------------------------------------------------------------------
// Setup for all the tests
	auto cabs = [](complex<double>& c) { return std::sqrt(real(c) * real(c) + imag(c) * imag(c)); };
	const size_t DataSize {1 << 10};
	const double Scale {1. / DataSize};

	const double SamplingFreq {44100};	// in Hz
	const double SineFreq {200};		// in Hz
	const double StepSize {1.0 / SamplingFreq};	// in s

	double max {}, diff {};
	size_t maxind {};
	complex<double> cmplx;
	double sine[2 * FFT_LENGTH];

// ----------------------------------------------------------------------------
// Compute the forward complex FFT on a complex exponential signal
	for (size_t i {}; i < DataSize; ++i)
	{
		sine[2*i]   = cos(2*M_PI*SineFreq*i*StepSize);
		sine[2*i+1] = sin(2*M_PI*SineFreq*i*StepSize);
	}

	double fcsine[2 * DataSize];
	for (size_t i {}; i < 2 * DataSize; ++i) fcsine[i] = sine[i];
	cfft(fcsine, DataSize);
	for (size_t i {}; i < 2 * DataSize; ++i) fcsine[i] *= Scale;

	// Find the maximum frequency index
	cmplx.real(fcsine[0]);
	cmplx.imag(fcsine[1]);
	max = cabs(cmplx);
	for (size_t i {1}; i < DataSize; ++i)
	{
		cmplx.real(fcsine[2*i]);
		cmplx.imag(fcsine[2*i+1]);
		if (cabs(cmplx) > max) {
			maxind = i;
			max    = cabs(cmplx);
		}
	}

	// Check that the maximum frequency index matches the expectation
	if (round(SineFreq / (SamplingFreq / DataSize)) != maxind) {
		cout << "Maximum frequency index " << maxind << " does not match expected ";
		cout << round(SineFreq / (SamplingFreq / DataSize)) << endl;
		return 2;
	}

	double icsine[2 * DataSize];
	for (size_t i {}; i < 2 * DataSize; ++i) icsine[i] = fcsine[i];
	vitis_cifft(icsine, DataSize);

#ifdef FEXPORT
	ofs.open("cfft.csv");
	ofs << "Sine,CFFT,CIFFT" << endl;
	for (size_t i {}; i < 2 * DataSize; ++i)
		ofs << sine[i] << "," << fcsine[i] << "," << icsine[i] << endl;
	ofs.close();
#endif

	// Check that the waveform matches the expectation
	diff = 0;
	for (size_t i {}; i < 2 * DataSize; ++i) diff += abs(sine[i] - icsine[i]);
	if (diff > 1.0) {
		cout << "Waveform error " << diff << " exceeds the expected " << 1.0 << endl;
		return 2;
	}
}
