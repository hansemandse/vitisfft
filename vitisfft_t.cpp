
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
// Compute the forward real-valued FFT on a single generated sine wave (see PG109 p. 33)
	for (size_t i {}; i < DataSize; ++i) sine[i] = sin(2*M_PI*SineFreq*i*StepSize);

	double fsine[2 * DataSize];
	for (size_t i {}; i < DataSize; ++i) fsine[i] = sine[i];
	vitis_rfft(fsine, DataSize);
	for (size_t i {}; i < 2 * DataSize; ++i) fsine[i] *= Scale;

	// Find the maximum frequency index
	cmplx.real(fsine[0]);
	cmplx.imag(fsine[1]);
	max = cabs(cmplx);
	for (size_t i {1}; i < DataSize; ++i)
	{
		cmplx.real(fsine[2*i]);
		cmplx.imag(fsine[2*i+1]);
		if (cabs(cmplx) > max) {
			maxind = i;
			max	   = cabs(cmplx);
		}
	}

	// Check that the maximum frequency index matches the expectation
	if (round(SineFreq / (SamplingFreq / DataSize)) != maxind) {
		cout << "Maximum frequency index " << maxind << " does not match expected ";
		cout << round(SineFreq / (SamplingFreq / DataSize)) << endl;
		return 1;
	}

	double isine[2 * DataSize];
	for (size_t i {}; i < 2 * DataSize; ++i) isine[i] = fsine[i];
	vitis_rifft(isine, DataSize);

#ifdef FEXPORT
	ofstream ofs;
	ofs.open("rfft.csv");
	ofs << "Sine,RFFT,IRFFT" << endl;
	for (size_t i {}; i < 2 * DataSize; ++i)
		ofs << (i < DataSize ? sine[i] : 0) << "," << fsine[i] << "," << (i < DataSize ? isine[i] : 0) << endl;
	ofs.close();
#endif

	// Check that the waveform matches the expectation
	diff = 0;
	for (size_t i {}; i < DataSize; ++i) diff += abs(sine[i] - isine[i]);
	if (diff > 1.0) {
		cout << "Waveform error " << diff << " exceeds the expected " << 1.0 << endl;
		return 1;
	}

// ----------------------------------------------------------------------------
// Perform the same experiment using a complex signal and FFT
	for (size_t i {}; i < DataSize; ++i)
	{
		sine[2*i]   = cos(2*M_PI*SineFreq*i*StepSize);
		sine[2*i+1] = sin(2*M_PI*SineFreq*i*StepSize);
	}

	double fcsine[2 * DataSize];
	for (size_t i {}; i < 2 * DataSize; ++i) fcsine[i] = sine[i];
	vitis_cfft(fcsine, DataSize);
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

// ----------------------------------------------------------------------------
// Perform an experiment on three generated sine waves and the real FFT
	const double Sine2Freq {750};
	const double Sine3Freq {1023};
	for (size_t i {}; i < DataSize; ++i)
		sine[i] = (
			sin(2*M_PI*SineFreq *i*StepSize) +
			sin(2*M_PI*Sine2Freq*i*StepSize) +
			sin(2*M_PI*Sine3Freq*i*StepSize)) / 3;

	for (size_t i {}; i < DataSize; ++i) fsine[i] = sine[i];
	vitis_rfft(fsine, DataSize);
	for (size_t i {}; i < 2 * DataSize; ++i) fsine[i] *= Scale;

	// Find the maximum frequency indices
	cmplx.real(fsine[0]);
	cmplx.imag(fsine[1]);
	complex<double> cmplx2 (fsine[2], fsine[3]);
	complex<double> cmplx3 (fsine[4], fsine[5]);
	tuple<double, size_t> maxes[3];
	if (cabs(cmplx) > cabs(cmplx2) && cabs(cmplx) > cabs(cmplx3)) {
		maxes[0] = {cabs(cmplx), 0};
		if (cabs(cmplx2) > cabs(cmplx3)) {
			maxes[1] = {cabs(cmplx2), 1};
			maxes[2] = {cabs(cmplx3), 2};
		} else {
			maxes[1] = {cabs(cmplx3), 2};
			maxes[2] = {cabs(cmplx2), 1};
		}
	} else if (cabs(cmplx) > cabs(cmplx2)) {
		maxes[0] = {cabs(cmplx3), 2};
		maxes[1] = {cabs(cmplx),  0};
		maxes[2] = {cabs(cmplx2), 1};
	} else if (cabs(cmplx) > cabs(cmplx3)) {
		maxes[0] = {cabs(cmplx2), 1};
		maxes[1] = {cabs(cmplx),  0};
		maxes[2] = {cabs(cmplx3), 2};
	} else {
		maxes[2] = {cabs(cmplx), 0};
		if (cabs(cmplx2) > cabs(cmplx3)) {
			maxes[0] = {cabs(cmplx2), 1};
			maxes[1] = {cabs(cmplx3), 2};
		} else {
			maxes[0] = {cabs(cmplx3), 2};
			maxes[1] = {cabs(cmplx2), 1};
		}
	}
	for (size_t i {3}; i < DataSize / 2; ++i)
	{
		cmplx.real(fsine[2*i]);
		cmplx.imag(fsine[2*i+1]);
		for (size_t j {}; j < 3; ++j)
		{
			if (cabs(cmplx) > std::get<0>(maxes[j])) {
				for (size_t k {2}; k > j; --k) maxes[k] = maxes[k-1];
				maxes[j] = {cabs(cmplx), i};
				break;
			}
		}
	}

	// Check that the maximum frequency indices match the expectation
	auto freq_found = [maxes, SamplingFreq, DataSize](const double freq) {
		bool res {false};
		for (size_t i {}; i < 3; ++i) res = res || (get<1>(maxes[i]) == round(freq / (SamplingFreq / DataSize)));
		return res;
	};
	bool freq1found = freq_found(SineFreq);
	bool freq2found = freq_found(Sine2Freq);
	bool freq3found = freq_found(Sine3Freq);
	if (!(freq1found && freq2found && freq3found))
	{
		cout << "Maximum frequency indices " << std::get<1>(maxes[0]) << ", " << get<1>(maxes[1]);
		cout << ", and " << get<1>(maxes[2]) << " do not match expected ";
		cout << round(SineFreq / (SamplingFreq / DataSize)) << ", " << round(Sine2Freq / (SamplingFreq / DataSize));
		cout << ", and " << round(Sine3Freq / (SamplingFreq / DataSize)) << endl;
		return 3;
	}

	for (size_t i {}; i < 2 * DataSize; ++i) isine[i] = fsine[i];
	vitis_rifft(isine, DataSize);

#ifdef FEXPORT
	ofs.open("rfft3.csv");
	ofs << "Sine,RFFT,IRFFT" << endl;
	for (size_t i {}; i < 2 * DataSize; ++i)
		ofs << (i < DataSize ? sine[i] : 0) << "," << fsine[i] << "," << (i < DataSize ? isine[i] : 0) << endl;
	ofs.close();
#endif

	// Check that the waveforms match the expectation
	diff = 0;
	for (size_t i {}; i < DataSize; ++i) diff += abs(sine[i] - isine[i]);
	if (diff > 1.0) {
		cout << "Waveform error " << diff << " exceeds the expected " << 1.0 << endl;
		return 3;
	}

	return 0;
}
