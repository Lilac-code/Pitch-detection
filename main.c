#include <stdio.h>

// Pi (M_PI) is defined in this header file
#include <math.h>

// The library to read and write wav files
#include "tinywav.h"

// The library for fourier transform
#include "kiss_fftr.h"

// The output of the FFT method is an array of complex numbers
#include <complex.h>

float equation6(int tau,
                int windowLength,
                float* signal) {
    float result = 0;
    for (int index=0; index<windowLength; index++) {
        result = result + pow(signal[index] - signal[index+tau], 2);
    }
    return result;
}

float equation8(int tau,
                int windowLength,
                float* signal,
                float* sumSqrDiffs) {
    if (tau == 0) {
        return 1.0;
    }
    else {
        sumSqrDiffs[tau] = equation6(tau, windowLength, signal);

        float mean = 0;
        for (int index=1; index<=tau; index++) {
            mean = mean + sumSqrDiffs[index];
        }
        mean = mean / tau;

        return sumSqrDiffs[tau] / mean;
    }
}

float yin(float* signal, float samplingRate) {
    float minF0 = 50; // Hz                 // the lower the frequency, the longer it takes to calculate, Yin chose 30hs in his code
    float maxF0 = floor(samplingRate/4); // divide by 4 is taken from the YIN MATLAB code
    float threshold = 0.1;               // for step 4
    int windowLength = ceil(samplingRate/minF0);   // taken from YIN MATLAB code

    printf("windowLength %d", windowLength);

    float yins[windowLength];   //this will contain the amps in fig 3
    float sumSqrDiffs[windowLength]; // this is an array that memorizes the values of equation 6 for previous taus
    for (int tau=0; tau<windowLength; tau++) {
        yins[tau] = equation8(tau, windowLength, signal, sumSqrDiffs);

        // to store in a csv file
        printf("%d, %lf\n", tau, yins[tau]);
    }

    //fig 3 (b)

    // Find the minimum dip
    int periodInNumSamples = 0; // the tau of the selected dip
    float globalMinimum = 1.0;
    for (int tau=0; tau<windowLength-1; tau++) {
        // Step 4 in the YIN algorithm
        if (yins[tau] < threshold) {
            // Local minimum: lower than the values immediately before and after
            if ((yins[tau] <= yins[tau-1]) &&
                (yins[tau] <= yins[tau+1])) {
                periodInNumSamples = tau;
                break;
            }
        }
        if (yins[tau] < globalMinimum) {
            globalMinimum = yins[tau];
            periodInNumSamples = tau;
        }
    }

    // convert period to frequency
    float f0 = samplingRate / periodInNumSamples;;
    return f0;
}

int main() {
    ////////////////////
    /// Generate the tone
    ////////////////////
    int samplingRate = 44100; // Hertz
    float frequency = 500.0; // Hertz
    float timeLength = 0.5; // Seconds
    float amplitude = 1.0;

    // The number of samples the code generates for the signal
    int numberOfSamples = timeLength * samplingRate;

    float radiansPerSecond = 2 * M_PI * frequency;

    float signal[numberOfSamples];

    for (int i=0; i<numberOfSamples; i++) {
        signal[i] = amplitude * sin(radiansPerSecond * i / samplingRate);
        float sampleTime = (float)i / samplingRate;
//        printf("%lf, %lf \n", sampleTime, signal[i]);
    }

    ////////////////////
    /// Store the signal in a wav file
    ////////////////////
    TinyWav tw;
    tinywav_open_write(&tw,
                       1,
                       samplingRate,
                       TW_FLOAT32, // the output samples will be 32-bit floats. TW_INT16 is also supported
                       TW_INLINE,  // the samples will be presented inlined in a single buffer.
                       "pitch.wav" // the output path
    );
    tinywav_write_f(&tw, signal, numberOfSamples);

    tinywav_close_write(&tw);

    ////////////////////
    /// YIN
    ////////////////////
    float yinF0 = yin(signal, samplingRate);
    printf("YIN F0 estimation: %lf\n", yinF0);

    ////////////////////
    /// Auto-correlation
    ////////////////////
    int tau = 0;
    float lowestFrequency = 100; // Hz
    float longestPeriodInSeconds = 1.0 / lowestFrequency;
    float longestPeriodInSamples = longestPeriodInSeconds * samplingRate;

    for (tau=1; tau<longestPeriodInSamples; tau++) {
        float autocorr = 0;
        for (int index=0; index<numberOfSamples-tau; index++) {
            autocorr = autocorr + signal[index] * signal[index+tau];
        }
        // normalize the dot product
        autocorr = autocorr / (numberOfSamples - tau);
        // calculate the  corresponding frequency to tau
        float frequency = samplingRate / tau;
//        printf("%lf, %lf \n", frequency, autocorr);
    }

    ////////////////////
    /// Spectrum
    ////////////////////
    int nfft = numberOfSamples;
    // FFT of a real signal
    kiss_fftr_cfg cfg = kiss_fftr_alloc(nfft, 0, 0, 0);
    // The array to store the FFT's output
    float complex ftOutput[nfft/2+1];
    // Compute the Fourier Transform of the signal
    kiss_fftr(cfg, signal, (kiss_fft_cpx*)ftOutput);

    for (int bin=0; bin<nfft/2+1; bin++) {
        float magnitude = cabsf(ftOutput[bin]) / nfft * 2.0;
        float binFrequency = bin / (float) (nfft / 2) * (samplingRate / 2);
//        printf("%0.2lf, %0.2lf \n", binFrequency, magnitude);
    }

    kiss_fft_free(cfg);

    return 0;
}