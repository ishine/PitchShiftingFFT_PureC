#include "stdint.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "ZykeePitchShiftFFT.h"
#include "assert.h"
#include "fftw3.h"


void ZykeePitchShiftInit(struct ZykeeFourier* Zy,double pitch)
{
    memset(Zy->gInFIFO, 0, MAX_FRAME_LENGTH * sizeof(SAMPLETYPE));
    memset(Zy->gOutFIFO, 0, MAX_FRAME_LENGTH * sizeof(SAMPLETYPE));
    memset(Zy->gLastPhase, 0, (MAX_FRAME_LENGTH / 2 + 1) * sizeof(SAMPLETYPE));
    memset(Zy->gSumPhase, 0, (MAX_FRAME_LENGTH / 2 + 1) * sizeof(SAMPLETYPE));
    memset(Zy->gOutputAccum, 0, 2 * MAX_FRAME_LENGTH * sizeof(SAMPLETYPE));
    memset(Zy->gAnaFreq, 0, MAX_FRAME_LENGTH * sizeof(SAMPLETYPE));
    memset(Zy->gAnaMagn, 0, MAX_FRAME_LENGTH * sizeof(SAMPLETYPE));
    for (int k = 0; k < win_len; k++) {
        Zy->win[k] = -.5 * cos(2. * M_PI * (double)k / (double)win_len) + .5;
    }
    Zy->gRover = FALSE;
    Zy->gInit = TRUE;
    Zy->rate = exp(0.69314718056 * (pitch / 12.0));
    Zy->tempo = 1 / Zy->rate;
    Zy->sampleRate = 44100;
}

void ZykeePitchShift(double pitchShift, long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, SAMPLETYPE* indata, SAMPLETYPE* outdata, struct ZykeeFourier* Zy)
{
    double magn, phase, tmp, window, real, imag;
    double freqPerBin, expct;
    long i, k, qpd, index, inFifoLatency, stepSize, fftFrameSize2;
    fftw_complex* fft_in, * fft_out;
    fftw_plan fft, ifft;
    fft_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * win_len);
    fft_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * win_len);
    fft = fftw_plan_dft_1d(win_len, fft_in, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    ifft = fftw_plan_dft_1d(win_len, fft_in, fft_out, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftFrameSize2 = fftFrameSize / 2;
    stepSize = fftFrameSize / osamp;
    freqPerBin = sampleRate / (double)fftFrameSize;
    expct = 2. * M_PI * (double)stepSize / (double)fftFrameSize;
    inFifoLatency = fftFrameSize - stepSize;
    if (Zy->gRover == FALSE)
        Zy->gRover = inFifoLatency;
    /* main processing loop */
    for (i = 0; i < numSampsToProcess; i++) {
        Zy->gInFIFO[Zy->gRover] = indata[i];
        outdata[i] = Zy->gOutFIFO[Zy->gRover - inFifoLatency];
        Zy->gRover++;
        /* now we have enough data for processing */
        if (Zy->gRover >= fftFrameSize) {
            Zy->gRover = inFifoLatency;
            for (k = 0; k < fftFrameSize; k++) {
                fft_in[k][0] = Zy->gInFIFO[k] * Zy->win[k];
                fft_in[k][1] = 0;
            }
            /* ***************** ANALYSIS ******************* */
            /* do transform */
            fftw_execute(fft);
            /* this is the analysis step */
            for (k = 0; k <= fftFrameSize2; k++) {
                real = fft_out[k][0];
                imag = fft_out[k][1];
                magn = 2. * sqrt(real * real + imag * imag);
                phase = atan2(imag, real);
                tmp = phase - Zy->gLastPhase[k];
                Zy->gLastPhase[k] = phase;
                tmp -= (double)k * expct;
                qpd = tmp / M_PI;
                if (qpd >= 0) qpd += qpd & 1;
                else qpd -= qpd & 1;
                tmp -= M_PI * (double)qpd;
                tmp = osamp * tmp / (2. * M_PI);
                tmp = (double)k * freqPerBin + tmp * freqPerBin;
                Zy->gAnaMagn[k] = magn;
                Zy->gAnaFreq[k] = tmp;

            }
            /* ***************** PROCESSING ******************* */
            /* this does the actual pitch shifting */
            memset(Zy->gSynMagn, 0, fftFrameSize * sizeof(float));
            memset(Zy->gSynFreq, 0, fftFrameSize * sizeof(float));
            for (k = 0; k <= fftFrameSize2; k++) {
                index = k * pitchShift;
                if (index <= fftFrameSize2) {
                    Zy->gSynMagn[index] += Zy->gAnaMagn[k];
                    Zy->gSynFreq[index] = Zy->gAnaFreq[k] * pitchShift;
                }
            }
            /* ***************** SYNTHESIS ******************* */
            /* this is the synthesis step */
            for (k = 0; k <= fftFrameSize2; k++) {
                magn = Zy->gSynMagn[k];
                tmp = Zy->gSynFreq[k];
                tmp -= (double)k * freqPerBin;
                tmp /= freqPerBin;
                tmp = 2. * M_PI * tmp / osamp;
                tmp += (double)k * expct;
                Zy->gSumPhase[k] += tmp;
                phase = Zy->gSumPhase[k];
                fft_in[k][0] = magn * cos(phase);
                fft_in[k][1] = magn * sin(phase);
            }
            for (k = fftFrameSize2 + 1; k < fftFrameSize; k++)
            {
                fft_in[k][0] = 0.;
                fft_in[k][1] = 0.;
            }
            fftw_execute(ifft);
            for (k = 0; k < fftFrameSize; k++)
            {
                Zy->gOutputAccum[k] += 2. * Zy->win[k] * fft_out[k][0] / (fftFrameSize2 * osamp);
            }
            memcpy(Zy->gOutFIFO, Zy->gOutputAccum, stepSize * sizeof(float));
            memmove(Zy->gOutputAccum, Zy->gOutputAccum + stepSize, fftFrameSize * sizeof(float));
            memcpy(Zy->gInFIFO, Zy->gInFIFO + stepSize, sizeof(float) * inFifoLatency);
        }
    }
}