#pragma once


typedef float SAMPLETYPE;
#define M_PI 3.14159265358979323846
#define win_len 1280
#define MAX_FRAME_LENGTH 8192
#define TRUE 1
#define FALSE 0
struct ZykeeFourier {

    SAMPLETYPE gInFIFO[MAX_FRAME_LENGTH];
    SAMPLETYPE gOutFIFO[MAX_FRAME_LENGTH];
    SAMPLETYPE gLastPhase[MAX_FRAME_LENGTH / 2 + 1];
    SAMPLETYPE gSumPhase[MAX_FRAME_LENGTH / 2 + 1];
    SAMPLETYPE gOutputAccum[2 * MAX_FRAME_LENGTH];
    SAMPLETYPE gAnaFreq[MAX_FRAME_LENGTH];
    SAMPLETYPE gAnaMagn[MAX_FRAME_LENGTH];
    SAMPLETYPE gSynFreq[MAX_FRAME_LENGTH];
    SAMPLETYPE gSynMagn[MAX_FRAME_LENGTH];
    SAMPLETYPE win[win_len];
    long gRover;
    long gInit;
    double tempo;
    double rate;
    int sampleRate;
};
void ZykeePitchShiftInit(struct ZykeeFourier* Zy, double pitch);
void ZykeePitchShift(double pitchShift, long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, SAMPLETYPE* indata, SAMPLETYPE* outdata, struct ZykeeFourier* Zy);