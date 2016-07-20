#ifndef CDFT_H
#define CDFT_H

#include "CWaveFile.h"
//--------------------------------------------------------------------------
//   Interface da classe CDFT
//   Abstração para calculo da Transformada discreta de fourier
//--------------------------------------------------------------------------
class CDFT
{
private:
    vector<Complexdouble> vDFT;
    bool IsComputed;
public:
    CDFT()
    {
        vDFT.clear();
    };
    ~CDFT() {};
    bool ComputeDFT(CWaveFile mWaveFile, bool ZeroPad, int idxWindow);
    bool ComputeDFT(vector<double> mdouble, bool ZeroPad, int idxWindow);
    vector<Complexdouble> GetDFT()
    {
        return vDFT;
    };
    vector<double> BuildFrequencyVector(double mFs);
    vector<double> GetMagnitudeDFT();
    vector<double> GetLogMagnitudeDFT();
    vector<double> GetPhaseDFT();
    vector<double> GetRealDFT();
    vector<double> GetImaginaryDFT();
    vector<double> GetAutoCorrelation();
    vector<Complexdouble> ComputeComplexIDFT(vector<Complexdouble> mDFT);
    vector<double> ComputeIDFT(vector<Complexdouble> mDFT);
    vector<double> ComputeIDFT(vector<Complexdouble> mDFT, double fIni, double fFim);
    vector<double> ComputeIDFT();
    vector<double> ComputeIDFT(double fIni, double fFim);
};

#endif // CDFT_H
