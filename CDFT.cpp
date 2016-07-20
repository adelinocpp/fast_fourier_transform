#include "CDFT.h"
using namespace std;
using namespace math;
//--------------------------------------------------------------------------
bool CDFT::ComputeDFT(CWaveFile mWaveFile, bool ZeroPad, int idxWindow)
{
    ComputeDFT(mWaveFile.GetDataVector(), ZeroPad, idxWindow);
    return true;
};
//--------------------------------------------------------------------------
bool CDFT::ComputeDFT(vector<double> mdouble, bool ZeroPad, int idxWindow)
{
    int r1, r2, n0, n1, k0, k1, nPoints;
    int k, n, nPow2;
    matrix<double> RealX0, ImagX0, RealX1, ImagX1, RealX2, ImagX2;
    double Hreal, Himag, ArgConst;
    Complexdouble TempComp;

    vDFT.clear();
    nPoints = mdouble.size();
    ArgConst = 2*M_PI/(nPoints-1);

    if (idxWindow == 0)
        for (k = 0; k < nPoints; k++)
            vDFT.push_back(Complexdouble(0,0));
    else
    {
        for (k = 0; k < nPoints; k++)
        {
            mdouble[k] = mdouble[k]*(0.5*(1-cos(k*ArgConst)));
            vDFT.push_back(Complexdouble(0,0));
        }
    }


    if (ZeroPad)
    {
        nPow2 = ceil(log(nPoints)/log(2));
        nPow2 = pow(2,nPow2);
        for(k = nPoints; k < nPow2; k++)
        {
            mdouble.push_back(0);
            vDFT.push_back(Complexdouble(0,0));
        }
    }

    nPoints = mdouble.size();
    ArgConst = 2*M_PI/(nPoints-1);

    r1 = floor(sqrt(nPoints));
    r2 = ceil(sqrt(nPoints));
    while((r1*r2) != nPoints)
    {
        if ((r1*r2) >= nPoints)
            r2--;
        if ((r1*r2) <= nPoints)
            r1++;
        if ((r1*r2) == nPoints)
            break;
        if ((r1 == nPoints) || (r2 == 1))
        {
            r1 = nPoints;
            r2 = 1;
            break;
        }
    };
    RealX0.SetSize(r1,r2);
    ImagX0.SetSize(r1,r2);
    RealX1.SetSize(r1,r2);
    ImagX1.SetSize(r1,r2);
    RealX2.SetSize(r1,r2);
    ImagX2.SetSize(r1,r2);
    for(k0 = 0; k0 < r2; k0++)
    {
        for(k1 = 0; k1 < r1; k1++)
        {
            RealX0(k1,k0) = mdouble[k1*r2+k0];
            ImagX0(k1,k0) = 0.0;
        };
    };

    for(n0 = 0; n0 < r1; n0++)
    {
        for(k0 = 0; k0 < r2; k0++)
        {
            Hreal = 0;
            Himag = 0;
            ArgConst = (2*M_PI/nPoints)*n0*r2;
            for(k1 = 0; k1 < r1; k1++)
            {
                ArgConst = (2*M_PI/nPoints)*n0*r2*k1;
                Hreal = Hreal + RealX0(k1,k0)*cos(ArgConst) + ImagX0(k1,k0)*sin(ArgConst);
                Himag = Himag + ImagX0(k1,k0)*cos(ArgConst) - RealX0(k1,k0)*sin(ArgConst);
            }
            RealX1(n0,k0) = Hreal;
            ImagX1(n0,k0) = Himag;
        };
    };

    for(n0 = 0; n0 < r1; n0++)
    {
        for(n1 = 0; n1 < r2; n1++)
        {
            Hreal = 0;
            Himag = 0;
            for(k0 = 0; k0 < r2; k0++)
            {
                ArgConst = (2*M_PI/nPoints)*(n1*r1+n0)*k0;
                Hreal = Hreal + RealX1(n0,k0)*cos(ArgConst) + ImagX1(n0,k0)*sin(ArgConst);
                Himag = Himag + ImagX1(n0,k0)*cos(ArgConst) - RealX1(n0,k0)*sin(ArgConst);
            }
            RealX2(n0,n1) = Hreal;
            ImagX2(n0,n1) = Himag;
        };
    };
    for(n0 = 0; n0 < r1; n0++)
    {
        for(n1 = 0; n1 < r2; n1++)
        {
            TempComp = Complexdouble(RealX2(n0,n1),ImagX2(n0,n1));
            n = n1*r1 + n0;
            vDFT[n] = TempComp;
        };
    };
    IsComputed = true;
    return true;
};
//--------------------------------------------------------------------------
vector<double> CDFT::BuildFrequencyVector(double mFs)
{
    int i, nSamples;
    vector<double> tF;
    double ArgConst;

    nSamples = vDFT.size();
    ArgConst = mFs/nSamples;
    for(i = 0; i < nSamples; i++)
        tF.push_back(ArgConst*(i+0.5));

    return tF;
};
//--------------------------------------------------------------------------
vector<double> CDFT::GetMagnitudeDFT()
{
    int i, nSamples;
    vector<double> tMag;

    nSamples = vDFT.size();
    for(i = 0; i < nSamples; i++)
        tMag.push_back(abs(vDFT[i]));

    return tMag;
};
//--------------------------------------------------------------------------
vector<double> CDFT::GetLogMagnitudeDFT()
{
    int i, nSamples;
    vector<double> tMag;

    nSamples = vDFT.size();
    for(i = 0; i < nSamples; i++)
        tMag.push_back(log10(abs(vDFT[i])));

    return tMag;
};
//--------------------------------------------------------------------------
vector<double> CDFT::GetPhaseDFT()
{
    int i, nSamples;
    vector<double> tAng;

    nSamples = vDFT.size();
    for(i = 0; i < nSamples; i++)
        tAng.push_back(arg(vDFT[i]));

    return tAng;
};
//--------------------------------------------------------------------------
vector<double> CDFT::GetRealDFT()
{
    int i, nSamples;
    vector<double> tReal;

    nSamples = vDFT.size();
    for(i = 0; i < nSamples; i++)
        tReal.push_back(real(vDFT[i]));

    return tReal;
};
//--------------------------------------------------------------------------
vector<double> CDFT::GetImaginaryDFT()
{
    int i, nSamples;
    vector<double> tImag;

    nSamples = vDFT.size();
    for(i = 0; i < nSamples; i++)
        tImag.push_back(imag(vDFT[i]));

    return tImag;
};
//--------------------------------------------------------------------------
vector<double> CDFT::ComputeIDFT(vector<Complexdouble> mDFT)
{
    int r1, r2, n0, n1, k0, k1, nPoints;
    int k, n;
    matrix<double> RealX0, ImagX0, RealX1, ImagX1, RealX2, ImagX2;
    double Hreal, Himag, ArgConst;
    Complexdouble TempComp;
    vector<double> vVector;

    if(IsComputed == false)
        return vVector;

    nPoints = mDFT.size();

    for (k = 0; k < nPoints; k++)
        vVector.push_back(0.0);

    r1 = floor(sqrt(nPoints));
    r2 = ceil(sqrt(nPoints));
    while((r1*r2) != nPoints)
    {
        if ((r1*r2) >= nPoints)
            r2--;
        if ((r1*r2) <= nPoints)
            r1++;
        if ((r1*r2) == nPoints)
            break;
        if ((r1 == nPoints) || (r2 == 1))
        {
            r1 = nPoints;
            r2 = 1;
            break;
        }
    };
    RealX0.SetSize(r1,r2);
    ImagX0.SetSize(r1,r2);
    RealX1.SetSize(r1,r2);
    ImagX1.SetSize(r1,r2);
    RealX2.SetSize(r1,r2);
    ImagX2.SetSize(r1,r2);
    for(k0 = 0; k0 < r2; k0++)
    {
        for(k1 = 0; k1 < r1; k1++)
        {
            TempComp = mDFT[k1*r2+k0];
            RealX0(k1,k0) = TempComp.real();
            ImagX0(k1,k0) = TempComp.imag();
        };
    };

    for(n0 = 0; n0 < r1; n0++)
    {
        for(k0 = 0; k0 < r2; k0++)
        {
            Hreal = 0;
            Himag = 0;
            ArgConst = (2*M_PI/nPoints)*n0*r2;
            for(k1 = 0; k1 < r1; k1++)
            {
                ArgConst = (2*M_PI/nPoints)*n0*r2*k1;
                Hreal = Hreal + RealX0(k1,k0)*cos(ArgConst) - ImagX0(k1,k0)*sin(ArgConst);
                Himag = Himag + ImagX0(k1,k0)*cos(ArgConst) + RealX0(k1,k0)*sin(ArgConst);
            }
            RealX1(n0,k0) = Hreal;
            ImagX1(n0,k0) = Himag;
        };
    };

    for(n0 = 0; n0 < r1; n0++)
    {
        for(n1 = 0; n1 < r2; n1++)
        {
            Hreal = 0;
            Himag = 0;
            for(k0 = 0; k0 < r2; k0++)
            {
                ArgConst = (2*M_PI/nPoints)*(n1*r1+n0)*k0;
                Hreal = Hreal + RealX1(n0,k0)*cos(ArgConst) - ImagX1(n0,k0)*sin(ArgConst);
                Himag = Himag + ImagX1(n0,k0)*cos(ArgConst) + RealX1(n0,k0)*sin(ArgConst);
            }
            RealX2(n0,n1) = Hreal;
            ImagX2(n0,n1) = Himag;
        };
    };
    //ArgConst = 2*M_PI/(nPoints-1);
    for(n0 = 0; n0 < r1; n0++)
    {
        for(n1 = 0; n1 < r2; n1++)
        {
            TempComp = Complexdouble(RealX2(n0,n1),ImagX2(n0,n1));
            n = n1*r1 + n0;
            vVector[n] = TempComp.real()/nPoints;
            //vVector[n] = TempComp.real()/(nPoints*);
        };
    };
    return vVector;
};
//--------------------------------------------------------------------------
vector<double> CDFT::GetAutoCorrelation()
{
    int i, nSamples;
    Complexdouble tempComplex;
    vector<Complexdouble> ConjComplex;
    vector<double> tempAutoCorr;
    if(IsComputed == false)
        return tempAutoCorr;
    nSamples = vDFT.size();
    for(i = 0; i < nSamples; i++)
    {
        tempComplex = vDFT[i];
        ConjComplex.push_back(Complexdouble((pow(tempComplex.real(),2) + pow(tempComplex.imag(),2)),0.0));
    }
    return this->ComputeIDFT(ConjComplex);
};
//--------------------------------------------------------------------------
vector<double> CDFT::ComputeIDFT()
{
    return this->ComputeIDFT(vDFT);
};
//--------------------------------------------------------------------------
vector<double> CDFT::ComputeIDFT(vector<Complexdouble> mDFT, double fIni, double fFim)
{
    int i, nSamples;
    bool bTestF1, bTestF2, bTestF3;
    if(fIni < 0)
        fIni = 0;
    if(fFim > 1.0)
        fFim = 1.0;
    nSamples = mDFT.size();
    for(i = 0; i < nSamples; i++)
    {
        bTestF1 = (i <= floor((fIni/2)*(nSamples-1)) );
        bTestF2 = ( (i >= ceil((fFim/2)*(nSamples-1))) & (i <= floor((1 - fFim/2)*(nSamples-1))) );
        bTestF3 = (i >= ceil((1 - fFim/2)*(nSamples-1)/2));

        if (bTestF1 | bTestF2 | bTestF3)
            mDFT[i] = Complexdouble(0.0,0.0);
    }
    return this->ComputeIDFT(mDFT);
};
//--------------------------------------------------------------------------
vector<double> CDFT::ComputeIDFT(double fIni, double fFim)
{
    vector<Complexdouble> mDFT;
    mDFT = vDFT;
    return this->ComputeIDFT(mDFT, fIni, fFim);
};
//--------------------------------------------------------------------------
vector<Complexdouble> CDFT::ComputeComplexIDFT(vector<Complexdouble> mDFT)
{
    int r1, r2, n0, n1, k0, k1, nPoints;
    int k, n;
    matrix<double> RealX0, ImagX0, RealX1, ImagX1, RealX2, ImagX2;
    double Hreal, Himag, ArgConst;
    Complexdouble TempComp;
    vector<Complexdouble> vVector;

    if(IsComputed == false)
        return vVector;

    nPoints = mDFT.size();

    for (k = 0; k < nPoints; k++)
        vVector.push_back(0.0);

    r1 = floor(sqrt(nPoints));
    r2 = ceil(sqrt(nPoints));
    while((r1*r2) != nPoints)
    {
        if ((r1*r2) >= nPoints)
            r2--;
        if ((r1*r2) <= nPoints)
            r1++;
        if ((r1*r2) == nPoints)
            break;
        if ((r1 == nPoints) || (r2 == 1))
        {
            r1 = nPoints;
            r2 = 1;
            break;
        }
    };
    RealX0.SetSize(r1,r2);
    ImagX0.SetSize(r1,r2);
    RealX1.SetSize(r1,r2);
    ImagX1.SetSize(r1,r2);
    RealX2.SetSize(r1,r2);
    ImagX2.SetSize(r1,r2);
    for(k0 = 0; k0 < r2; k0++)
    {
        for(k1 = 0; k1 < r1; k1++)
        {
            TempComp = mDFT[k1*r2+k0];
            RealX0(k1,k0) = TempComp.real();
            ImagX0(k1,k0) = TempComp.imag();
        };
    };

    for(n0 = 0; n0 < r1; n0++)
    {
        for(k0 = 0; k0 < r2; k0++)
        {
            Hreal = 0;
            Himag = 0;
            ArgConst = (2*M_PI/nPoints)*n0*r2;
            for(k1 = 0; k1 < r1; k1++)
            {
                ArgConst = (2*M_PI/nPoints)*n0*r2*k1;
                Hreal = Hreal + RealX0(k1,k0)*cos(ArgConst) - ImagX0(k1,k0)*sin(ArgConst);
                Himag = Himag + ImagX0(k1,k0)*cos(ArgConst) + RealX0(k1,k0)*sin(ArgConst);
            }
            RealX1(n0,k0) = Hreal;
            ImagX1(n0,k0) = Himag;
        };
    };

    for(n0 = 0; n0 < r1; n0++)
    {
        for(n1 = 0; n1 < r2; n1++)
        {
            Hreal = 0;
            Himag = 0;
            for(k0 = 0; k0 < r2; k0++)
            {
                ArgConst = (2*M_PI/nPoints)*(n1*r1+n0)*k0;
                Hreal = Hreal + RealX1(n0,k0)*cos(ArgConst) - ImagX1(n0,k0)*sin(ArgConst);
                Himag = Himag + ImagX1(n0,k0)*cos(ArgConst) + RealX1(n0,k0)*sin(ArgConst);
            }
            RealX2(n0,n1) = Hreal;
            ImagX2(n0,n1) = Himag;
        };
    };
    ArgConst = 2*M_PI/(nPoints-1);
    for(n0 = 0; n0 < r1; n0++)
    {
        for(n1 = 0; n1 < r2; n1++)
        {
            TempComp = Complexdouble(RealX2(n0,n1),ImagX2(n0,n1));
            n = n1*r1 + n0;
            vVector[n] = Complexdouble(TempComp.real()/nPoints,TempComp.imag()/nPoints);
        };
    };
    return vVector;
};
//--------------------------------------------------------------------------
