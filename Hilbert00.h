class CFft00;

class CHilbert00
{
public:
    CHilbert00();
    ~CHilbert00();

    bool Prepare(int maxN);
    bool Hilbert(double *outputReal, double *outputImag,
                 double *input, int inputN, double freq0, double freq1);
    bool InstantPower(double* output, double* inputReal, double* inputImag, int inputN);
    bool InstantFreq(double* output, double* inputReal, double* inputImag, int inputN);

protected:
    CFft00* m_pFft;

    double* m_pWork[2];
    double* m_pReal;
    double* m_pImag;

    int m_nMaxEx;
    int m_nMaxN;
    
};

