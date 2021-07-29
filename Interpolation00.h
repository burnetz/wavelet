class CFft00;

class CInterpolation00{
    public:
    CInterpolation00();
    ~CInterpolation00();

    int m_nMaxVol;
    int m_nMaxEx;
    int m_nOrthoN;

    int m_nMLevel00;
    int m_nILevel00;
    
    bool Prepare();
    bool Inter(double* output, double* input, int inputN, double* scal, int scalN, int scalOffset, int mode);
    bool IInter(double* output, double* input, int inputN, double* scal, int scalN, int scalOffset, int mode);

    protected:
    CFft00* m_pFft;

    double* m_pBuf;
    double* m_pBuf0;
    double* m_pBuf1;
    double* m_pBuf2;
    double* m_pBuf3;
    double* m_pBuf4;
    double* m_pBuf5;
    double* m_pBuf6;
};