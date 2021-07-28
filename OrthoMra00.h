class COrthoWavelet00;
class CInterpolation00;

class COrthoMra00{
    public:
    COrthoMra00();
    ~COrthoMra00();

    bool Prepare(int wName, int maxN);

    bool Mra(double* output, double* input, int level, int mul, int mode);
    bool IMra(double* output, double* input, int level, int mul, int mode);

    protected:
    double Read(int n, int lpsz, int offset);

    COrthoWavelet00* m_pOrtho;
    CInterpolation00* m_pInter;

    int m_nWName;
    int m_nMaxN;

    int m_nLevel;
    int m_nN;

    int m_nScalOffset;
    int m_nScalN;

    double* m_pBuf;
    double* m_pBuf0;
    double* m_pBuf1;
    double* m_pScal; 
};