class COrthoMra00;

class CGraphOrthoMra00{
    public:
    CGraphOrthoMra00();
    ~CGraphOrthoMra00();

    bool Prepare(int wName, int maxX, int maxY);
    bool Mra(double* output, double* input, int nX, int nY, int level);
    bool IMra(double* output, double* input, int nX, int nY, int level, double delta);

    protected:
    COrthoMra00* m_pMra;

    int m_nWName;
    int m_nMaxX;
    int m_nMaxY;
    int m_nMaxLen;
    int m_nMaxN;
    int m_nX;
    int m_nY;
    int m_nLevel;
    double* m_pBuf[3];
    double* m_pWork;
};