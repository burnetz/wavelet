class CDwvd00
{
public:
    CDwvd00();
    ~CDwvd00();

    bool Prepare(double *input, int inputN);
    bool SetParameters(double sigma, double damp, double omega);
    double Dwvd(int x);
    bool PrepareSp(double sigma, double damp);
    double Sp(double *input, int inputN, int x);

protected:
    double* m_pWork[8];
    double* m_pReal;
    double* m_pImag;
    double* m_pWindowR;
    double* m_pWindowI;

    double* m_pGauss;
    double* m_pInput;
    int m_nInputN;

    int m_nWinLength;
    int m_nWinOffset;

    double* m_pSpInput;
    int m_nSpInputN;
    int m_nSpGaussLength;
    int m_nSpGaussOffset;

    double Real(int x);
    double Imag(int x);
    double ReadWindowR(int x);
    double ReadWindowI(int x);
    void ComplexMultiply(double &rOutput, double &iOutput, double rInput0, double iInput0, double rInput1, double iInput1);
    double SpReadData(int x);
    double SpReadGauss(int x);
};