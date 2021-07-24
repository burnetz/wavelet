class CGabor00
{
public:
    CGabor00();
    ~CGabor00();

    bool Prepare(double *input, int inputN);
    bool SetParameter(double sigma, double omega, double a);
    double Gabor(int b);

protected:
    double *m_pWork[8];
    double *m_pInput;
    int m_nInputN;

    double m_nSigma;
    double m_nOmega;
    double m_nA;

    double *m_pGaborR;
    double *m_pGaborI;
    int m_nGaborOffset;
    int m_nGaborN;

    double ReadInput(int x);
    double ReadGaborR(int x);
    double ReadGaborI(int x);
    void ComplexMultiply(double &rOutput, double &iOutput, double rInput0, double iInput0, double rInput1, double iInput1);
};