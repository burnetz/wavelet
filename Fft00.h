class CFft00{
    public:
    CFft00();
    ~CFft00();

    bool Prepare(int ex);
    bool Fft(double* real, double* imag, int inv);

    protected:
    void BitReverce(double* buf, double* a, int ex);

    int m_nEx;
    int m_nLength;
    int m_nBufLen;
    double* m_pBuf;
    double* m_pFftSin;
    double* m_pFftCos;
    double* m_pIfftSin;
    double* m_pIfftCos;
    double* m_pSubBuf;
};