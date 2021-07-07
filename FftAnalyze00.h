class CFft00;

class CFftAnalyze00{
    public:
    CFftAnalyze00();
    ~CFftAnalyze00();

    bool Prepare(int ex, int winType);
    bool Fft(double* outAmp, double* outPhase, double* input);
    protected:
    void Window(double* output, int winLen, int winType, double& scale);

    CFft00* m_pFft;

    double* m_pWork[16];
    double* m_pWindow;
    double* m_pReal;
    double* m_pImag;

    int m_nEx;
    int m_nFftLen;
    double m_nScale;
};