class COrthoWavelet00{
    public:
    COrthoWavelet00();
    ~COrthoWavelet00();

    //分解数列H0開始
    int H0Begin(){
        return -m_nPEnd;
    }

    //分解数列H0終了
    int H0End(){
        return -m_nPBegin;
    }

    //分解数列H1開始
    int H1Begin(){
        return -m_nQEnd;
    }

    //分解数列H1終了
    int H1End(){
        return -m_nQBegin;
    }

    //再構成数列H0開始
    int G0Begin(){
        return m_nPEnd;
    }

    //再構成数列H0終了
    int G0End(){
        return m_nPBegin;
    }

    //再構成数列H1開始
    int G1Begin(){
        return m_nQEnd;
    }

    //再構成数列H1終了
    int G1End(){
        return m_nQBegin;
    }

    bool Prepare(int N, int level);
    double Scaling(int level, int n);
    int ScalingBegin(int level);
    int ScalingEnd(int level);
    double Mother(int level, int n);
    int MotherBegin(int level);
    int MotherEnd(int level);
    double H0(int n);
    double H1(int n);
    double G0(int n);
    double G1(int n);

    protected:
    int m_nMtxMaxN;
    int m_nN;
    int m_nLevel;
    int m_nSLevel;
    int m_nMLevel;

    int m_nPBegin;
    int m_nPEnd;

    int m_nQBegin;
    int m_nQEnd;

    int m_nScalingBegin;
    int m_nScalingEnd;
    int m_nScalingDenomi;
    int m_nScalingSz;
    int m_nMotherBegin;
    int m_nMotherEnd;
    int m_nMotherDenomi;
    int m_nMotherSz;

    double* m_pBuf0;
    double* m_pBuf1;

    double* m_pScalingBuf;
    double* m_pMotherBuf;

    int m_nPNum;
    int m_nPShift;
    double* m_pP0;
    double* m_pQ0;
    double* m_pSc0;
    double* m_pPBuf0;
    double* m_pPBuf1;

    void PrepareSub(int N);


};