class COrthoWavelet00;
class COrthoMra00;
class CInterpolation00;

class CMatchingPursuit00
{
public:
	CMatchingPursuit00();
	~CMatchingPursuit00();

	bool Prepare(int wName, int maxN, int maxLevelt);

	bool Pursuit(double *output[], double *input,
					int inputN, int level, double alpha, double epsilon,
					double rho, int mode);

	bool IPursuit(double *output, double *input[],
				int inputN, int level, int mode);

	double Epsilon();

protected:
	COrthoWavelet00 *m_pOrtho;
	COrthoMra00 *m_pMra;
	CInterpolation00 *m_pInter;

	double ReadCZero(int x);
	double ReadB(int level, int x);
	void UpdateBuf(int L0, int x0);

	int m_nMaxN;
	int m_nMaxLevel;

	int m_nN;
	int m_nLevel;
	int m_nCount;
	int m_nBOffset;
	int m_nUOffset;

	double m_nEpsilon;

	double *m_pBuf;
	double *m_pDictBuf[20];
	double *m_pBBuf[20];
	int m_nBBegin[20];
	int m_nBEnd[20];
	double *m_pUBuf[10][10];
	int m_nUBegin[10][10];
	int m_nUEnd[10][10];
	double *m_pScal;
	int m_nScalN;
	int m_nScalOffset;

	double *m_pLowPassBuf;
};
