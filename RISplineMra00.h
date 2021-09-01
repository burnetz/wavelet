class CBSpline00;

class CRISplineMra00
{
public:

	CRISplineMra00();
	~CRISplineMra00();

protected:

	CBSpline00*		m_pRealSpline;
	CBSpline00*		m_pImagSpline;

	int				m_nMaxN;

	int				m_nLevel;
	int				m_nN;

	double*			m_pBuf;
	double*			m_pBetaR;
	double*			m_pBetaI;
	double*			m_pWork[4];

	int				m_nBetaRBegin;
	int				m_nBetaREnd;
	int				m_nBetaROffset;
	int				m_nBetaIBegin;
	int				m_nBetaIEnd;
	int				m_nBetaIOffset;

public:

	bool Prepare( int maxN );

	// 複素多重解像度解析
	bool Mra( double** output, double* input, int level, int mul, int lowMode );
	bool IMra( double* output, double** input, int level, int mul, int lowMode );

	// 実数部多重解像度解析
	bool RealMra( double* output, double* input, int level, int mul, int lowMode );
	bool ImagMra( double* output, double* input, int level, int mul, int lowMode );

	// 虚数部多重解像度解析
	bool RealIMra( double* output, double* input, int level, int mul, int lowMode );
	bool ImagIMra( double* output, double* input, int level, int mul, int lowMode );

protected:

	double BetaR(int n);
	double BetaI(int n);

	double Read( double* buf, int x, int sz, int ofst );
	void Write( double* buf, double d, int x, int sz, int ofst );
};
