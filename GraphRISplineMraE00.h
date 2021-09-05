class CRISplineMra00;

class CGraphRISplineMraE00
{
public:

	CGraphRISplineMraE00();
	~CGraphRISplineMraE00();

	bool Prepare( int maxX, int maxY );

	bool Mra(double* output, double* input, int nX, int nY, int level);
	bool IMra(double* output, double* input, int nX, int nY, int level);

protected:

	void SelectMra( double* output, double* input, int n, int ri );
	void SelectIMra( double* output, double* input, int n, int ri );

	CRISplineMra00* m_pMra;

	int m_nMaxX;
	int m_nMaxY;
	int m_nMaxLen;
	int m_nMaxN;
	int m_nX;
	int m_nY;
	int m_nLevel;

	double*	m_pBuf[3];
	double* m_pBox0[2][2];
	double* m_pBox1[2][2];
	double* m_pBox2[2][2];
	double* m_pWork;
};
