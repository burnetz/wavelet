class CGraphMatchingPursuit00
{
public:

	CGraphMatchingPursuit00();
	~CGraphMatchingPursuit00();

	bool Pursuit( double* output, double* input, int nX, int nY,
	int wName, int level, double alpha, double epsilon, double rho );

	bool IPursuit( double* output, double* input, int nX, int nY,
	int wName, int level );

protected:

	double*	m_pBuf;
};
