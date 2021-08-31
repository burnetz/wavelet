class CBSpline00
{
public:

	CBSpline00();
	~CBSpline00();

	int H0Begin()
		{ return m_nH0Begin; }

	int H0End()
		{ return m_nH0End; }

	int H1Begin()
		{ return m_nH1Begin; }

	int H1End()
		{ return m_nH1End; }

	int G0Begin()
		{ return m_nG0Begin; }

	int G0End()
		{ return m_nG0End; }

	int G1Begin()
		{ return m_nG1Begin; }

	int G1End()
		{ return m_nG1End; }

	int Sc1Begin()
		{ return m_nSc1Begin; }

	int Sc1End()
		{ return m_nSc1End; }

	int Sc2Begin()
		{ return m_nSc2Begin; }

	int Sc2End()
		{ return m_nSc2End; }

	int BetaSc2Begin()
		{ return m_nBetaSc2Begin; }

	int BetaSc2End()
		{ return m_nBetaSc2End; }

	bool		Prepare( int m );
	double		H0( int n );	// 分解数列
	double		H1( int n );	// 分解数列
	double		G0( int n );	// 再構成数列
	double		G1( int n );	// 再構成数列
	double		Sc1( int n );	// オリジナルのΦ(x)
	double		Sc2( int n );	// 原点を頂点とするΦ(x)
	double		BetaSc2( int n );	// 原点を頂点とするΦ(x)による補間

protected:

	int			m_nH0Begin;
	int			m_nH0End;
	int			m_nH0Offset;

	int			m_nH1Begin;
	int			m_nH1End;
	int			m_nH1Offset;

	int			m_nG0Begin;
	int			m_nG0End;
	int			m_nG0Offset;

	int			m_nG1Begin;
	int			m_nG1End;
	int			m_nG1Offset;

	int			m_nSc1Begin;
	int			m_nSc1End;
	int			m_nSc1Offset;

	int			m_nSc2Begin;
	int			m_nSc2End;
	int			m_nSc2Offset;

	int			m_nBetaSc2Begin;
	int			m_nBetaSc2End;
	int			m_nBetaSc2Offset;

	double*		m_pBuf;

	double*		m_pH0;
	double*		m_pH1;
	double*		m_pG0;
	double*		m_pG1;
	double*		m_pSc1;
	double*		m_pSc2;
	double*		m_pBetaSc2;
};
