#include <stdio.h>
#include <math.h>
#include "RISplineMra00.h"
#include "BSpline00.h"

#define MAXLEVEL 16



CRISplineMra00::CRISplineMra00()
{
	m_pRealSpline = new CBSpline00;
	m_pImagSpline = new CBSpline00;

	m_pBuf = NULL;
}

CRISplineMra00::~CRISplineMra00()
{
	if( m_pRealSpline != NULL )
		delete m_pRealSpline;

	if( m_pImagSpline != NULL )
		delete m_pImagSpline;

	if( m_pBuf != NULL )
		delete m_pBuf;
}

bool CRISplineMra00::Prepare( int maxN )
{
	int i;

	m_nMaxN = maxN;

	if( m_pBuf != NULL )
		delete m_pBuf;
	m_pBuf = new double[ 256*2 + m_nMaxN*4 + 256 ];

	if( m_pBuf != NULL )
	{
		m_pBetaR = m_pBuf;			// 補間の係数のバッファ
		m_pBetaI = m_pBetaR + 256;

		m_pWork[0] = m_pBetaI + 256;	// ワークエリア
		m_pWork[1] = m_pWork[0] + m_nMaxN;
		m_pWork[2] = m_pWork[1] + m_nMaxN;
		m_pWork[3] = m_pWork[2] + m_nMaxN;
	}
	else
	{
		return false;
	}

	// 実数部には4階Splineウェーブレット
	if( m_pRealSpline->Prepare( 4 ) == false )
		return false;

	// 虚数部には3階Splineウェーブレット
	if( m_pImagSpline->Prepare( 3 ) == false )
		return false;

	// 実数部の補間の数列βをバッファに設定
	m_nBetaRBegin = -31;
	m_nBetaREnd = 31;
	m_nBetaROffset = 128;
	m_pBetaR[m_nBetaROffset +  0] =  0.753107867808;
	m_pBetaR[m_nBetaROffset +  1] = -0.298756981267;
	m_pBetaR[m_nBetaROffset +  2] =  0.147442214619;
	m_pBetaR[m_nBetaROffset +  3] = -0.083695103622;
	m_pBetaR[m_nBetaROffset +  4] =  0.051488596200;
	m_pBetaR[m_nBetaROffset +  5] = -0.032375989433;
	m_pBetaR[m_nBetaROffset +  6] =  0.019730500034;
	m_pBetaR[m_nBetaROffset +  7] = -0.011583529427;
	m_pBetaR[m_nBetaROffset +  8] =  0.006830545907;
	m_pBetaR[m_nBetaROffset +  9] = -0.004063234733;
	m_pBetaR[m_nBetaROffset + 10] =  0.002428424925;
	m_pBetaR[m_nBetaROffset + 11] = -0.001451381897;
	m_pBetaR[m_nBetaROffset + 12] =  0.000865696123;
	m_pBetaR[m_nBetaROffset + 13] = -0.000515751651;
	m_pBetaR[m_nBetaROffset + 14] =  0.000307321675;
	m_pBetaR[m_nBetaROffset + 15] = -0.000183205927;
	m_pBetaR[m_nBetaROffset + 16] =  0.000109238414;
	m_pBetaR[m_nBetaROffset + 17] = -0.000065132096;
	m_pBetaR[m_nBetaROffset + 18] =  0.000038830353;
	m_pBetaR[m_nBetaROffset + 19] = -0.000023148825;
	m_pBetaR[m_nBetaROffset + 20] =  0.000013800394;
	m_pBetaR[m_nBetaROffset + 21] = -0.000008227391;
	m_pBetaR[m_nBetaROffset + 22] =  0.000004904930;
	m_pBetaR[m_nBetaROffset + 23] = -0.000002924101;
	m_pBetaR[m_nBetaROffset + 24] =  0.000001743100;
	m_pBetaR[m_nBetaROffset + 25] = -0.000001038896;
	m_pBetaR[m_nBetaROffset + 26] =  0.000000618862;
	m_pBetaR[m_nBetaROffset + 27] = -0.000000368140;
	m_pBetaR[m_nBetaROffset + 28] =  0.000000218203;
	m_pBetaR[m_nBetaROffset + 29] = -0.000000128005;
	m_pBetaR[m_nBetaROffset + 30] =  0.000000072549;
	m_pBetaR[m_nBetaROffset + 31] = -0.000000035311;
	for( i=1; i<32; i++ )
		m_pBetaR[m_nBetaROffset - i] = m_pBetaR[m_nBetaROffset + i];

	// 虚数部の補間の数列βをバッファに設定
	m_nBetaIBegin = -37;
	m_nBetaIEnd = 36;
	m_nBetaIOffset = 128;
	m_pBetaI[m_nBetaIOffset +  0] =  0.310498828942;
	m_pBetaI[m_nBetaIOffset +  1] = -0.205356034324;
	m_pBetaI[m_nBetaIOffset +  2] =  0.131333945767;
	m_pBetaI[m_nBetaIOffset +  3] = -0.082829084249;
	m_pBetaI[m_nBetaIOffset +  4] =  0.050736413374;
	m_pBetaI[m_nBetaIOffset +  5] = -0.029925906319;
	m_pBetaI[m_nBetaIOffset +  6] =  0.017442614011;
	m_pBetaI[m_nBetaIOffset +  7] = -0.010382673592;
	m_pBetaI[m_nBetaIOffset +  8] =  0.006213985653;
	m_pBetaI[m_nBetaIOffset +  9] = -0.003716801670;
	m_pBetaI[m_nBetaIOffset + 10] =  0.002217526462;
	m_pBetaI[m_nBetaIOffset + 11] = -0.001320833184;
	m_pBetaI[m_nBetaIOffset + 12] =  0.000786813652;
	m_pBetaI[m_nBetaIOffset + 13] = -0.000469044718;
	m_pBetaI[m_nBetaIOffset + 14] =  0.000279691640;
	m_pBetaI[m_nBetaIOffset + 15] = -0.000166770272;
	m_pBetaI[m_nBetaIOffset + 16] =  0.000099425392;
	m_pBetaI[m_nBetaIOffset + 17] = -0.000059271885;
	m_pBetaI[m_nBetaIOffset + 18] =  0.000035335196;
	m_pBetaI[m_nBetaIOffset + 19] = -0.000021065923;
	m_pBetaI[m_nBetaIOffset + 20] =  0.000012559120;
	m_pBetaI[m_nBetaIOffset + 21] = -0.000007487521;
	m_pBetaI[m_nBetaIOffset + 22] =  0.000004463955;
	m_pBetaI[m_nBetaIOffset + 23] = -0.000002661442;
	m_pBetaI[m_nBetaIOffset + 24] =  0.000001586940;
	m_pBetaI[m_nBetaIOffset + 25] = -0.000000946531;
	m_pBetaI[m_nBetaIOffset + 26] =  0.000000565056;
	m_pBetaI[m_nBetaIOffset + 27] = -0.000000338154;
	m_pBetaI[m_nBetaIOffset + 28] =  0.000000203666;
	m_pBetaI[m_nBetaIOffset + 29] = -0.000000124662;
	m_pBetaI[m_nBetaIOffset + 30] =  0.000000079360;
	m_pBetaI[m_nBetaIOffset + 31] = -0.000000054833;
	m_pBetaI[m_nBetaIOffset + 32] =  0.000000018101;
	m_pBetaI[m_nBetaIOffset + 33] = -0.000000007299;
	m_pBetaI[m_nBetaIOffset + 34] =  0.000000002659;
	m_pBetaI[m_nBetaIOffset + 35] = -0.000000000765;
	m_pBetaI[m_nBetaIOffset + 36] =  0.000000000128;
	for( i=0; i<37; i++ )
		m_pBetaI[m_nBetaIOffset - 1 - i] = m_pBetaI[m_nBetaIOffset + i];

	return true;
}

// 実数部の補間の係数を読み取る関数
double CRISplineMra00::BetaR( int n )
{
	if( n<m_nBetaRBegin || n>m_nBetaREnd )
		return 0.0;
	else
		return m_pBetaR[n + m_nBetaROffset];
}

// 虚数部の補間の係数を読み取る関数
double CRISplineMra00::BetaI( int n )
{
	if( n<m_nBetaIBegin || n>m_nBetaIEnd )
		return 0.0;
	else
		return m_pBetaI[n + m_nBetaIOffset];
}

// 複素多重解像度解析
bool CRISplineMra00::Mra( double** output, double* input, int level, int mul, int lowMode )
{
	int i, j, k, L, n;
	int begin, end, sz;
	double sum;

	m_nN = (0x01 << level)*mul;
	m_nLevel = level;

	if( level < 1 || level > MAXLEVEL || m_nN > m_nMaxN || m_pBuf == NULL )
		return false;

	// 実数部
	begin = m_nBetaRBegin;		// 補間して入力
	end = m_nBetaREnd;
	for( j=0; j<m_nN; j++ )
	{
		sum = 0.0;
		for( k=begin; k<=end; k++ )
		{
			sum += BetaR(k) * Read( input, j-k, m_nN, 0 );
		}
		m_pWork[0][j] = sum;
	}

	for( i=1; i<=level; i++ )	// 分解アルゴリズムのループ
	{
		n = 0;
		sz = m_nN >> i;
		begin = m_pRealSpline->H0Begin();
		end = m_pRealSpline->H0End();
		for( k=0; k<sz; k++ )
		{
			sum = 0;
			for( L=2*k - end; L<=2*k - begin; L++ )
			{
				sum += m_pRealSpline->H0( 2*k - L )*Read(m_pWork[0], L, sz*2, 0);
			}
			m_pWork[1][n] = sum;
			n++;
		}

		begin = m_pRealSpline->H1Begin();
		end = m_pRealSpline->H1End();
		for( k=0; k<sz; k++ )
		{
			sum = 0;
			for( L=2*k - end; L<=2*k - begin; L++ )
			{
				double temp = m_pRealSpline->H1( 2*k - L );
				sum += m_pRealSpline->H1( 2*k - L )*Read(m_pWork[0], L, sz*2, 0);
			}
			m_pWork[1][n] = sum;
			n++;
		}

		for( j=0; j<sz; j++ )
		{
			m_pWork[0][j] = m_pWork[1][j];
		}
	}

	if( lowMode == 0 )
	{
		sz = m_nN >> level;		// 低域成分はそのまま出力
		for( j=0; j<sz; j++ )
		{
			output[0][j] = m_pWork[1][j];
		}
	}
	else
	{
		sz = m_nN >> level;		// 低域成分は逆補間して離散信号で出力
		begin = m_pRealSpline->Sc2Begin();
		end = m_pRealSpline->Sc2End();
		for( j=0; j<sz; j++ )
		{
			sum = 0.0;
			for( k=begin; k<=end; k++ )
			{
				sum += m_pRealSpline->Sc2(k) * Read( m_pWork[1], j-k, sz, 0 );
			}
			output[0][j] = sum;
		}
	}

	for( i=1; i<=level; i++ )
	{
		sz = m_nN >> i;
		for( j=0; j<sz; j++ )
		{
			output[0][j + sz] = m_pWork[1][j + sz];
		}
	}

	// 虚数部
	begin = m_nBetaIBegin;		// 補間して入力
	end = m_nBetaIEnd;
	for( j=0; j<m_nN; j++ )
	{
		sum = 0.0;
		for( k=begin; k<=end; k++ )
		{
			sum += BetaI(k) * Read( input, j-k, m_nN, 0 );
		}
		m_pWork[0][j] = sum;
	}

	for( i=1; i<=level; i++ )	// 分解アルゴリズムのループ
	{
		n = 0;
		sz = m_nN >> i;
		begin = m_pImagSpline->H0Begin();
		end = m_pImagSpline->H0End();
		for( k=0; k<sz; k++ )
		{
			sum = 0;
			for( L=2*k - end; L<=2*k - begin; L++ )
			{
				sum += m_pImagSpline->H0( 2*k - L )*Read(m_pWork[0], L, sz*2, 0);
			}
			m_pWork[1][n] = sum;
			n++;
		}

		begin = m_pRealSpline->H1Begin();
		end = m_pRealSpline->H1End();
		for( k=0; k<sz; k++ )
		{
			sum = 0 ;
			for( L=2*k - end; L<=2*k - begin; L++ )
			{
				sum += m_pImagSpline->H1( 2*k - L )*Read(m_pWork[0], L, sz*2, 0);
			}
			m_pWork[1][n] = sum;
			n++;
		}

		for( j=0; j<sz; j++ )
		{
			m_pWork[0][j] = m_pWork[1][j];
		}
	}

	if( lowMode == 0 )
	{
		sz = m_nN >> level;		// 低域成分はそのまま出力
		for( j=0; j<sz; j++ )
		{
			output[1][j] = m_pWork[1][j];
		}
	}
	else
	{
		sz = m_nN >> level;		// 低域成分は逆補間して離散信号で出力
		begin = m_pImagSpline->Sc2Begin();
		end = m_pImagSpline->Sc2End();
		for( j=0; j<sz; j++ )
		{
			sum = 0.0;
			for( k=begin; k<=end; k++ )
			{
				sum += m_pImagSpline->Sc2(k) * Read( m_pWork[1], j-k, sz, 0 );
			}
			output[1][j] = sum;
		}
	}

	for( i=1; i<=level; i++ )
	{
		sz = m_nN >> i;
		for( j=0; j<sz; j++ )
		{
			output[1][j + sz] = m_pWork[1][j + sz];
		}
	}

	return true;
}

// 複素数多重解像度解析の逆変換
bool CRISplineMra00::IMra( double* output, double** input, int level, int mul, int lowMode )
{
	int i,j,k,L, n;
	int begin1, end1, begin2, end2, sz;
	double sum;

	m_nN = (0x01 << level)*mul;
	m_nLevel = level;

	if( level < 1 || level > MAXLEVEL || m_nN > m_nMaxN || m_pBuf == NULL )
		return false;

	// 実数部
	if( lowMode == 0 )
	{
		sz = m_nN >> level;		// 低域成分はそのまま入力
		for( j=0; j<sz; j++ )
		{
			m_pWork[0][j] = input[0][j];
		}
	}
	else
	{
		sz = m_nN >> level;		// 低域成分は補間して入力
		begin1 = m_pRealSpline->BetaSc2Begin();
		end1 = m_pRealSpline->BetaSc2End();
		for( j=0; j<sz; j++ )
		{
			sum = 0.0;
			for( k=begin1; k<=end1; k++ )
			{
				sum += m_pRealSpline->BetaSc2(k) * Read( input[0], j-k, sz, 0 );
			}
			m_pWork[0][j] = sum;
		}
	}

	for( i=1; i<=level; i++ )
	{
		sz = m_nN >> i;
		for( j=0; j<sz; j++ )
		{
			m_pWork[0][j + sz] = input[0][j + sz];
		}
	}

	begin1 = m_pRealSpline->G0Begin();
	end1 = m_pRealSpline->G0End();
	begin2 = m_pRealSpline->G1Begin();
	end2 = m_pRealSpline->G1End();

	for( i=level ; i>0 ; i-- )	// 再構成アルゴリズムのループ
	{
		n = 0 ;
		sz = m_nN >> i;
		if( sz > 0 )
		{
			for( k=0 ; k<sz*2 ; k++ )
			{
				sum = 0 ;
				for( L=k/2 - end1/2 -5; L<=k/2 - begin1/2 +5; L++ )
					sum += m_pRealSpline->G0( k-2*L )*Read(m_pWork[0], L, sz, 0);

				for( L=k/2 - end2/2 -5; L<=k/2 - begin2/2 +5; L++ )
					sum += m_pRealSpline->G1( k-2*L )*Read(m_pWork[0], L, sz, sz);

				m_pWork[1][n] = sum;
				n++;
			}

			for( j=0 ; j<sz*2 ; j++ )
			{
				m_pWork[0][j] = m_pWork[1][j];
			}
		}
	}

	begin1 = m_pRealSpline->Sc1Begin();	// 逆補間して離散信号で出力
	end1 = m_pRealSpline->Sc1End();
	for( j=0; j<m_nN; j++ )
	{
		sum = 0.0;
		for( k=begin1; k<=end1; k++ )
		{
			sum += m_pRealSpline->Sc1(k) * Read( m_pWork[0], j-k, m_nN, 0 );
		}
		output[j] = sum;
	}

	// 虚数部
	if( lowMode == 0 )
	{
		sz = m_nN >> level;		// 低域成分はそのまま入力
		for( j=0; j<sz; j++ )
		{
			m_pWork[0][j] = input[1][j];
		}
	}
	else
	{
		sz = m_nN >> level;		// 低域成分は補間して入力
		begin1 = m_pImagSpline->BetaSc2Begin();
		end1 = m_pImagSpline->BetaSc2End();
		for( j=0; j<sz; j++ )
		{
			sum = 0.0;
			for( k=begin1; k<=end1; k++ )
			{
				sum += m_pImagSpline->BetaSc2(k) * Read( input[1], j-k, sz, 0 );
			}
			m_pWork[0][j] = sum;
		}
	}

	for( i=1; i<=level; i++ )
	{
		sz = m_nN >> i;
		for( j=0; j<sz; j++ )
		{
			m_pWork[0][j + sz] = input[1][j + sz];
		}
	}

	begin1 = m_pImagSpline->G0Begin();
	end1 = m_pImagSpline->G0End();
	begin2 = m_pImagSpline->G1Begin();
	end2 = m_pImagSpline->G1End();

	for( i=level ; i>0 ; i-- )	// 再構成アルゴリズムのループ
	{
		n = 0 ;
		sz = m_nN >> i;
		if( sz > 0 )
		{
			for( k=0 ; k<sz*2 ; k++ )
			{
				sum = 0 ;
				for( L=k/2 - end1/2 -5; L<=k/2 - begin1/2 +5; L++ )
					sum += m_pImagSpline->G0( k-2*L )*Read(m_pWork[0], L, sz, 0);

				for( L=k/2 - end2/2 -5; L<=k/2 - begin2/2 +5; L++ )
					sum += m_pImagSpline->G1( k-2*L )*Read(m_pWork[0], L, sz, sz);

				m_pWork[1][n] = sum;
				n++;
			}

			for( j=0 ; j<sz*2 ; j++ )
			{
				m_pWork[0][j] = m_pWork[1][j];
			}
		}
	}

	begin1 = m_pImagSpline->Sc1Begin();	// 逆補間して離散信号で出力
	end1 = m_pImagSpline->Sc1End();
	for( j=0; j<m_nN; j++ )
	{
		sum = 0.0;
		for( k=begin1; k<=end1; k++ )
		{
			sum += m_pImagSpline->Sc1(k) * Read( m_pWork[0], j-k, m_nN, 0 );
		}
		output[j] += sum;
	}

	return true;
}

// 実数部多重解像度解析
bool CRISplineMra00::RealMra( double* output, double* input, int level, int mul, int lowMode )
{
	int i, j, k, L, n;
	int begin, end, sz;
	double sum;

	m_nN = (0x01 << level)*mul;
	m_nLevel = level;

	if( level < 1 || level > MAXLEVEL || m_nN > m_nMaxN || m_pBuf == NULL )
		return false;

	if( lowMode == 0 )
	{
		for( j=0; j<m_nN; j++ )				// そのまま入力
		{
			m_pWork[0][j] = input[j];
		}
	}
	else
	{
		begin = m_pRealSpline->BetaSc2Begin();	// 補間して入力
		end = m_pRealSpline->BetaSc2End();
		for( j=0; j<m_nN; j++ )
		{
			sum = 0.0;
			for( k=begin; k<=end; k++ )
			{
				sum += m_pRealSpline->BetaSc2(k) * Read( input, j-k, m_nN, 0 );
			}
			m_pWork[0][j] = sum;
		}
	}

	for( i=1; i<=level; i++ )	// 分解アルゴリズムのループ
	{
		n = 0;
		sz = m_nN >> i;
		begin = m_pRealSpline->H0Begin();
		end = m_pRealSpline->H0End();
		for( k=0; k<sz; k++ )
		{
			sum = 0;
			for( L=2*k - end; L<=2*k - begin; L++ )
			{
				sum += m_pRealSpline->H0( 2*k - L )*Read(m_pWork[0], L, sz*2, 0);
			}
			m_pWork[1][n] = sum;
			n++;
		}

		begin = m_pRealSpline->H1Begin();
		end = m_pRealSpline->H1End();
		for( k=0; k<sz; k++ )
		{
			sum = 0 ;
			for( L=2*k - end; L<=2*k - begin; L++ )
			{
				sum += m_pRealSpline->H1( 2*k - L )*Read(m_pWork[0], L, sz*2, 0);
			}
			m_pWork[1][n] = sum;
			n++;
		}

		for( j=0; j<sz; j++ )
		{
			m_pWork[0][j] = m_pWork[1][j];
		}
	}

	if( lowMode == 0 )
	{
		sz = m_nN >> level;		// 低域成分はそのまま出力
		for( j=0; j<sz; j++ )
		{
			output[j] = m_pWork[1][j];
		}
	}
	else
	{
		sz = m_nN >> level;		// 低域成分は逆補間して離散信号で出力
		begin = m_pRealSpline->Sc2Begin();
		end = m_pRealSpline->Sc2End();
		for( j=0; j<sz; j++ )
		{
			sum = 0.0;
			for( k=begin; k<=end; k++ )
			{
				sum += m_pRealSpline->Sc2(k) * Read( m_pWork[1], j-k, sz, 0 );
			}
			output[j] = sum;
		}
	}

	for( i=1; i<=level; i++ )
	{
		sz = m_nN >> i;
		for( j=0; j<sz; j++ )
		{
			output[j + sz] = m_pWork[1][j + sz];
		}
	}

	return true;
}

// 実数部多重解像度解析の逆変換
bool CRISplineMra00::RealIMra( double* output, double* input, int level, int mul, int lowMode )
{
	int i,j,k,L, n;
	int begin1, end1, begin2, end2, sz;
	double sum;

	m_nN = (0x01 << level)*mul;
	m_nLevel = level;

	if( level < 1 || level > MAXLEVEL || m_nN > m_nMaxN || m_pBuf == NULL )
		return false;

	if( lowMode == 0 )
	{
		sz = m_nN >> level;		// 低域成分はそのまま入力
		for( j=0; j<sz; j++ )
		{
			m_pWork[0][j] = input[j];
		}
	}
	else
	{
		sz = m_nN >> level;		// 低域成分は補間して入力
		begin1 = m_pRealSpline->BetaSc2Begin();
		end1 = m_pRealSpline->BetaSc2End();
		for( j=0; j<sz; j++ )
		{
			sum = 0.0;
			for( k=begin1; k<=end1; k++ )
			{
				sum += m_pRealSpline->BetaSc2(k) * Read( input, j-k, sz, 0 );
			}
			m_pWork[0][j] = sum;
		}
	}

	for( i=1; i<=level; i++ )
	{
		sz = m_nN >> i;
		for( j=0; j<sz; j++ )
		{
			m_pWork[0][j + sz] = input[j + sz];
		}
	}

	begin1 = m_pRealSpline->G0Begin();
	end1 = m_pRealSpline->G0End();
	begin2 = m_pRealSpline->G1Begin();
	end2 = m_pRealSpline->G1End();

	for( i=level ; i>0 ; i-- )	// 再構成アルゴリズムのループ
	{
		n = 0 ;
		sz = m_nN >> i;
		if( sz > 0 )
		{
			for( k=0 ; k<sz*2 ; k++ )
			{
				sum = 0 ;
				for( L=k/2 - end1/2 -5; L<=k/2 - begin1/2 +5; L++ )
					sum += m_pRealSpline->G0( k-2*L )*Read(m_pWork[0], L, sz, 0);

				for( L=k/2 - end2/2 -5; L<=k/2 - begin2/2 +5; L++ )
					sum += m_pRealSpline->G1( k-2*L )*Read(m_pWork[0], L, sz, sz);

				m_pWork[1][n] = sum;
				n++;
			}

			for( j=0 ; j<sz*2 ; j++ )
			{
				m_pWork[0][j] = m_pWork[1][j];
			}
		}
	}

	if( lowMode == 0 )
	{
		for( j=0; j<m_nN; j++ )				// そのまま出力
		{
			output[j] = m_pWork[0][j];
		}
	}
	else
	{
		begin1 = m_pRealSpline->Sc2Begin();	// 逆補間して離散信号で出力
		end1 = m_pRealSpline->Sc2End();
		for( j=0; j<m_nN; j++ )
		{
			sum = 0.0;
			for( k=begin1; k<=end1; k++ )
			{
				sum += m_pRealSpline->Sc2(k) * Read( m_pWork[0], j-k, m_nN, 0 );
			}
			output[j] = sum;
		}
	}

	return true;
}

// 虚数部多重解像度解析
bool CRISplineMra00::ImagMra( double* output, double* input, int level, int mul, int lowMode )
{
	int i, j, k, L, n;
	int begin, end, sz;
	double sum;

	m_nN = (0x01 << level)*mul;
	m_nLevel = level;

	if( level < 1 || level > MAXLEVEL || m_nN > m_nMaxN || m_pBuf == NULL )
		return false;

	if( lowMode == 0 )
	{
		for( j=0; j<m_nN; j++ )				// そのまま入力
		{
			m_pWork[0][j] = input[j];
		}
	}
	else
	{
		begin = m_pImagSpline->BetaSc2Begin();	// 補間して入力
		end = m_pImagSpline->BetaSc2End();
		for( j=0; j<m_nN; j++ )
		{
			sum = 0.0;
			for( k=begin; k<=end; k++ )
			{
				sum += m_pImagSpline->BetaSc2(k) * Read( input, j-k, m_nN, 0 );
			}
			m_pWork[0][j] = sum;
		}
	}

	for( i=1; i<=level; i++ )	// 分解アルゴリズムのループ
	{
		n = 0;
		sz = m_nN >> i;
		begin = m_pImagSpline->H0Begin();
		end = m_pImagSpline->H0End();
		for( k=0; k<sz; k++ )
		{
			sum = 0;
			for( L=2*k - end; L<=2*k - begin; L++ )
			{
				sum += m_pImagSpline->H0( 2*k - L )*Read(m_pWork[0], L, sz*2, 0);
			}
			m_pWork[1][n] = sum;
			n++;
		}

		begin = m_pImagSpline->H1Begin();
		end = m_pImagSpline->H1End();
		for( k=0; k<sz; k++ )
		{
			sum = 0 ;
			for( L=2*k - end; L<=2*k - begin; L++ )
			{
				sum += m_pImagSpline->H1( 2*k - L )*Read(m_pWork[0], L, sz*2, 0);
			}
			m_pWork[1][n] = sum;
			n++;
		}

		for( j=0; j<sz; j++ )
		{
			m_pWork[0][j] = m_pWork[1][j];
		}
	}

	if( lowMode == 0 )
	{
		sz = m_nN >> level;		// 低域成分はそのまま出力
		for( j=0; j<sz; j++ )
		{
			output[j] = m_pWork[1][j];
		}
	}
	else
	{
		sz = m_nN >> level;		// 低域成分は逆補間して離散信号で出力
		begin = m_pImagSpline->Sc2Begin();
		end = m_pImagSpline->Sc2End();
		for( j=0; j<sz; j++ )
		{
			sum = 0.0;
			for( k=begin; k<=end; k++ )
			{
				sum += m_pImagSpline->Sc2(k) * Read( m_pWork[1], j-k, sz, 0 );
			}
			output[j] = sum;
		}
	}

	for( i=1; i<=level; i++ )
	{
		sz = m_nN >> i;
		for( j=0; j<sz; j++ )
		{
			output[j + sz] = m_pWork[1][j + sz];
		}
	}

	return true;
}

// 虚数部多重解像度解析の逆変換
bool CRISplineMra00::ImagIMra( double* output, double* input, int level, int mul, int lowMode )
{
	int i,j,k,L, n;
	int begin1, end1, begin2, end2, sz;
	double sum;

	m_nN = (0x01 << level)*mul;
	m_nLevel = level;

	if( level < 1 || level > MAXLEVEL || m_nN > m_nMaxN || m_pBuf == NULL )
		return false;

	if( lowMode == 0 )
	{
		sz = m_nN >> level;		// 低域成分はそのまま入力
		for( j=0; j<sz; j++ )
		{
			m_pWork[0][j] = input[j];
		}
	}
	else
	{
		sz = m_nN >> level;		// 低域成分は補間して入力
		begin1 = m_pImagSpline->BetaSc2Begin();
		end1 = m_pImagSpline->BetaSc2End();
		for( j=0; j<sz; j++ )
		{
			sum = 0.0;
			for( k=begin1; k<=end1; k++ )
			{
				sum += m_pImagSpline->BetaSc2(k) * Read( input, j-k, sz, 0 );
			}
			m_pWork[0][j] = sum;
		}
	}

	for( i=1; i<=level; i++ )
	{
		sz = m_nN >> i;
		for( j=0; j<sz; j++ )
		{
			m_pWork[0][j + sz] = input[j + sz];
		}
	}

	begin1 = m_pImagSpline->G0Begin();
	end1 = m_pImagSpline->G0End();
	begin2 = m_pImagSpline->G1Begin();
	end2 = m_pImagSpline->G1End();

	for( i=level ; i>0 ; i-- )	// 再構成アルゴリズムのループ
	{
		n = 0 ;
		sz = m_nN >> i;
		if( sz > 0 )
		{
			for( k=0 ; k<sz*2 ; k++ )
			{
				sum = 0 ;
				for( L=k/2 - end1/2 -5; L<=k/2 - begin1/2 +5; L++ )
					sum += m_pImagSpline->G0( k-2*L )*Read(m_pWork[0], L, sz, 0);

				for( L=k/2 - end2/2 -5; L<=k/2 - begin2/2 +5; L++ )
					sum += m_pImagSpline->G1( k-2*L )*Read(m_pWork[0], L, sz, sz);

				m_pWork[1][n] = sum;
				n++;
			}

			for( j=0 ; j<sz*2 ; j++ )
			{
				m_pWork[0][j] = m_pWork[1][j];
			}
		}
	}

	if( lowMode == 0 )
	{
		for( j=0; j<m_nN; j++ )				// そのまま出力
		{
			output[j] = m_pWork[0][j];
		}
	}
	else
	{
		begin1 = m_pImagSpline->Sc2Begin();	// 逆変換して離散信号で出力
		end1 = m_pImagSpline->Sc2End();
		for( j=0; j<m_nN; j++ )
		{
			sum = 0.0;
			for( k=begin1; k<=end1; k++ )
			{
				sum += m_pImagSpline->Sc2(k) * Read( m_pWork[0], j-k, m_nN, 0 );
			}
			output[j] = sum;
		}
	}

	return true;
}

// ループ状にデータを読み込む
double CRISplineMra00::Read( double* buf, int x, int sz, int ofst )
{
	while( x < 0)
		x += sz;
	while( x >= sz)
		x -= sz;

	return buf[x + ofst];
}

// ループ状にデータを書き込む
void CRISplineMra00::Write( double* buf, double d, int x, int sz, int ofst )
{
	while( x < 0 )
		x += sz;
	while( x >= sz )
		x -= sz;

	buf[x + ofst ] = d;
}
