#include <math.h>
#include "GraphRISplineMraE00.h"
#include "RISplineMra00.h"

#define MAXLEVEL 8

CGraphRISplineMraE00::CGraphRISplineMraE00()
{
	int i;

	m_pMra = new CRISplineMra00;
	
	for( i=0; i<3; i++ )
		m_pBuf[i] = NULL;

	m_pWork = NULL;
}

CGraphRISplineMraE00::~CGraphRISplineMraE00()
{
	int i;

	if( m_pMra != NULL )
		delete m_pMra;

	for( i=0; i<3; i++ )
	{
		if( m_pBuf[i] != NULL )
			delete m_pBuf[i];
	}

	if( m_pWork != NULL )
		delete m_pWork;
}

bool CGraphRISplineMraE00::Prepare( int maxX, int maxY )
{
	int i;

	m_nMaxX = maxX;
	m_nMaxY = maxY;

	m_nMaxLen = m_nMaxX;
	if( m_nMaxY > m_nMaxX )
		m_nMaxLen = m_nMaxY;

	m_nMaxN = m_nMaxX * m_nMaxY;

	m_pMra->Prepare( m_nMaxLen );

	for( i=0; i<3; i++ )
	{
		m_pBuf[i] = new double[ m_nMaxN*4 + 256 ];
		if( m_pBuf[i] == NULL )
		{
			return false;
		}
	}

	m_pBox0[0][0] = m_pBuf[0];
	m_pBox0[0][1] = m_pBox0[0][0] + m_nMaxN;
	m_pBox0[1][0] = m_pBox0[0][1] + m_nMaxN;
	m_pBox0[1][1] = m_pBox0[1][0] + m_nMaxN;

	m_pBox1[0][0] = m_pBuf[1];
	m_pBox1[0][1] = m_pBox1[0][0] + m_nMaxN;
	m_pBox1[1][0] = m_pBox1[0][1] + m_nMaxN;
	m_pBox1[1][1] = m_pBox1[1][0] + m_nMaxN;

	m_pBox2[0][0] = m_pBuf[2];
	m_pBox2[0][1] = m_pBox2[0][0] + m_nMaxN;
	m_pBox2[1][0] = m_pBox2[0][1] + m_nMaxN;
	m_pBox2[1][1] = m_pBox2[1][0] + m_nMaxN;

	m_pWork = new double[ m_nMaxLen * 16 ];
	if( m_pWork == NULL )
		return false;

	return true;
}

// 2次元複素数多重解像度解析
bool CGraphRISplineMraE00::Mra(double* output, double* input, int nX, int nY, int level)
{
	int i, j, L, x, y;
	int xsz, ysz;

	if( nX > m_nMaxX || nY > m_nMaxY || level < 1 || level > MAXLEVEL )
		return false;

	double* work0;
	double* work1;
	double* work2[2];

	work0 = m_pWork;
	work1 = work0 + m_nMaxLen;
	work2[0] = work1 + m_nMaxLen;
	work2[1] = work2[0] + m_nMaxLen;

	m_nX = nX;
	m_nY = nY;
	m_nLevel = level;

	// 入力
	for( i=0; i<m_nX*m_nY; i++ )
	{
		m_pBox0[0][0][i] = input[i];
	}

	// level -1 の処理
	xsz = m_nX;
	ysz = m_nY;
	for( y=0; y<ysz; y++ )
	{
		for( x=0; x<xsz; x++ )
		{
			work0[x] = m_pBox0[0][0][ x + y*m_nX ];
		}
		m_pMra->Mra( work2, work0, 1, xsz >> 1, 1 );
		for( x=0; x<xsz; x++ )
		{
			m_pBox1[0][0][ x + y*m_nX ] = work2[0][x];
			m_pBox1[0][1][ x + y*m_nX ] = work2[1][x];
		}
	}

	for( x=0; x<xsz; x++ )
	{
		for( y=0; y<ysz; y++ )
		{
			work0[y] = m_pBox1[0][0][ x + y*m_nX ];
		}
		m_pMra->Mra( work2, work0, 1, ysz >> 1, 1 );
		for( y=0; y<ysz; y++ )
		{
			m_pBox2[0][0][ x + y*m_nX ] = work2[0][y];
			m_pBox2[0][1][ x + y*m_nX ] = work2[1][y];
		}

		for( y=0; y<ysz; y++ )
		{
			work0[y] = m_pBox1[0][1][ x + y*m_nX ];
		}
		m_pMra->Mra( work2, work0, 1, ysz >> 1, 1 );
		for( y=0; y<ysz; y++ )
		{
			m_pBox2[1][0][ x + y*m_nX ] = work2[0][y];
			m_pBox2[1][1][ x + y*m_nX ] = work2[1][y];
		}
	}

	for( i=0; i<2; i++ )
	{
		for( j=0; j<2; j++ )
		{
			for( y=0; y<ysz; y++ )
			{
				for( x=0; x<xsz; x++ )
				{
					m_pBox0[i][j][ x + y*m_nX ] = m_pBox2[i][j][x + y*m_nX];
				}
			}
		}
	}

	// level -2 以降の処理
	for( L=2; L<=m_nLevel; L++ )
	{
		xsz = m_nX >> (L - 1);
		ysz = m_nY >> (L - 1);
		for( i=0; i<2; i++ )
		{
			for( j=0; j<2; j++ )
			{
				for( y=0; y<ysz; y++ )
				{
					for( x=0; x<xsz; x++ )
					{
						work0[x] = m_pBox0[i][j][ x + y*m_nX ];
					}
					SelectMra( work1, work0, xsz, i );
					for( x=0; x<xsz; x++ )
					{
						m_pBox1[i][j][ x + y*m_nX ] = work1[x];
					}
				}
				for( x=0; x<xsz; x++ )
				{
					for( y=0; y<ysz; y++ )
					{
						work0[y] = m_pBox1[i][j][ x + y*m_nX ];
					}
					SelectMra( work1, work0, ysz, j );
					for( y=0; y<ysz; y++ )
					{
						m_pBox2[i][j][ x + y*m_nX ] = work1[y];
					}
				}
			}
		}

		for( i=0; i<2; i++ )
		{
			for( j=0; j<2; j++ )
			{
				for( y=0; y<ysz; y++ )
				{
					for( x=0; x<xsz; x++ )
					{
						m_pBox0[i][j][ x + y*m_nX ] = m_pBox2[i][j][x + y*m_nX];
					}
				}
			}
		}
	}

	// 出力
	for( i=0; i<2; i++ )
	{
		for( j=0; j<2; j++ )
		{
			for( y=0; y<m_nY; y++ )
			{
				for( x=0; x<m_nX; x++ )
				{
					output[ (m_nX*j + x) + (m_nY*i + y)*m_nX*2 ]
						= m_pBox2[i][j][x + y*m_nX];
				}
			}
		}
	}

	return true;
}

// 2次元複素数多重解像度解析の逆変換
bool CGraphRISplineMraE00::IMra(double* output, double* input, int nX, int nY, int level)
{
	int i, j, L, x, y;
	int xsz, ysz;

	if( nX > m_nMaxX || nY > m_nMaxY || level < 1 || level > MAXLEVEL )
		return false;

	double* work0;
	double* work1;
	double* work2[2];

	work0 = m_pWork;
	work1 = work0 + m_nMaxLen;
	work2[0] = work1 + m_nMaxLen;
	work2[1] = work2[0] + m_nMaxLen;

	m_nX = nX;
	m_nY = nY;
	m_nLevel = level;

	// 入力
	for( i=0; i<2; i++ )
	{
		for( j=0; j<2; j++ )
		{
			for( y=0; y<m_nY; y++ )
			{
				for( x=0; x<m_nX; x++ )
				{
					m_pBox0[i][j][ x + y*m_nX ]
						= input[ (m_nX*j + x) + (m_nY*i + y)*m_nX*2 ];
				}
			}
		}
	}

	// level -2 までの処理
	for( L=m_nLevel; L>=2; L-- )
	{
		xsz = m_nX >> (L - 1);
		ysz = m_nY >> (L - 1);
		for( i=0; i<2; i++ )
		{
			for( j=0; j<2; j++ )
			{
				for( x=0; x<xsz; x++ )
				{
					for( y=0; y<ysz; y++ )
					{
						work0[y] = m_pBox0[i][j][ x + y*m_nX ];
					}
					SelectIMra( work1, work0, ysz, j );
					for( y=0; y<ysz; y++ )
					{
						m_pBox1[i][j][ x + y*m_nX ] = work1[y];
					}
				}
				for( y=0; y<ysz; y++ )
				{
					for( x=0; x<xsz; x++ )
					{
						work0[x] = m_pBox1[i][j][ x + y*m_nX ];
					}
					SelectIMra( work1, work0, xsz, i );
					for( x=0; x<xsz; x++ )
					{
						m_pBox2[i][j][ x + y*m_nX ] = work1[x];
					}
				}
			}
		}

		for( i=0; i<2; i++ )
		{
			for( j=0; j<2; j++ )
			{
				for( y=0; y<ysz; y++ )
				{
					for( x=0; x<xsz; x++ )
					{
						m_pBox0[i][j][ x + y*m_nX ] = m_pBox2[i][j][x + y*m_nX];
					}
				}
			}
		}
	}

	// level -1 の処理
	xsz = m_nX;
	ysz = m_nY;
	for( x=0; x<xsz; x++ )
	{
		for( y=0; y<ysz; y++ )
		{
			work2[0][y] = m_pBox0[0][0][ x + y*m_nX ];
			work2[1][y] = m_pBox0[0][1][ x + y*m_nX ];
		}
		m_pMra->IMra( work1, work2, 1, ysz >> 1, 1 );
		for( y=0; y<ysz; y++ )
		{
			m_pBox1[0][0][ x + y*m_nX ] = work1[y];
		}

		for( y=0; y<ysz; y++ )
		{
			work2[0][y] = m_pBox0[1][0][ x + y*m_nX ];
			work2[1][y] = m_pBox0[1][1][ x + y*m_nX ];
		}
		m_pMra->IMra( work1, work2, 1, ysz >> 1, 1 );
		for( y=0; y<ysz; y++ )
		{
			m_pBox1[0][1][ x + y*m_nX ] = work1[y];
		}
	}

	for( y=0; y<ysz; y++ )
	{
		for( x=0; x<xsz; x++ )
		{
			work2[0][x] = m_pBox1[0][0][ x + y*m_nX ];
			work2[1][x] = m_pBox1[0][1][ x + y*m_nX ];
		}
		m_pMra->IMra( work1, work2, 1, xsz >> 1, 1 );
		for( x=0; x<xsz; x++ )
		{
			m_pBox2[0][0][ x + y*m_nX ] = work1[x];
		}
	}

	// 出力
	for( i=0; i<nX*nY; i++ )
	{
		output[i] = m_pBox2[0][0][i];
	}

	return true;
}

// 実数部、虚数部を選択して多重解像度解析を実行
void CGraphRISplineMraE00::SelectMra(double* output, double* input, int n, int ri)
{
	if( ri == 0 )
		m_pMra->RealMra( output, input, 1, n >> 1, 1 );
	else
		m_pMra->ImagMra( output, input, 1, n >> 1, 1 );
}

// 実数部、虚数部を選択して多重解像度解析の逆変換を実行
void CGraphRISplineMraE00::SelectIMra(double* output, double* input, int n, int ri)
{
	if( ri == 0 )
		m_pMra->RealIMra( output, input, 1, n >> 1, 1 );
	else
		m_pMra->ImagIMra( output, input, 1, n >> 1, 1 );
}
