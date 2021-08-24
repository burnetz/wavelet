#include <Math.h>
#include "MatchingPursuit00.h"
#include "OrthoWavelet00.h"
#include "OrthoMra00.h"
#include "Interpolation00.h"

#define MAXLEVEL 8

CMatchingPursuit00::CMatchingPursuit00()
{
	m_pOrtho = new COrthoWavelet00;
	m_pMra = new COrthoMra00;
	m_pInter = new CInterpolation00;

	m_nEpsilon = 0.0;

	m_pBuf = NULL;
	m_nN = 0;
	m_nBOffset = 0;
	m_nUOffset = 0;
}

CMatchingPursuit00::~CMatchingPursuit00()
{
	if (m_pOrtho != NULL)
		delete m_pOrtho;
	if (m_pMra != NULL)
		delete m_pMra;
	if (m_pInter != NULL)
		delete m_pInter;
	if (m_pBuf != NULL)
		delete m_pBuf;
}

bool CMatchingPursuit00::Prepare(int wName, int maxN, int maxLevel)
{
	double *imraIn;
	double *imraOut;
	int i, j, x, x1;
	double temp;

	m_nMaxN = maxN;
	m_nMaxLevel = maxLevel;

	if (maxLevel < 1 || maxLevel > MAXLEVEL)
		return false;
	if (m_pOrtho->Prepare(wName, maxLevel) == false)
		return false;
	if (m_pInter->Prepare() == false)
		return false;

	int scalbegin = m_pOrtho->ScalingBegin(0);
	int scalend = m_pOrtho->ScalingEnd(0);
	m_nScalOffset = -scalbegin;
	m_nScalN = scalend - scalbegin;

	int mSz = m_pOrtho->MotherEnd(maxLevel) - m_pOrtho->MotherBegin(maxLevel);

	int bEx = 0;
	for (i = 0; i < 16; i++)
	{
		if ((0x01 << i) > mSz)
			break;
	}
	if (i >= 16)
		return false;
	bEx = i + 1;
	int bSz = 0x01 << bEx;

	// 多重解像度解析の準備
	if (m_pMra->Prepare(wName, bSz) == false)
		return false;
	int uSz = mSz * 2;

	// バッファ確保
	if (m_pBuf != NULL)
		delete m_pBuf;
	if ((m_pBuf = new double[bSz * 2 + (maxLevel + 1) * bSz + (maxLevel + 1) * (maxLevel + 1) * uSz + (maxLevel + 1) * maxN + m_nScalN + maxN + 256]) != NULL)
	{
		imraIn = m_pBuf;
		imraOut = imraIn + bSz;
		m_pBBuf[0] = imraOut + bSz;
		for (i = 1; i < maxLevel + 1; i++)
			m_pBBuf[i] = m_pBBuf[0] + i * bSz;
		m_pUBuf[0][0] = m_pBBuf[0] + (maxLevel + 1) * bSz;
		for (i = 0; i < maxLevel + 1; i++)
			for (j = 0; j < maxLevel + 1; j++)
				m_pUBuf[i][j] = m_pUBuf[0][0] + (i * (maxLevel + 1) + j) * uSz;
		m_pDictBuf[0] = m_pUBuf[0][0] + (maxLevel + 1) * (maxLevel + 1) * uSz;
		for (i = 1; i < maxLevel + 1; i++)
			m_pDictBuf[i] = m_pDictBuf[0] + i * maxN;
		m_pScal = m_pDictBuf[0] + (maxLevel + 1) * maxN;
		m_pLowPassBuf = m_pScal + m_nScalN;
	}
	else
	{
		return false;
	}

	// 再構成アルゴリズムを実行して係数を収集
	for (i = 0; i <= maxLevel; i++)
	{
		for (j = 0; j < bSz; j++)
		{
			imraIn[j] = 0.0;
		}

		if (i)
		{
			// ウェーブレット係数を設定
			imraIn[(bSz >> i) + (bSz >> i) / 2] = 1.0;
		}
		else
		{
			// m_pBBuf[0]には低周波フィルタ係数を入れるために
			// iが0のときはスケーリング関数係数を設定
			imraIn[(bSz >> maxLevel) / 2] = 1.0;
		}

		m_pMra->IMra(imraOut, imraIn, maxLevel, bSz >> maxLevel, 0);
		for (j = 0; j < bSz; j++)
		{
			m_pBBuf[i][j] = imraOut[j];
		}
	}

	// 係数bの開始と終了を調査
	m_nBOffset = bSz / 2;
	for (i = 0; i <= maxLevel; i++)
	{
		m_nBBegin[i] = 0;
		for (j = 0; j < bSz; j++)
		{
			if (m_pBBuf[i][j] != 0.0)
			{
				m_nBBegin[i] = j;
				break;
			}
		}
		m_nBEnd[i] = bSz - 1;
		for (j = bSz - 1; j >= 0; j--)
		{
			if (m_pBBuf[i][j] != 0.0)
			{
				m_nBEnd[i] = j;
				break;
			}
		}
		if (m_nBEnd[i] < bSz - 1)
			m_nBEnd[i] += 1;
		m_nBBegin[i] -= m_nBOffset;
		m_nBEnd[i] -= m_nBOffset;
	}

	// 係数uの計算
	m_nUOffset = uSz / 2;
	for (i = 0; i <= maxLevel; i++)
	{
		int begin = m_nBBegin[i];
		int end = m_nBEnd[i];
		for (j = 0; j <= maxLevel; j++)
		{
			int begin1 = m_nBBegin[j];
			int end1 = m_nBEnd[j];
			for (x = -m_nUOffset; x < m_nUOffset; x++)
			{
				int begin2 = begin1 + x;
				int end2 = end1 + x;
				int maxBegin = begin;
				if (begin2 > begin)
					maxBegin = begin2;
				int minEnd = end;
				if (end2 < end)
					minEnd = end2;
				if (maxBegin <= minEnd)
				{
					temp = 0.0;
					for (x1 = maxBegin; x1 < minEnd; x1++)
					{
						temp += ReadB(i, x1) * ReadB(j, x1 - x);
					}
					m_pUBuf[i][j][x + m_nUOffset] = temp;
				}
				else
				{
					m_pUBuf[i][j][x + m_nUOffset] = 0.0;
				}
			}
		}
	}

	// 係数uの開始と終了を調査
	for (i = 0; i <= maxLevel; i++)
	{
		for (j = 0; j <= maxLevel; j++)
		{
			m_nUBegin[i][j] = 0;
			for (x = 0; x < uSz; x++)
			{
				if (m_pUBuf[i][j][x] != 0.0)
				{
					m_nUBegin[i][j] = x;
					break;
				}
			}
			m_nUEnd[i][j] = uSz - 1;
			for (x = uSz - 1; x >= 0; x--)
			{
				if (m_pUBuf[i][j][x] != 0.0)
				{
					m_nUEnd[i][j] = x;
					break;
				}
			}
			if (m_nUEnd[i][j] < uSz - 1)
				m_nUEnd[i][j] += 1;
			m_nUBegin[i][j] -= m_nUOffset;
			m_nUEnd[i][j] -= m_nUOffset;
		}
	}

	for (i = 0; i < m_nScalN; i++)
		m_pScal[i] = m_pOrtho->Scaling(0, i - m_nScalOffset);

	return true;
}

// マッチング追跡 (double* output[]はlevel+2個確保すること)
bool CMatchingPursuit00::Pursuit(double *output[], double *input,
								int inputN, int level, double alpha, double epsilon, double rho, int mode)
{
	int i, j, k, L, x;
	double temp;

	m_nEpsilon = 0.0;
	m_nN = inputN;
	m_nLevel = level;

	if (m_pBuf == NULL || level > m_nMaxLevel || m_nN > m_nMaxN)
		return false;

	// rhoよりループのカウンター数を設定
	m_nCount = m_nN - (m_nN >> level);
	m_nCount = (int)(rho * m_nCount);

	if (mode == 0) // 補間なし
	{
		for (i = 0; i < m_nN; i++)
			m_pDictBuf[0][i] = input[i];
	}
	else // 補間あり
	{
		m_pInter->Inter(m_pDictBuf[0], input, m_nN,
						m_pScal, m_nScalN, m_nScalOffset, 1);
	}

	// 低域分離
	for (i = 0; i < m_nN; i++)
	{
		temp = 0.0;
		for (j = m_nBBegin[0]; j <= m_nBEnd[0]; j++)
		{
			x = i + j;
			while (x < 0)
				x += inputN;
			while (x >= inputN)
				x -= inputN;

			temp += m_pBBuf[0][j + m_nBOffset] * m_pDictBuf[0][x];
		}
		m_pDictBuf[1][i] = temp / (double)(0x01 << level);
	} // 未使用のm_pDictBuf[1]を一時的に借用

	for (i = 0; i < m_nN; i++) // 低域分離2回目の操作
	{
		temp = 0.0;
		for (j = m_nBBegin[0]; j <= m_nBEnd[0]; j++)
		{
			x = i - j;
			while (x < 0)
				x += inputN;
			while (x >= inputN)
				x -= inputN;

			temp += m_pBBuf[0][j + m_nBOffset] * m_pDictBuf[1][x];
		}
		m_pLowPassBuf[i] = temp; // cl0
	}

	for (i = 0; i < m_nN; i++)
	{
		m_pDictBuf[0][i] -= m_pLowPassBuf[i]; // ch0
	}

	// エネルギー計算
	double f2 = 0.0;
	for (i = 0; i < m_nN; i++)
	{
		temp = ReadCZero(i);
		f2 += temp * temp;
	}
	double rf2 = f2;

	for (L = 1; L <= level; L++)
	{
		for (x = 0; x < m_nN; x++)
		{
			output[L][x] = 0.0;
		}
	}

	// 辞書作成
	for (i = 1; i <= level; i++)
	{
		for (j = 0; j < m_nN; j++)
		{
			temp = 0.0;
			for (k = m_nBBegin[i]; k < m_nBEnd[i]; k++)
			{
				temp += ReadB(i, k) * ReadCZero(j + k);
			}
			m_pDictBuf[i][j] = temp;
		}
	}

	// マッチング追跡のメインループ
	int count = 0;
	while (count < m_nCount && rf2 > epsilon * epsilon * f2)
	{
		int maxL = 0;
		int maxX = 0;
		double maxData = 0.0;
		for (L = 1; L <= level; L++)
		{
			for (x = 0; x < m_nN; x++)
			{
				temp = m_pDictBuf[L][x];
				if (temp < 0)
					temp = -temp;
				if (maxData < temp)
				{
					maxL = L;
					maxX = x;
					maxData = temp;
				}
			}
		}
		if (maxData == 0.0)
			break;

		for (L = 1; L <= level; L++)
		{
			for (x = 0; x < m_nN; x++)
			{
				temp = m_pDictBuf[L][x];
				if (temp < 0)
					temp = -temp;
				if (temp > maxData * alpha)
				{
					double temp1 = m_pDictBuf[L][x];
					output[L][x] += temp1;
					rf2 -= temp1 * temp1;
					UpdateBuf(L, x);
					count++;
					if (count >= m_nCount || rf2 <= epsilon * epsilon * f2)
						break;
				}
			}
			if (count >= m_nCount || rf2 <= epsilon * epsilon * f2)
				break;
		}
	}

	// 残留成分を計算
	for (L = 1; L <= level; L++)
	{
		for (i = 0; i < m_nN; i++)
		{
			for (j = m_nBBegin[L]; j < m_nBEnd[L]; j++)
			{
				x = i + j;
				while (x < 0)
					x += m_nN;
				while (x >= m_nN)
					x -= m_nN;
				m_pDictBuf[0][x] -= output[L][i] * ReadB(L, j);
			}
		}
	}

	// 真のεを計算
	double rf22 = 0.0;
	for (x = 0; x < m_nN; x++)
	{
		rf22 += m_pDictBuf[0][x] * m_pDictBuf[0][x];
	}
	m_nEpsilon = sqrt(rf22 / f2);

	// 低周波成分と残留成分の出力
	for (x = 0; x < m_nN; x++)
	{
		output[0][x] = m_pLowPassBuf[x];
		output[level + 1][x] = m_pDictBuf[0][x];
	}

	return true;
}

// 係数の読み込み
double CMatchingPursuit00::ReadB(int level, int n)
{
	if (n < m_nBBegin[level] || n >= m_nBEnd[level])
		return 0.0;
	else
		return m_pBBuf[level][n + m_nBOffset];
}

// ch0の読み込み
double CMatchingPursuit00::ReadCZero(int x)
{
	while (x < 0)
		x += m_nN;

	while (x >= m_nN)
		x -= m_nN;

	return m_pDictBuf[0][x];
}

// 辞書のアップデート
void CMatchingPursuit00::UpdateBuf(int L0, int x0)
{
	int L, x, xx;
	double temp;

	double d0 = m_pDictBuf[L0][x0];

	for (L = 1; L <= m_nLevel; L++)
	{
		for (x = m_nUBegin[L0][L]; x < m_nUEnd[L0][L]; x++)
		{
			temp = m_pUBuf[L0][L][x + m_nUOffset];
			xx = x + x0;
			while (xx < 0)
				xx += m_nN;
			while (xx >= m_nN)
				xx -= m_nN;

			m_pDictBuf[L][xx] -= d0 * temp;
		}
	}
}

// マッチング追跡の逆変換 (double* input[]はlevel+2個あること)
bool CMatchingPursuit00::IPursuit(double *output, double *input[],
								int inputN, int level, int mode)
{
	int i, j, L, x;

	m_nN = inputN;

	if (m_pBuf == NULL || level > m_nMaxLevel || m_nN > m_nMaxN)
		return false;

	for (i = 0; i < m_nN; i++)
	{
		m_pDictBuf[0][i] = input[0][i] + input[level + 1][i];
	}

	for (L = 1; L <= level; L++)
	{
		for (i = 0; i < m_nN; i++)
		{
			for (j = m_nBBegin[L]; j < m_nBEnd[L]; j++)
			{
				x = i + j;
				while (x < 0)
					x += m_nN;
				while (x >= m_nN)
					x -= m_nN;
				m_pDictBuf[0][x] += input[L][i] * ReadB(L, j);
			}
		}
	}

	if (mode == 0)
	{
		for (i = 0; i < m_nN; i++)
			output[i] = m_pDictBuf[0][i];
	}
	else
	{
		m_pInter->IInter(output, m_pDictBuf[0], m_nN,
						m_pScal, m_nScalN, m_nScalOffset, 1);
	}

	return true;
}

// 計算後に真のεを確認する関数
double CMatchingPursuit00::Epsilon()
{
	return m_nEpsilon;
}
