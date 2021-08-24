#include <math.h>
#include "OrthoWavelet00.h"
#include "Matrix.h"

COrthoWavelet00::COrthoWavelet00()
{
    m_pBuf0 = NULL;
    m_pBuf1 = NULL;
    m_pScalingBuf = NULL;
    m_pMotherBuf = NULL;

    m_nN = 2;
    m_nLevel = 8;

    m_nPBegin = 0;
    m_nPEnd = 0;
    m_nQBegin = 0;
    m_nQEnd = 0;

    m_nScalingBegin = 0;
    m_nScalingEnd = 0;
    m_nScalingDenomi = 256;
    m_nScalingSz = 0;

    m_nMotherBegin = 0;
    m_nMotherEnd = 0;
    m_nMotherDenomi = 512;
    m_nMotherSz = 0;
}

COrthoWavelet00::~COrthoWavelet00()
{
    if (m_pBuf0 != NULL)
    {
        delete m_pBuf0;
    }
    if (m_pBuf1 != NULL)
    {
        delete m_pBuf1;
    }
}

double COrthoWavelet00::Scaling(int level, int n)
{
    int lv = m_nSLevel - level;

    if (lv < 0)
    {
        lv = 0;
    }
    if (lv > m_nSLevel)
    {
        lv = m_nSLevel;
    }

    n = n << lv;

    if (m_pScalingBuf == NULL)
    {
        return 0.0;
    }

    if (n < m_nScalingBegin || n >= m_nScalingEnd)
    {
        return 0.0;
    }
    else
    {
        return m_pScalingBuf[n - m_nScalingBegin];
    }
}

int COrthoWavelet00::ScalingBegin(int level)
{
    int lv = m_nSLevel - level;

    if (lv < 0)
    {
        lv = 0;
    }
    if (lv > m_nSLevel)
    {
        lv = m_nSLevel;
    }

    return m_nScalingBegin >> lv;
}

int COrthoWavelet00::ScalingEnd(int level)
{
    int lv = m_nSLevel - level;

    if (lv < 0)
    {
        lv = 0;
    }
    if (lv > m_nSLevel)
    {
        lv = m_nSLevel;
    }

    return m_nScalingEnd >> lv;
}

double COrthoWavelet00::Mother(int level, int n)
{
    int lv = m_nMLevel - level;

    if (lv < 0)
    {
        lv = 0;
    }
    if (lv > m_nMLevel)
    {
        lv = m_nMLevel;
    }

    n = n << lv;

    if (m_pMotherBuf == NULL)
    {
        return 0.0;
    }

    if (n < m_nMotherBegin || n >= m_nMotherEnd)
    {
        return 0.0;
    }
    else
    {
        return m_pMotherBuf[n - m_nMotherBegin];
    }
}

int COrthoWavelet00::MotherBegin(int level)
{
    int lv = m_nLevel - level;

    if (lv < 0)
    {
        lv = 0;
    }
    if (lv > m_nMLevel)
    {
        lv = m_nMLevel;
    }

    return m_nMotherBegin >> lv;
}

int COrthoWavelet00::MotherEnd(int level)
{
    int lv = m_nLevel - level;

    if (lv < 0)
    {
        lv = 0;
    }
    if (lv > m_nMLevel)
    {
        lv = m_nMLevel;
    }

    return m_nMotherEnd >> lv;
}

//分解数列H0
double COrthoWavelet00::H0(int n)
{
    if (-n < m_nPBegin || -n > m_nPEnd)
    {
        return 0.0;
    }
    else
    {
        return m_pP0[-n - m_nPBegin] / sqrt(2.0);
    }
}

//分解数列H1
double COrthoWavelet00::H1(int n)
{
    if (-n < m_nQBegin || -n > m_nQEnd)
    {
        return 0.0;
    }
    else
    {
        return m_pQ0[-n - m_nQBegin] / sqrt(2.0);
    }
}

//再構成数列G0
double COrthoWavelet00::G0(int n)
{
    if (n < m_nPBegin || n > m_nPEnd)
    {
        return 0.0;
    }
    else
    {
        return m_pP0[n - m_nPBegin] / sqrt(2.0);
    }
}

//再構成数列G1
double COrthoWavelet00::G1(int n)
{
    if (n < m_nQBegin || n > m_nQEnd)
    {
        return 0.0;
    }
    else
    {
        return m_pQ0[n - m_nQBegin] / sqrt(2.0);
    }
}

bool COrthoWavelet00::Prepare(int N, int level)
{
    CMatrix matrix;

    m_nMtxMaxN = matrix.MaxN();

    if (level < 0 || level > 12)
    {
        return false;
    }

    m_nN = N;
    m_nLevel = level;

    if (m_pBuf0 != NULL)
    {
        delete m_pBuf0;
    }
    m_pBuf0 = new double[(0x01 << (level + 1)) * m_nMtxMaxN * 2 + m_nMtxMaxN * 3];

    if (m_pBuf0 != NULL)
    {
        m_pPBuf0 = m_pBuf0;
        m_pPBuf1 = &m_pBuf0[(0x01 << (level + 1)) * m_nMtxMaxN];
        m_pP0 = &m_pPBuf1[(0x01 << (level + 1)) * m_nMtxMaxN];
        m_pQ0 = &m_pP0[m_nMtxMaxN];
        m_pSc0 = &m_pQ0[m_nMtxMaxN];
    }
    else
    {
        m_pScalingBuf = NULL;
        m_pMotherBuf = NULL;
        return false;
    }

    for (int i = 0; i < m_nMtxMaxN; i++)
    {
        m_pP0[i] = 0.0;
        m_pQ0[i] = 0.0;
    }

    PrepareSub(N);

    if (matrix.InitCk() == true)
    {
        for (int i = 0; i < m_nMtxMaxN; i++)
        {
            for (int j = 0; j < m_nMtxMaxN; j++)
            {
                if (i * 2 + 1 - j >= 0 && i * 2 + 1 - j < m_nMtxMaxN)
                {
                    matrix.m_nInpMatrix[i][j] = m_pP0[i * 2 + 1 - j];
                    if (i == j)
                    {
                        matrix.m_nInpMatrix[i][j] -= 1.0;
                    }
                }
                else
                {
                    matrix.m_nInpMatrix[i][j] = 0.0;
                }
            }
        }

        for (int i = 0; i < m_nMtxMaxN; i++)
        {
            matrix.m_nInpMatrix[0][i] = 1.0;
        }

        matrix.InverseMatrix(m_nPNum - 2);
    }
    else
    {
        return false;
    }

    for (int i = 0; i < m_nMtxMaxN; i++)
    {
        m_pSc0[i] = 0.0;
    }

    for (int i = 1; i < m_nPNum - 1; i++)
    {
        m_pSc0[i] = matrix.m_nOutMatrix[i - 1][0];
    }

    m_nPBegin = 0;
    m_nPEnd = m_nPNum - 1;

    m_nQBegin = -m_nPEnd + 1;
    m_nQEnd = 1;

    for (int i = 0; i <= m_nPEnd; i++)
    {
        m_pPBuf0[i] = m_pP0[i];
    }
    int area = m_nPEnd;

    for (int i = 2; i <= m_nLevel; i++)
    {
        for (int k = 0; k <= (area * 2 + m_nPEnd); k++)
        {
            double temp = 0.0;

            for (int L = 0; L <= area; L++)
            {
                if (k - L * 2 >= 0 && k - L * 2 <= m_nPEnd)
                {
                    temp += m_pP0[k - L * 2] * m_pPBuf0[L];
                }
            }
            m_pPBuf1[k] = temp;
        }
        area = area * 2 + m_nPEnd;

        for (int j = 0; j <= area; j++)
        {
            m_pPBuf0[j] = m_pPBuf1[j];
        }
    }

    //スケーリング関数の領域は計算の都合上0からにしておく
    //最後に本来の範囲に修正する
    m_nScalingDenomi = 0x01 << m_nLevel;
    m_nScalingSz = m_nScalingDenomi * m_nPEnd;
    m_nScalingBegin = 0;
    m_nScalingEnd = m_nScalingDenomi * m_nPEnd;

    m_nMotherDenomi = 0x01 << (m_nLevel + 1);
    m_nMotherSz = m_nMotherDenomi * m_nPEnd;
    m_nMotherBegin = m_nMotherDenomi * (1 - m_nPEnd) / 2;
    m_nMotherEnd = m_nMotherDenomi * (1 + m_nPEnd) / 2;

    if (m_pBuf1 != NULL)
    {
        m_pScalingBuf = NULL;
        m_pMotherBuf = NULL;
        delete m_pBuf1;
    }
    m_pBuf1 = new double[m_nScalingSz + m_nMotherSz];
    if (m_pBuf1 == NULL)
    {
        return false;
    }

    m_pScalingBuf = m_pBuf1;
    m_pMotherBuf = &m_pScalingBuf[m_nScalingSz];

    if (m_nLevel > 0)
    {
        for (int i = 0; i < m_nScalingSz; i++)
        {
            double temp = 0.0;

            for (int k = 0; k <= area; k++)
            {
                if (i - k >= 0 && i - k <= m_nPEnd)
                {
                    temp += m_pPBuf0[k] * m_pSc0[i - k];
                }
            }
            m_pScalingBuf[i] = temp;
        }
    }
    else
    {
        for (int i = 0; i < m_nScalingSz; i++)
        {
            m_pScalingBuf[i] = m_pSc0[i];
        }
    }

    for (int i = m_nQBegin; i <= m_nQBegin + m_nPEnd; i++)
    {
        double temp = m_pP0[1 - i];
        if (i & 0x01)
        {
            temp = -temp;
        }
        m_pQ0[i - m_nQBegin] = -temp;
    }

    for (int i = m_nMotherBegin; i < m_nMotherEnd; i++)
    {
        double temp = 0.0;
        for (int k = m_nQBegin; k <= m_nQBegin + m_nPEnd; k++)
        {
            int x = 2 * i / 2 - k * m_nScalingDenomi;
            if (x >= 0 && x < m_nScalingSz)
            {
                temp += m_pQ0[k - m_nQBegin] * m_pScalingBuf[x];
            }
        }
        m_pMotherBuf[i - m_nMotherBegin] = temp;
    }

    m_nPBegin += m_nPShift;
    m_nPEnd += m_nPShift;
    m_nQBegin -= m_nPShift;
    m_nQEnd -= m_nPShift;

    //スケーリング関数の領域を修正
    m_nScalingBegin = m_nScalingDenomi * m_nPBegin;
    m_nScalingEnd = m_nScalingDenomi * m_nPEnd;

    m_nSLevel = m_nLevel;
    m_nMLevel = m_nLevel + 1;

    return true;
}

void COrthoWavelet00::PrepareSub(int N)
{
    switch (N)
    {
    case 11:
        m_nPShift = -6; // Coiflet 6
        m_nPNum = 18;
        m_pP0[0] = -0.005364837342;
        m_pP0[1] = 0.011006253418;
        m_pP0[2] = 0.033167120958;
        m_pP0[3] = -0.093015528958;
        m_pP0[4] = -0.086441527120;
        m_pP0[5] = 0.573006670548;
        m_pP0[6] = 1.122570513740;
        m_pP0[7] = 0.605967143546;
        m_pP0[8] = -0.101540281510;
        m_pP0[9] = -0.116392501524;
        m_pP0[10] = 0.048868188642;
        m_pP0[11] = 0.022458481924;
        m_pP0[12] = -0.012739202022;
        m_pP0[13] = -0.003640917832;
        m_pP0[14] = 0.001580410202;
        m_pP0[15] = 0.000659330348;
        m_pP0[16] = -0.000100385550;
        m_pP0[17] = -0.000048931468;
        break;
    case 22:
        m_nPShift = -30; // O Spline 3
        m_nPNum = 61;
        m_pP0[0] = -0.00000757744;
        m_pP0[1] = 0.00001144676;
        m_pP0[2] = 0.00001460739;
        m_pP0[3] = -0.00002222839;
        m_pP0[4] = -0.00002821716;
        m_pP0[5] = 0.00004330400;
        m_pP0[6] = 0.00005463413;
        m_pP0[7] = -0.00008468218;
        m_pP0[8] = -0.00010606379;
        m_pP0[9] = 0.00016635055;
        m_pP0[10] = 0.00020653439;
        m_pP0[11] = -0.00032858869;
        m_pP0[12] = -0.00040359353;
        m_pP0[13] = 0.00065352962;
        m_pP0[14] = 0.00079187000;
        m_pP0[15] = -0.00131125702;
        m_pP0[16] = -0.00156092382;
        m_pP0[17] = 0.00266173876;
        m_pP0[18] = 0.00309247829;
        m_pP0[19] = -0.00549057847;
        m_pP0[20] = -0.00615725881;
        m_pP0[21] = 0.01159864030;
        m_pP0[22] = 0.01228286172;
        m_pP0[23] = -0.02543084221;
        m_pP0[24] = -0.02429097832;
        m_pP0[25] = 0.05949363315;
        m_pP0[26] = 0.04536924030;
        m_pP0[27] = -0.15561584377;
        m_pP0[28] = -0.07099595988;
        m_pP0[29] = 0.61365927344;
        m_pP0[30] = 1.08347151257;
        m_pP0[31] = 0.61365927344;
        m_pP0[32] = -0.07099595988;
        m_pP0[33] = -0.15561584377;
        m_pP0[34] = 0.04536924030;
        m_pP0[35] = 0.05949363315;
        m_pP0[36] = -0.02429097832;
        m_pP0[37] = -0.02543084221;
        m_pP0[38] = 0.01228286172;
        m_pP0[39] = 0.01159864030;
        m_pP0[40] = -0.00615725881;
        m_pP0[41] = -0.00549057847;
        m_pP0[42] = 0.00309247829;
        m_pP0[43] = 0.00266173876;
        m_pP0[44] = -0.00156092382;
        m_pP0[45] = -0.00131125702;
        m_pP0[46] = 0.00079187000;
        m_pP0[47] = 0.00065352962;
        m_pP0[48] = -0.00040359353;
        m_pP0[49] = -0.00032858869;
        m_pP0[50] = 0.00020653439;
        m_pP0[51] = 0.00016635055;
        m_pP0[52] = -0.00010606379;
        m_pP0[53] = -0.00008468218;
        m_pP0[54] = 0.00005463413;
        m_pP0[55] = 0.00004330400;
        m_pP0[56] = -0.00002821716;
        m_pP0[57] = -0.00002222839;
        m_pP0[58] = 0.00001460739;
        m_pP0[59] = 0.00001144676;
        m_pP0[60] = -0.00000757744;
        break;
    }
}
