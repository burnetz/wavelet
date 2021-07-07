#include <math.h>
#include "FftAnalyze00.h"
#include "Fft00.h"

#define _PI 3.14159265358979323846264338327950288419716939937510

CFftAnalyze00::CFftAnalyze00()
{
    m_pFft = new CFft00;

    for (int i = 0; i < 16; i++)
    {
        m_pWork[i] = NULL;
    }
}

CFftAnalyze00::~CFftAnalyze00()
{
    if (m_pFft != NULL)
    {
        delete m_pFft;
    }

    for (int i = 0; i < 16; i++)
    {
        if (m_pWork[i] != NULL)
        {
            delete m_pWork[i];
        }
    }
}

bool CFftAnalyze00::Prepare(int ex, int winType)
{
    m_nEx = ex;
    m_nFftLen = 0x01 << m_nEx;

    if (m_pFft->Prepare(m_nEx) == false)
    {
        return false;
    }

    for (int i = 0; i < 3; i++)
    {
        if (m_pWork[i] != NULL)
        {
            delete m_pWork[i];
        }
        m_pWork[i] = NULL;
    }

    int i;
    for (i = 0; i < 3; i++)
    {
        m_pWork[i] = new double[m_nFftLen];
        if (m_pWork[i] == NULL)
        {
            break;
        }
    }

    if (i < 3)
    {
        for (i = 0; i < 3; i++)
        {
            if (m_pWork[i] != NULL)
            {
                delete m_pWork[i];
            }
            m_pWork[i] = NULL;
        }
        return false;
    }

    m_pWindow = m_pWork[0];
    m_pReal = m_pWork[1];
    m_pImag = m_pWork[2];

    Window(m_pWindow, m_nFftLen, winType, m_nScale);

    return true;
}

void CFftAnalyze00::Window(double *output, int winLen, int winType, double &scale)
{
    scale = 0.0;

    for (int i = 0; i < winLen; i++)
    {
        switch (winType)
        {
        case 00: //ハニング
            output[i] = .5 - .5 * cos(2 * _PI * i / winLen);
            break;
        case 01: //ハミング
            output[i] = .54 - .46 * cos(2 * _PI * i / winLen);
            break;
        case 02: //ブラックマン
            output[i] = .42 - .5 * cos(2 * _PI * i / winLen) + .08 * cos(4*_PI*i / winLen);
            break;
        case 03: //ガウス
            output[i] = exp(-4.5*(2.0*i - winLen)/winLen * (2.0*i - winLen)/winLen);
            break;
        case 04: //方形
            output[i] = 1.0;
            break;
        }

        scale += output[i]*output[i];
    }

    scale /= winLen;
    scale = 1.0/sqrt(scale);
}

bool CFftAnalyze00::Fft(double* outAmp, double* outPhase, double* input){
    if(m_pWork[0] == NULL){
        return false;
    }

    for(int i = 0; i < m_nFftLen; i++){
        m_pReal[i] = input[i];
        m_pImag[i] = 0.0;
    }

    for(int i = 0; i < m_nFftLen; i++){
        m_pReal[i] *= m_pWindow[i];
        m_pImag[i] *= m_pWindow[i];
    }

    m_pFft->Fft(m_pReal, m_pImag, 1);

    for(int i = 0; i < m_nFftLen; i++){
        m_pReal[i] *= m_nScale;
        m_pImag[i] *= m_nScale;
    }

    for(int i = 0; i < m_nFftLen; i++){
        outAmp[i] = sqrt(m_pReal[i]*m_pReal[i] + m_pImag[i]*m_pImag[i]);
        outPhase[i] = atan2(m_pImag[i], m_pReal[i]);
    }

    return true;
}