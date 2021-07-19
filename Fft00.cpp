#include <math.h>
#include "Fft00.h"
#include <stdio.h>

#define _PI 3.14159265358979323846264338327950288419716939937510

CFft00::CFft00()
{
    m_pBuf = NULL;
    m_nEx = -1;
    m_nLength = -1;
    m_nBufLen = -1;
}

CFft00::~CFft00()
{
    if (m_pBuf != NULL)
    {
        delete m_pBuf;
    }
}

//ex: サンプル数の指数（N=2^ex）
bool CFft00::Prepare(int ex)
{
    int len;

    if (m_nEx == ex && m_pBuf != NULL)
    {
        return true;
    }

    m_nEx = ex;
    len = 0x01 << ex;
    if (len > m_nBufLen)
    {
        m_nBufLen = 0;
        if (m_pBuf != NULL)
        {
            delete m_pBuf;
        }
        m_pBuf = new double[len * 5 + 256];
        if (m_pBuf == NULL)
        {
            return false;
        }
        m_nBufLen = len;
    }
    m_nLength = len;

    m_pFftSin = m_pBuf;
    m_pFftCos = m_pFftSin + m_nLength;
    m_pIfftSin = m_pFftCos + m_nLength;
    m_pIfftCos = m_pIfftSin + m_nLength;
    m_pSubBuf = m_pIfftCos + m_nLength;

    for (int i = 0; i < m_nLength; i++)
    {
        //本とは符号が逆だが、こうするとwolfram等のデフォルトと同じ
        //ただしwikipedia等の逆
        m_pFftSin[i] = sin(_PI * 2 * i / m_nLength);
        m_pFftCos[i] = cos(_PI * 2 * i / m_nLength);
        m_pIfftSin[i] = sin(-_PI * 2 * i / m_nLength);
        m_pIfftCos[i] = cos(-_PI * 2 * i / m_nLength);
    }

    return true;
}

//real: 入力信号の実数部
//imag: 入力信号の虚数部
//inv: 1でFFT, -1でIFFTの計算を行う
bool CFft00::Fft(double *real, double *imag, int inv)
{
    double *sinTable;
    double *cosTable;
    if (inv == 1)
    {
        sinTable = m_pFftSin;
        cosTable = m_pFftCos;
    }
    else
    {
        sinTable = m_pIfftSin;
        cosTable = m_pIfftCos;
    }

    int count = 1;

    int length2 = m_nLength;

    for (int i = 0; i < m_nEx; i++)
    {
        length2 >>= 1;
        int temp = 0;
        for (int j = 0; j < count; j++)
        {
            int w = 0;
            for (int k = 0; k < length2; k++)
            {
                int j1 = temp + k;
                int j2 = j1 + length2;
                double ar = real[j1];
                double ai = imag[j1];
                double br = real[j2];
                double bi = imag[j2];
                real[j1] = ar + br;
                imag[j1] = ai + bi;
                ar = ar - br;
                ai = ai - bi;
                real[j2] = ar * cosTable[w] - ai * sinTable[w];
                imag[j2] = ar * sinTable[w] + ai * cosTable[w];
                w += count;
            }
            temp += (2 * length2);
        }
        count <<= 1;
    }
    BitReverce(m_pSubBuf, real, m_nEx);
    BitReverce(m_pSubBuf, imag, m_nEx);

    if (inv != 1)
    {
        double nrml = (double)m_nLength;

        for (int i = 0; i < m_nLength; i++)
        {
            real[i] /= nrml;
            imag[i] /= nrml;
        }
    }

    return true;
}

//real: 入力信号の実数部
//imag: 入力信号の虚数部
//inv: 1でFFT, -1でIFFTの計算を行う
bool CFft00::Dft(double *real, double *imag, int inv)
{
    double *sinTable;
    double *cosTable;
    if (inv == 1)
    {
        sinTable = m_pFftSin;
        cosTable = m_pFftCos;
    }
    else
    {
        sinTable = m_pIfftSin;
        cosTable = m_pIfftCos;
    }

    double* tmpReal = new double[m_nLength];
    double* tmpImag = new double[m_nLength];
    for(int t = 0; t < m_nLength; t++){
        tmpReal[t] = tmpImag[t] = 0;
        for(int x = 0; x < m_nLength; x++){
            int w = t*x%m_nLength;
            tmpReal[t] += real[x]*cosTable[w] - imag[x]*sinTable[w];
            tmpImag[t] += real[x]*sinTable[w] + imag[x]*cosTable[w];
        }
    }

    for(int i = 0; i < m_nLength; i++){
        real[i] = tmpReal[i];
        imag[i] = tmpImag[i];
    }
    if (inv != 1)
    {
        double nrml = (double)m_nLength;

        for (int i = 0; i < m_nLength; i++)
        {
            real[i] /= nrml;
            imag[i] /= nrml;
        }
    }

    return true;
}

void CFft00::BitReverce(double *buf, double *a, int ex)
{
    int length = 1 << ex;

    for (int i = 0; i < length; i++)
    {
        int bit = 0;
        for (int j = 0; j < ex; j++)
        {
            if (i & (1 << j))
            {
                bit |= (1 << (ex - j - 1));
            }
        }
        buf[i] = a[bit];
    }
    for (int i = 0; i < length; i++)
    {
        a[i] = buf[i];
    }
}

#ifdef FFT00_TEST 

int main(){
    double input[8] = {1,1,0.5,0.5,0,0,0,0};

    int n = 8;

    CFft00* fft = new CFft00;

    fft->Prepare(3);

    double imag[8] = {0,0,0,0,0,0,0,0};//{0,0,1,1,0,0,-1,-1};
    fft->Fft(input, imag, 1);

    printf("FFT : \n");
    for(int i = 0; i < 8; i++){
        printf("%f %f\n", input[i], imag[i]);
    }

    fft->Fft(input, imag, -1);

    printf("IFFT : \n");
    for(int i = 0; i < 8; i++){
        printf("%f %f\n", input[i], imag[i]);
    }

    return 0;
}

#endif

#ifdef DFT00_TEST 

int main(){
    double input[8] = {1,1,0.5,0.5,0,0,0,0};

    int n = 8;

    CFft00* fft = new CFft00;

    fft->Prepare(3);

    double imag[8] = {0,0,1,1,0,0,-1,-1};
    fft->Dft(input, imag, 1);

    printf("DFT : \n");
    for(int i = 0; i < 8; i++){
        printf("%f %f\n", input[i], imag[i]);
    }

    fft->Dft(input, imag, -1);

    printf("IDFT : \n");
    for(int i = 0; i < 8; i++){
        printf("%f %f\n", input[i], imag[i]);
    }

    return 0;
}

#endif