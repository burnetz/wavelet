#include <math.h>
#include "Hilbert00.h"
#include "Fft00.h"
#include <stdio.h>
#include <string.h>

#define _PI 3.14159265358979323846264338327950288419716939937510

CHilbert00::CHilbert00()
{
    m_pFft = new CFft00;

    for (int i = 0; i < 2; i++)
    {
        m_pWork[i] = NULL;
    }
}

CHilbert00::~CHilbert00()
{
    if (m_pFft != NULL)
    {
        delete m_pFft;
    }

    for (int i = 0; i < 2; i++)
    {
        if (m_pWork[i] != NULL)
        {
            delete m_pWork[i];
        }
    }
}

bool CHilbert00::Prepare(int maxN)
{
    int i;
    for (i = 0; i < 32; i++)
    {
        //4096は巡回畳み込み緩衝のため
        if ((0x01 << i) > (maxN + 4096))
            break;
    }

    m_nMaxEx = i;
    //printf("%d\n", m_nMaxEx);
    m_nMaxN = 0x01 << m_nMaxEx;
    //printf("%d\n", m_nMaxN);

    if (m_pFft->Prepare(m_nMaxEx) == false)
    {
        return false;
    }

    for (i = 0; i < 2; i++)
    {
        if (m_pWork[i] != NULL)
        {
            delete m_pWork[i];
        }
        m_pWork[i] = NULL;
    }
    for (i = 0; i < 2; i++)
    {
        m_pWork[i] = new double[m_nMaxN];
        if (m_pWork[i] == NULL)
        {
            break;
        }
    }

    if (i < 2)
    {
        for (i = 0; i < 2; i++)
        {
            if (m_pWork[i] != NULL)
            {
                delete m_pWork[i];
            }
            m_pWork[i] = NULL;
        }

        return false;
    }
    m_pReal = m_pWork[0];
    m_pImag = m_pWork[1];

    return true;
}

bool CHilbert00::Hilbert(double *outputReal, double *outputImag,
             double *input, int inputN, double freq0, double freq1)
{
    int i;

    if(inputN > m_nMaxN || m_pWork[0] == NULL){
        printf("error in hilbert 0\n");
        printf("inputN %d m_nMaxN %d\n", inputN, m_nMaxN);
        return false;
    }

    if(freq1 < freq0){
        double d = freq0;
        freq0 = freq1;
        freq1 = d; //多分誤植。やりたいことはおそらくスワップ
    }
    if(freq0 < 0.0)freq0 = 0.0;
    if(freq0 > 0.5)freq0 = 0.5;
    if(freq1 < 0.0)freq1 = 0.0;
    if(freq1 > 0.5)freq1 = 0.5;

    for(i = 0; i < inputN; i++){
        m_pReal[i] = input[i];
    }
    for(i = inputN; i < m_nMaxN; i++){
        m_pReal[i] = 0.0;
    }
    for(i = 0; i < m_nMaxN; i++){
        m_pImag[i] = 0.0;
    }

    //FFT
    m_pFft->Fft(m_pReal, m_pImag, 1);

    //周波数0~0.5を2倍, DCはそのまま
    for(i = 1; i < m_nMaxN/2; i++){
        m_pReal[i] *= 2.0;
        m_pImag[i] *= 2.0;
    }
    for(i = m_nMaxN/2; i < m_nMaxN; i++){
        m_pReal[i] = 0.0;
        m_pImag[i] = 0.0;
    }

    //帯域制限（freq0 = 0, freq1 = 0.5の場合は影響なし）
    for(i = 0; i < m_nMaxN*freq0; i++){
        m_pReal[i] = 0.0;
        m_pImag[i] = 0.0;
    }
    for(i = (int)(m_nMaxN*freq1); i < m_nMaxN/2; i++){
        m_pReal[i] = 0.0;
        m_pImag[i] = 0.0;
    }
    //IFFT
    m_pFft->Fft(m_pReal, m_pImag, -1);

    for(i = 0; i < inputN; i++){
        outputReal[i] = m_pReal[i];
        outputImag[i] = m_pImag[i];
    }

    return true;

}

//瞬時パワー
bool CHilbert00::InstantPower(double* output, double* inputReal, double* inputImag, int inputN){
    int i;

    for(i = 0; i < inputN; i++){
        output[i] = inputReal[i]*inputReal[i] + inputImag[i]*inputImag[i];
    }

    return true;
}

//瞬時周波数
bool CHilbert00::InstantFreq(double* output, double* inputReal, double* inputImag, int inputN){
    int i;

    for(i = 0; i < inputN - 1; i++){
        double phase0 = atan2(inputImag[i], inputReal[i]);
        double phase1 = atan2(inputImag[i + 1], inputReal[i + 1]);
        double omega = phase1 - phase0;

        while(omega <= -_PI){
            omega += 2*_PI;
        }
        while(omega > _PI){
            omega -= 2*_PI;
        }

        output[i] = omega/(2*_PI);

    }

    return true;
}

#ifdef HILBERT_TEST

int main(int argc, char* argv[])
{
    int n;

    if(argc != 2){
        return 0;
    }
    scanf("%d", &n);

    double input[n];
    for (int i = 0; i < n; i++){
        scanf("%lf", &input[i]);
    }

    CHilbert00* chilb = new CHilbert00();

    if(!chilb->Prepare(n)){
        printf("Can't prepare!\n");
        return 0;
    }

    double oreal[n], oimag[n];

    if(!chilb->Hilbert(oreal, oimag, input, n, 0.0, 0.5)){
        printf("failed!\n");
        return 0;
    }

    for(int i = 0; i < n; i++){
        if(strcmp(argv[1], "real") == 0){
            printf("%d %f\n", i, oreal[i]);
        }
        if(strcmp(argv[1], "imag") == 0){
            printf("%d %f\n", i, oimag[i]);
        }
    }
}

#endif