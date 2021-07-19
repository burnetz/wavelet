#include <math.h>
#include "Dwvd00.h"
#include "Hilbert00.h"

#define _PI 3.14159265358979323846264338327950288419716939937510

CDwvd00::CDwvd00(){
    for(int i = 0; i < 8; i++){
        m_pWork[i] = NULL;
    }
}

CDwvd00::~CDwvd00(){
    for(int i = 0; i < 8; i++){
        if(m_pWork[i] != NULL){
            delete m_pWork[i];
        }
    }
}

bool CDwvd00::Prepare(double* input, int inputN){
    for(int i = 0; i < 2; i++){
        if(m_pWork[i] != NULL){
            delete m_pWork[i];
        }
        m_pWork[i] = NULL;
    }

    m_pInput = input;
    m_nInputN = inputN;

    int i;
    for(i = 0; i < 2; i++){
        m_pWork[i] = new double[m_nInputN];
        if(m_pWork[i] == NULL){
            break;
        }
    }

    bool ret;
    if(i == 2){
        m_pReal = m_pWork[0];
        m_pImag = m_pWork[1];
        ret = true;
    }
    if(ret == false){
        for(i = 0; i < 2; i++){
            if(m_pWork[i] != NULL){
                delete m_pWork[i];
            }
            m_pWork[i] = NULL;
        }
        return ret;
    }

    //Hilbert
    CHilbert00 hilbert;

    if(ret == true){
        ret = hilbert.Prepare(m_nInputN);
    }
    if(ret == true){
        ret = hilbert.Hilbert(m_pReal, m_pImag, m_pInput, m_nInputN, 0.0, 0.5);
    }
    return ret;
}

bool CDwvd00::SetParameters(double sigma, double damp, double omega){
    bool ret = false;

    for(int i = 2; i < 4; i++){
        if(m_pWork[i] != NULL){
            delete m_pWork[i];
        }
        m_pWork[i] = NULL;
    }

    if(sigma <= 0.0 || damp <= 0.0 || damp >= 1.0 || omega <= 0.0){
        return ret;
    }

    m_nWinOffset = (int)sqrt(-2.0 * sigma * sigma * log(damp));
    m_nWinLength = m_nWinOffset*2 + 1;

    int i;
    for(i = 0; i < 4; i++){
        m_pWork[i] = new double[m_nWinLength];
        if(m_pWork[i] == NULL){
            break;
        }
    }

    if(i == 4){
        m_pWindowR = m_pWork[2];
        m_pWindowI = m_pWork[3];
        ret = true;
    }

    if(ret == true){
        //e^i w xに窓をかけて収納
        for(i = 0; i < m_nWinLength; i++){
            double x = (double)(i - m_nWinOffset);
            double d = exp(-x*x/(2.0*sigma*sigma))/sqrt(2.0*_PI*sigma*sigma);
            m_pWindowR[i] = d*cos(2.0*omega*x);
            m_pWindowI[i] = d*sin(2.0*omega*x);
        }
    }
    return ret;
}

//座標xのウィグナー・ビレ分布
double CDwvd00::Dwvd(int x){
    double real2, imag2, real4, imag4;

    double d = 0.0;
    for(int i = -m_nWinOffset; i <= m_nWinOffset; i++){
        double real0 = Real(x + i);
        double imag0 = Imag(x + i);
        double real1 = Real(x - i);
        double imag1 = Imag(x - i);

        ComplexMultiply(real2, imag2, real0, imag0, real1, -imag1);

        double real3 = ReadWindowR(i);
        double imag3 = ReadWindowI(i);

        ComplexMultiply(real4, imag4, real2, imag2, real3, -imag3);

        d += 2*real4;
    }

    return d;
}

//解析信号の実数部読み込み
double CDwvd00::Real(int x){
    if(x < 0 || x >= m_nInputN){
        return 0.0;
    }
    return m_pReal[x];
}

//解析信号の虚数部読み込み
double CDwvd00::Imag(int x){
    if(x < 0 || x >= m_nInputN){
        return 0.0;
    }
    return m_pImag[x];
}

//複素数sin * 窓の実数部読み込み
double CDwvd00::ReadWindowR(int x){
    x += m_nWinOffset;
    if(x < 0 || x >= m_nWinLength){
        return 0.0;
    }
    return m_pWindowR[x];
}

//複素数sin * 窓の虚数部読み込み
double CDwvd00::ReadWindowI(int x){
    x += m_nWinOffset;
    if(x < 0 || x >= m_nWinLength){
        return 0.0;
    }
    return m_pWindowI[x];
}

//複素数の掛け算
void CDwvd00::ComplexMultiply(double &rOutput, double &iOutput, double rInput0, double iInput0, double rInput1, double iInput1){
    rOutput = rInput0*rInput1 - iInput0*iInput1;
    iOutput = iInput0*rInput1 + rInput0*iInput1;
}

//平滑化の準備
bool CDwvd00::PrepareSp(double sigma, double damp){
    bool ret = false;

    if(m_pWork[4] != NULL){
        delete m_pWork[4];
    }
    m_pWork[4] = NULL;

    if(sigma < 0.0 || damp <= 0.0 || damp >= 1.0){
        return ret;
    }

    m_nSpGaussOffset = (int)sqrt(-2.0*sigma*sigma*log(damp));
    m_nSpGaussLength = m_nSpGaussOffset*2 + 1;

    m_pWork[4] = new double[m_nSpGaussLength];

    if(m_pWork[4] != NULL){
        m_pGauss = m_pWork[4];

        for(int i = 0; i < m_nSpGaussLength; i++){
            if(sigma > 0.0){
                double x = (double)(i - m_nSpGaussOffset);
                m_pGauss[i] = exp(-x*x/(2.0*sigma*sigma))/sqrt(2.0*_PI*sigma*sigma);
            }
            else{
                m_pGauss[i] = 1.0;
            }
        }
        ret = true;
    }
    return ret;
}

//平滑化
double CDwvd00::Sp(double* input, int inputN, int x){
    m_pSpInput = input;
    m_nSpInputN = inputN;

    double d = 0.0;
    for(int i = -m_nSpGaussOffset; i <= m_nSpGaussOffset; i++){
        d += SpReadGauss(i)*SpReadData(x - i);
    }
    return d;
}

//平滑化のサブルーチン
double CDwvd00::SpReadData(int x){
    if(x < 0){
        x = 0;
    }
    if(x >= m_nSpInputN){
        x = m_nSpInputN - 1;
    }

    return m_pSpInput[x];
}

//平滑化のガウス関数
double CDwvd00::SpReadGauss(int x){
    int n = x + m_nSpGaussOffset;
    if(n < 0 || n >= m_nSpGaussLength){
        return 0.0;
    }
    return m_pGauss[n];
}
