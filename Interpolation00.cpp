#include <math.h>
#include "Interpolation00.h"
#include "Fft00.h"

#define FFTBUFEX 9

CInterpolation00::CInterpolation00(){
    m_pFft = new CFft00;

    m_pBuf = NULL;
    m_pBuf0 = NULL;
    m_pBuf1 = NULL;
    m_pBuf2 = NULL;
    m_pBuf3 = NULL;
    m_pBuf4 = NULL;
    m_pBuf5 = NULL;
    m_pBuf6 = NULL;
}

CInterpolation00::~CInterpolation00(){
    if(m_pFft != NULL){
        delete m_pFft;
    }

    if(m_pBuf != NULL){
        delete m_pBuf;
    }
}

bool CInterpolation00::Prepare(){
    int fftbufn = 0x01 << FFTBUFEX;

    if(m_pBuf != NULL){
        delete m_pBuf;
    }
    m_pBuf = new double[fftbufn * 7 + 256];

    if(m_pBuf != NULL){
        m_pBuf0 = m_pBuf;
        m_pBuf1 = m_pBuf0 + fftbufn;
        m_pBuf2 = m_pBuf1 + fftbufn;
        m_pBuf3 = m_pBuf2 + fftbufn;
        m_pBuf4 = m_pBuf3 + fftbufn;
        m_pBuf5 = m_pBuf4 + fftbufn;
        m_pBuf6 = m_pBuf5 + fftbufn;
    }
    else{
        return false;
    }

    if(m_pFft->Prepare(FFTBUFEX) == false){
        return false;
    }

    return true;
}

//補間
bool CInterpolation00::Inter(double* output, double* input, int inputN, double* scal, int scalN, int scalOffset, int mode){
    int fftbufn = 0x01 << FFTBUFEX;
    int fftn = fftbufn - scalN*2;

    if(m_pBuf != NULL || fftn < 0){
        return false;
    }

    for(int i = 0; i < fftbufn; i++){
        m_pBuf5[i] = 0.0;
    }

    int n = 0;
    int n1;

    while(n < inputN){
        if(fftn < inputN - n){
            n1 = fftn;
        }
        else{
            n1 = inputN - n;
        }

        //データをバッファに
        for(int i = 0; i < n1; i++){
            m_pBuf0[i] = input[n + i];
        }
        for(int i = n1; i < fftbufn; i++){
            m_pBuf0[i] = 0.0;
        }
        for(int i = 0; i < fftbufn; i++){
            m_pBuf1[i] = 0.0;
        }
        for(int i = 0; i < fftbufn; i++){
            m_pBuf2[i] = 0.0;
        }
        for(int i = 0; i < (scalN - scalOffset); i++){
            m_pBuf2[i] = scal[scalOffset + i];
        }
        for(int i = 0; i < scalOffset; i++){
            m_pBuf2[fftbufn - scalOffset + i] = scal[i];
        }

        for(int i = 0; i < fftbufn; i++){
            m_pBuf3[i] = 0.0;
        }

        //FFT
        if(!m_pFft->Fft(m_pBuf0, m_pBuf1, 1)){
            return false;
        }

        if(!m_pFft->Fft(m_pBuf2, m_pBuf3, 1)){
            return false;
        }

        //周波数領域の操作
        for(int i = 0; i < fftbufn; i++){
            double temp1 = m_pBuf2[i]*m_pBuf2[i] + m_pBuf3[i]*m_pBuf3[i];
            double temp2 = m_pBuf0[i]*m_pBuf2[i] + m_pBuf1[i]*m_pBuf3[i];
            double temp3 = m_pBuf1[i]*m_pBuf2[i] + m_pBuf0[i]*m_pBuf3[i];

            if(temp1 == 0.0){
                return false;
            }

            m_pBuf0[i] = temp2/temp1;
            m_pBuf1[i] = temp3/temp1;
        }

        //IFFT
        if(!m_pFft->Fft(m_pBuf0, m_pBuf1, -1)){
            return false;
        }

        if(n > 0){
            for(int i = 0; i < scalN; i++){
                if(n - scalN + i >= 0){
                    output[n - scalN + i] += m_pBuf0[fftbufn - scalN + i];
                }
            }
        }
        //信号をループとみなす処理のための保存
        else{
            for(int i = 0; i < fftbufn; i++){
                m_pBuf6[i] = m_pBuf0[i];
            }
        }

        for(int i = 0; i < scalN; i++){
            m_pBuf0[i] += m_pBuf5[fftn + i];
        }

        for(int i = 0; i < fftbufn; i++){
            m_pBuf5[i] = m_pBuf0[i];
        }

        for(int i = 0; i < n1; i++){
            output[n + i] = m_pBuf0[i];
        }
        n += n1;
    }

    //信号をループとみなす処理
    if(mode == 1){
        for(int i = 0; i < scalN; i++){
            if(i < inputN){
                output[i] += m_pBuf0[n1 + i];
            }
        }

        for(int i = 0; i < scalN ; i++){
            if(inputN - scalN + i >= 0){
                output[inputN - scalN + i] += m_pBuf6[fftbufn - scalN + i];
            }
        }
    }

    return true;
}

bool CInterpolation00::IInter(double* output, double* input, int inputN, double* scal, int scalN, int scalOffset, int mode){
    //信号をループとみなす
    if(mode == 1){
        //畳み込みの計算
        for(int i = 0; i < inputN; i++){
            double temp = 0.0;
            for(int j = -scalOffset; j < scalN - scalOffset; j++){
                int x = i - j;
                while(x < 0){
                    x += inputN;
                }
                while(x >= inputN){
                    x -= inputN;
                }

                temp += scal[j + scalOffset] * input[x];
            }
            output[i] = temp;
        }
    }
    else{
        //畳み込みの計算
        for(int i = 0; i < inputN; i++){
            double temp = 0.0;
            for(int j = -scalOffset; j < scalN - scalOffset; j++){
                int x = i - j;
                if(x >= 0 && x < inputN){
                    temp += scal[j + scalOffset] * input[x];
                }
            }
            output[i] = temp;
        }
    }

    return true;
}