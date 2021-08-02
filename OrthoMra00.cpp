#include <math.h>
#include "OrthoMra00.h"
#include "OrthoWavelet00.h"
#include "Interpolation00.h"
#include <stdio.h>

#define MAXLEVEL 32

COrthoMra00::COrthoMra00(){
    m_pOrtho = new COrthoWavelet00;
    m_pInter = new CInterpolation00;

    m_pBuf = NULL;
    m_pBuf0 = NULL;
    m_pBuf1 = NULL;
}

COrthoMra00::~COrthoMra00(){
    if(m_pOrtho != NULL){
        delete m_pOrtho;
    }
    if(m_pInter != NULL){
        delete m_pInter;
    }
    if(m_pBuf != NULL){
        delete m_pBuf;
    }
}

bool COrthoMra00::Prepare(int wName, int maxN){
    m_nWName = wName;
    m_nMaxN = maxN;

    if(!m_pOrtho->Prepare(wName, 0)){
        return false;
    }
    
    if(!m_pInter->Prepare()){
        return false;
    }

    int scalBegin = m_pOrtho->ScalingBegin(0);
    int scalEnd = m_pOrtho->ScalingEnd(0);
    m_nScalOffset = -scalBegin;
    m_nScalN = scalEnd - scalBegin + 1;

    if(m_pBuf != NULL){
        delete m_pBuf;
    }
    m_pBuf = new double[maxN*2 + m_nScalN + 256];

    if(m_pBuf != NULL){
        m_pBuf0 = m_pBuf;
        m_pBuf1 = m_pBuf0 + maxN;
        m_pScal = m_pBuf1 + maxN;
    }
    else{
        return false;
    }

    for(int i = 0; i < m_nScalN; i++){
        m_pScal[i] = m_pOrtho->Scaling(0, i - m_nScalOffset);
    }

    return true;
}

//多重解像度解析
//output内に各レベルの解析結果が保存される
//内訳は例えばサンプル数512でレベル3まで解析する場合
//低域(64)、レベル-3(64)、レベル-2(128)、レベル-1(256)のようになる
bool COrthoMra00::Mra(double* output, double* input, int level, int mul, int mode){
    //int i, j, k, L, n;
    double sum;

    m_nN = (0x01 << level)*mul;
    m_nLevel = level;

    if(level > MAXLEVEL || m_nN > m_nMaxN || m_pBuf == NULL){
        return false;
    }

    //補間
    if(mode){
        if(!m_pInter->Inter(m_pBuf0, input, m_nN, m_pScal, m_nScalN, m_nScalOffset, 1)){
            printf("Interpolation Failed!\n");
        }
    }
    //補間を省略
    else{
        for(int i = 0; i < m_nN ; i++){
            m_pBuf0[i] = input[i];
        }
    }

    //分解アルゴリズムのループ
    for(int i = 1; i <= level; i++){
        int n = 0;
        int sz = m_nN >> i;

        int begin = m_pOrtho->H0Begin();
        int end = m_pOrtho->H0End();

        for(int k = 0; k < sz; k++){
            sum = 0;
            for(int L = 2*k - end; L <= 2*k - begin; L++){
                sum += m_pOrtho->H0(2*k - L)*Read(L, sz*2, 0);
            }
            m_pBuf1[n] = sum;
            n++;
        }

        begin = m_pOrtho->H1Begin();
        end = m_pOrtho->H1End();
        
        for(int k = 0; k < sz; k++){
            sum = 0;
            
            for(int L = 2*k - end; L <= 2*k - begin; L++){
                sum += m_pOrtho->H1(2*k - L)*Read(L, sz*2, 0);
            }
            m_pBuf1[n] = sum;
            n++;
        }

        for(int j = 0; j < sz*2; j++){
            m_pBuf0[j] = m_pBuf1[j];
        }
    }
    for(int i = 0; i < m_nN; i++){
        output[i] = m_pBuf0[i];
    }

    return true;
}

//多重解像度解析の逆変換
bool COrthoMra00::IMra(double* output, double* input, int level, int mul, int mode){
    double sum;

    m_nN = (0x01 << level)*mul;
    m_nLevel = level;

    if(level > MAXLEVEL || m_nN > m_nMaxN || m_pBuf == NULL){
        return false;
    }

    int begin1 = m_pOrtho->G0Begin();
    int end1 = m_pOrtho->G0End();
    int begin2 = m_pOrtho->G1Begin();
    int end2 = m_pOrtho->G1End();

    //分解アルゴリズムのループ
    for(int i = level; i > 0; i--){
        int n = 0;
        int sz = m_nN >> i;

        for(int k = 0; k < sz*2; k++){
            sum = 0;
            for(int L = k/2 - end1/2 - 5; L <= k/2 - begin1/2 + 5; L++){
                sum += m_pOrtho->G0(k - 2*L)*Read(L, sz, 0);
            }
            for(int L = k/2 - end2/2 - 5; L <= k/2 - begin2/2 + 5; L++){
                sum += m_pOrtho->G1(k - 2*L)*Read(L, sz, sz);
            }
            m_pBuf1[n] = sum;
            n++;
        }

        for(int j = 0; j < sz*2; j++){
            m_pBuf0[j] = m_pBuf1[j];
        }
    }

    if(mode){
        m_pInter->IInter(output, m_pBuf0, m_nN, m_pScal, m_nScalN, m_nScalOffset, 1);
    }
    else{
        for(int i = 0; i < m_nN; i++){
            output[i] = m_pBuf0[i];
        }
    }

    return true;
}

//ループ状にデータを読み込む
double COrthoMra00::Read(int n, int lpsz, int offset){
    while(n < 0){
        n += lpsz;
    }
    while(n >= lpsz){
        n -= lpsz;
    }

    return m_pBuf0[n + offset];
}


#ifdef ORTHO_MRA_TEST

int main(int argc, char* argv[]){
    int n;

    scanf("%d", &n);

    double input[n];
    for (int i = 0; i < n; i++){
        scanf("%lf", &input[i]);
    }

    COrthoMra00* cmra = new COrthoMra00();
    //現状、タイプは11しか選べない
    if(!cmra->Prepare(11, n)){
        printf("prepare failed\n");
        return 0;
    }

    double output[n];
    for(int i = 0; i < n; i++){
        output[i] = 0;
    }

    int level = 8;
    cmra->Mra(output, input, level, (int)n/pow(2, level), 1);
    printf("%d %d\n", n, level);
    //return 0;
    for(int i = 1; i <= level; i++){
        int num = (int)n/pow(2, i);
        int d = (int)pow(2, i);
        for(int j = 0; j < num; j++){
            printf("%d %d %.8f\n", i, j*d, output[j + num]);
        }
    }
    /*
    for(int i = 0; i < 512; i++){
        cgabor->Prepare(input, n);
        cgabor->SetParameter(2.0, _PI, (double)n/(i + 1));

        for(int j = 0; j < n; j++){
            printf("%.5f ", cgabor->Gabor(j));
        }
        printf("\n");
    }*/
}

#endif