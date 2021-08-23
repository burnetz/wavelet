#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "GraphOrthoMra00.h"
#include "OrthoMra00.h"

#define MAXLEVEL 8

CGraphOrthoMra00::CGraphOrthoMra00(){
    m_pMra = new COrthoMra00;

    for(int i = 0; i < 3; i++){
        m_pBuf[i] = NULL;
    }
    m_pWork = NULL;
}

CGraphOrthoMra00::~CGraphOrthoMra00(){
    if(m_pMra != NULL){
        delete m_pMra;
    }

    for(int i = 0; i < 3; i++){
        if(m_pBuf[i] != NULL){
            delete m_pBuf[i];
        }
    }

    if(m_pWork != NULL){
        delete m_pWork;
    }
}

bool CGraphOrthoMra00::Prepare(int wName, int maxX, int maxY){
    m_nWName = wName;
    m_nMaxX = maxX;
    m_nMaxY = maxY;

    m_nMaxLen = m_nMaxX;
    if(m_nMaxY > m_nMaxX){
        m_nMaxLen = m_nMaxY;
    }

    m_nMaxN = m_nMaxX * m_nMaxY;

    m_pMra->Prepare(m_nWName, m_nMaxLen);

    for(int i = 0; i < 3; i++){
        m_pBuf[i] = new double[m_nMaxN + 256];
        if(m_pBuf[i] == NULL){
            return false;
        }
    }

    m_pWork = new double[m_nMaxLen * 2 + 256];
    if(m_pWork == NULL){
        return false;
    }
    return true;
}

//2次元多重解像度解析
bool CGraphOrthoMra00::Mra(double* output, double* input, int nX, int nY, int level){

    if(nX > m_nMaxX || nY > m_nMaxY || level < 1 || level > MAXLEVEL ){
        return false;
    }

    double* work0;
    double* work1;

    work0 = m_pWork;
    work1 = work0 + m_nMaxLen;

    m_nX = nX;
    m_nY = nY;
    m_nLevel = level;

    for(int i = 0; i < m_nX*m_nY ; i++){
        m_pBuf[0][i] = input[i];
    }

    int xsz = m_nX;
    int ysz = m_nY;

    //レベルのループ
    for(int L = 1; L <= m_nLevel; L++){
        xsz = m_nX >> (L - 1);
        ysz = m_nY >> (L - 1);

        //縦のループ
        for(int y = 0; y < ysz ; y++){
            //横方向の処理
            for(int x = 0; x < xsz; x++){
                work0[x] = m_pBuf[0][x + y*m_nX];
            }
            m_pMra->Mra(work1, work0, 1, xsz >> 1, 0);
            for(int x = 0; x < xsz; x++){
                m_pBuf[1][x + y*m_nX] = work1[x];
            }
        }
        //横のループ
        for(int x = 0; x < xsz ; x++){
            //縦方向の処理
            for(int y = 0; y < ysz; y++){
                work0[y] = m_pBuf[1][x + y*m_nX];
            }
            m_pMra->Mra(work1, work0, 1, ysz >> 1, 0);
            for(int y = 0; y < ysz; y++){
                m_pBuf[2][x + y*m_nX] = work1[y];
            }
        }
        //低域成分LLの移動
        for(int y = 0; y < ysz; y++){
            for(int x = 0; x < xsz; x++){
                m_pBuf[0][x + y*m_nX] = m_pBuf[2][x + y*m_nX];
            }
        }
    }

    for(int y = 0; y < m_nY; y++){
        for(int x = 0; x < m_nX; x++){
            output[x + y*m_nX] = m_pBuf[2][x + y*m_nX];
        }
    }

    return true;
}

bool CGraphOrthoMra00::IMra(double* output, double* input, int nX, int nY, int level, double delta = 0){
    if(nX > m_nMaxX || nY > m_nMaxY || level < 1 || level > MAXLEVEL ){
        return false;
    }

    double* work0;
    double* work1;

    work0 = m_pWork;
    work1 = work0 + m_nMaxLen;

    m_nX = nX;
    m_nY = nY;
    m_nLevel = level;

    for(int y = 0; y < m_nY; y++){
        for(int x = 0; x < m_nX; x++){
            m_pBuf[0][x + y*m_nX] = input[x + y*m_nX];
        }
    }

    //レベルのループ
    for(int L = m_nLevel; L >= 1; L--){
        int xsz = m_nX >> (L - 1);
        int ysz = m_nY >> (L - 1);

        //横のループ
        for(int x = 0; x < xsz; x++){
            //縦方向の処理
            for(int y = 0; y < ysz; y++){
                work0[y] = m_pBuf[0][x + y*m_nX];
            }
            m_pMra->IMra(work1, work0, 1, ysz >> 1, 0);

            for(int y = 0; y < ysz; y++){
                m_pBuf[1][x + y*m_nX] = work1[y];
            }
        }

        //縦のループ
        for(int y = 0; y < ysz; y++){
            //横方向の処理
            for(int x = 0; x < xsz; x++){
                work0[x] = m_pBuf[1][x + y*m_nX];
            }
            m_pMra->IMra(work1, work0, 1, xsz >> 1, 0);

            for(int x = 0; x < xsz; x++){
                m_pBuf[2][x + y*m_nX] = work1[x];
            }
        }

        //低域成分LLの移動
        for(int y = 0; y < ysz; y++){
            for(int x = 0; x < xsz ;x++){
                m_pBuf[0][x + y*m_nX] = m_pBuf[2][x + y*m_nX];
            }
        }
    }

    double lambda = 3*delta;
    for(int i = 0; i < nX*nY ; i++){
        double tmp = m_pBuf[2][i];

        if(tmp < -lambda){
            tmp += lambda;
        }
        else if(tmp > lambda){
            tmp -= lambda;
        }
        else{
            tmp = 0;
        }
        output[i] = tmp;
    }

    return true;
    
}

#ifdef GRAPH_ORTHO_TEST

int main(int argc, char* argv[]){
    int n;

    //入力形式はPGMであること
    //GIMPで生成したPGMなら恐らく基本的に大丈夫
    char strBuf[128];
    fgets(strBuf, sizeof(strBuf), stdin);
    if(strBuf[0] != 'P' || strBuf[1] != '2'){
        printf("Invalid File Format!\n");
        exit(0);
    }

    int w = 0, h = 0;
    while(true){
        fgets(strBuf, sizeof(strBuf), stdin);
        if(strBuf[0] == '#'){
            continue;
        }
        sscanf(strBuf, "%d %d", &w, &h);
        break;
    }

    if(w == 0 || h == 0){
        printf("Invalid Value (width or height)\n");
        exit(0);
    }
    int whitePoint;
    scanf("%d", &whitePoint);

    double* input = new double[w*h];
    for (int i = 0; i < w*h; i++){
        scanf("%lf", &input[i]);
        input[i] /= 255.0;
    }

    CGraphOrthoMra00 *cgortho = new CGraphOrthoMra00;
    //現状、タイプは11しか選べない
    if(!cgortho->Prepare(11, w, h)){
        printf("prepare failed\n");
        return 0;
    }

    double* output = new double[w*h];
    for(int i = 0; i < w*h; i++){
        output[i] = 0;
    }

    int level = 3;
    if(cgortho->Mra(output, input, w, h, level)){
        printf("P2\n");
        printf("%d %d\n", w, h);
        printf("255\n");

        double max = -DBL_MAX;
        double min = DBL_MAX;
        for(int i = 0; i < w*h; i++){
            max = max < output[i] ? output[i] : max;
            min = min > output[i] ? output[i] : min;
        }

        if(min < 0){
            for(int i = 0; i < w*h; i++){
                output[i] -= min;
            }
            max -= min;
        }
        for(int i = 0; i < w*h; i++){
            int tmp = output[i]/max*255;
            if(tmp > 255 || tmp < 0){
                printf("error (%d)\n", tmp);
                return 0;
            }
            printf("%d\n", (int)(output[i]/max*255));
        }
    }
}

#endif

#ifdef GRAPH_ORTHO_I_TEST

int main(int argc, char* argv[]){
    int n;

    //入力形式はPGMであること
    //GIMPで生成したPGMなら恐らく基本的に大丈夫
    char strBuf[128];
    fgets(strBuf, sizeof(strBuf), stdin);
    if(strBuf[0] != 'P' || strBuf[1] != '2'){
        printf("Invalid File Format!\n");
        exit(0);
    }

    int w = 0, h = 0;
    while(true){
        fgets(strBuf, sizeof(strBuf), stdin);
        if(strBuf[0] == '#'){
            continue;
        }
        sscanf(strBuf, "%d %d", &w, &h);
        break;
    }

    if(w == 0 || h == 0){
        printf("Invalid Value (width or height)\n");
        exit(0);
    }
    int whitePoint;
    scanf("%d", &whitePoint);

    double* input = new double[w*h];
    for (int i = 0; i < w*h; i++){
        scanf("%lf", &input[i]);
        input[i];
    }

    CGraphOrthoMra00 *cgortho = new CGraphOrthoMra00;
    //現状、タイプは11しか選べない
    if(!cgortho->Prepare(11, w, h)){
        printf("prepare failed\n");
        return 0;
    }

    double* output = new double[w*h];
    for(int i = 0; i < w*h; i++){
        output[i] = 0;
    }

    int level = 3;
    if(cgortho->Mra(output, input, w, h, level)){
        if(cgortho->IMra(input, output, w, h, level)){
            printf("P2\n");
            printf("%d %d\n", w, h);
            printf("255\n");

            for(int i = 0; i < w*h; i++){
                int tmp = input[i];
                /*
                if(tmp > 255 || tmp < 0){
                    printf("error (%d)\n", tmp);
                    return 0;
                }*/
                printf("%d\n", tmp);
            }
        }
    }
}

#endif

#ifdef GRAPH_ORTHO_NOISE_TEST

int main(int argc, char* argv[]){
    if(argc != 2){
        printf("please set the value \"delta\"\n");
        return 0;
    }

    int n;

    //入力形式はPGMであること
    //GIMPで生成したPGMなら恐らく基本的に大丈夫
    char strBuf[128];
    fgets(strBuf, sizeof(strBuf), stdin);
    if(strBuf[0] != 'P' || strBuf[1] != '2'){
        printf("Invalid File Format!\n");
        exit(0);
    }

    int w = 0, h = 0;
    while(true){
        fgets(strBuf, sizeof(strBuf), stdin);
        if(strBuf[0] == '#'){
            continue;
        }
        sscanf(strBuf, "%d %d", &w, &h);
        break;
    }

    if(w == 0 || h == 0){
        printf("Invalid Value (width or height)\n");
        exit(0);
    }
    int whitePoint;
    scanf("%d", &whitePoint);

    double* input = new double[w*h];
    for (int i = 0; i < w*h; i++){
        scanf("%lf", &input[i]);
        input[i];
    }

    CGraphOrthoMra00 *cgortho = new CGraphOrthoMra00;
    //現状、タイプは11しか選べない
    if(!cgortho->Prepare(11, w, h)){
        printf("prepare failed\n");
        return 0;
    }

    double* output = new double[w*h];
    for(int i = 0; i < w*h; i++){
        output[i] = 0;
    }

    int level = 3;
    double delta = atof(argv[1]);
    if(delta < 0){
        printf("Invalid value \"delta\"\n");
        return 0;
    }
    if(cgortho->Mra(output, input, w, h, level)){
        if(cgortho->IMra(input, output, w, h, level, delta)){
            printf("P2\n");
            printf("%d %d\n", w, h);
            printf("255\n");

            for(int i = 0; i < w*h; i++){
                int tmp = input[i];

                printf("%d\n", tmp);
            }
        }
    }
}

#endif

#ifdef GRAPH_ORTHO_SHARP_TEST

int main(int argc, char* argv[]){

    int n;

    //入力形式はPGMであること
    //GIMPで生成したPGMなら恐らく基本的に大丈夫
    char strBuf[128];
    fgets(strBuf, sizeof(strBuf), stdin);
    if(strBuf[0] != 'P' || strBuf[1] != '2'){
        printf("Invalid File Format!\n");
        exit(0);
    }

    int w = 0, h = 0;
    while(true){
        fgets(strBuf, sizeof(strBuf), stdin);
        if(strBuf[0] == '#'){
            continue;
        }
        sscanf(strBuf, "%d %d", &w, &h);
        break;
    }

    if(w == 0 || h == 0){
        printf("Invalid Value (width or height)\n");
        exit(0);
    }
    int whitePoint;
    scanf("%d", &whitePoint);

    double* input = new double[w*h];
    for (int i = 0; i < w*h; i++){
        scanf("%lf", &input[i]);
        input[i];
    }

    CGraphOrthoMra00 *cgortho = new CGraphOrthoMra00;
    //現状、タイプは11しか選べない
    if(!cgortho->Prepare(11, w, h)){
        printf("prepare failed\n");
        return 0;
    }

    double* output = new double[w*h];
    for(int i = 0; i < w*h; i++){
        output[i] = 0;
    }

    int level = 2;

    if(cgortho->Mra(output, input, w, h, level)){
        //level2 LL以外の全てのウェーブレット係数を2倍にする
        for(int i = 0; i < w*h; i++){
            int x = i%w;
            int y = i/w;
            if(x >= w/4 || y >= h/4){
                //printf("%d %d\n", x, y);
                output[i] *= 2.0;
            }
        }
        
        if(cgortho->IMra(input, output, w, h, level)){
            printf("P2\n");
            printf("%d %d\n", w, h);
            printf("255\n");

            for(int i = 0; i < w*h; i++){
                int tmp = input[i];
                if(tmp < 0){
                    tmp = 0;
                }
                if(tmp > 255){
                    tmp = 255;
                }

                printf("%d\n", tmp);
            }
        }
    }
}

#endif

#ifdef GRAPH_ORTHO_CONTOUR_TEST

int main(int argc, char* argv[]){

    int n;

    //入力形式はPGMであること
    //GIMPで生成したPGMなら恐らく基本的に大丈夫
    char strBuf[128];
    fgets(strBuf, sizeof(strBuf), stdin);
    if(strBuf[0] != 'P' || strBuf[1] != '2'){
        printf("Invalid File Format!\n");
        exit(0);
    }

    int w = 0, h = 0;
    while(true){
        fgets(strBuf, sizeof(strBuf), stdin);
        if(strBuf[0] == '#'){
            continue;
        }
        sscanf(strBuf, "%d %d", &w, &h);
        break;
    }

    if(w == 0 || h == 0){
        printf("Invalid Value (width or height)\n");
        exit(0);
    }
    int whitePoint;
    scanf("%d", &whitePoint);

    double* input = new double[w*h];
    for (int i = 0; i < w*h; i++){
        scanf("%lf", &input[i]);
        input[i];
    }

    CGraphOrthoMra00 *cgortho = new CGraphOrthoMra00;
    //現状、タイプは11しか選べない
    if(!cgortho->Prepare(11, w, h)){
        printf("prepare failed\n");
        return 0;
    }

    double* output = new double[w*h];
    for(int i = 0; i < w*h; i++){
        output[i] = 0;
    }

    int level = 2;

    if(cgortho->Mra(output, input, w, h, level)){
        //level2 LL以外の全てのウェーブレット係数を2倍にする
        //level2 LLは灰色で置き換える
        for(int i = 0; i < w*h; i++){
            int x = i%w;
            int y = i/w;
            if(x >= w/4 || y >= h/4){
                //printf("%d %d\n", x, y);
                output[i] *= 2.0;
            }
            else {
                output[i] = 128;
            }
        }
        
        if(cgortho->IMra(input, output, w, h, level)){
            printf("P2\n");
            printf("%d %d\n", w, h);
            printf("255\n");

            for(int i = 0; i < w*h; i++){
                int tmp = input[i];
                if(tmp < 0){
                    tmp = 0;
                }
                if(tmp > 255){
                    tmp = 255;
                }

                printf("%d\n", tmp);
            }
        }
    }
}

#endif