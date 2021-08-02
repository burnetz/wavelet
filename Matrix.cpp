#include <math.h>
#include "Matrix.h"
#include <stdio.h>

CMatrix::CMatrix(){
    m_pBuf = new double[CMATRIX_MAX_N*CMATRIX_MAX_N*3];

    if(m_pBuf == NULL){
        return;
    }

    for(int i = 0; i < CMATRIX_MAX_N; i++){
        m_nInpMatrix[i] = &m_pBuf[CMATRIX_MAX_N*i];
    }

    for(int i = 0; i < CMATRIX_MAX_N; i++){
        m_nOutMatrix[i] = &m_pBuf[CMATRIX_MAX_N*i + CMATRIX_MAX_N*CMATRIX_MAX_N];
    }

    for(int i = 0; i < CMATRIX_MAX_N; i++){
        m_nBufMatrix[i] = &m_pBuf[CMATRIX_MAX_N*i + CMATRIX_MAX_N*CMATRIX_MAX_N*2];
    }

    m_nN = 3;
}

CMatrix::~CMatrix(){
    if(m_pBuf != NULL){
        delete m_pBuf;
    }
}

//バッファの存在を確認
bool CMatrix::InitCk(){
    return m_pBuf != NULL;
}

//行列式
double CMatrix::Determinant(int n){
    m_nN = n;

    if(n > CMATRIX_MAX_N || n < 1){
        return 0.0;
    }

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            m_nBufMatrix[i][j] = m_nInpMatrix[i][j];
        }
    }

    return DeterminantSub(n, m_nBufMatrix);
}

//行列式のサブルーチン
double CMatrix::DeterminantSub(int n, double* matrix[CMATRIX_MAX_N]){
    int i, j;
    double temp;

    if(n > CMATRIX_MAX_N || n < 1){
        return 0.0;
    }

    while(n > 2){
        if(matrix[0][0] == 0.0){
            for(i = 1; i < n; i++){
                if(matrix[0][i] != 0){
                    for(j = 0; j < n; j++){
                        temp = matrix[j][i];
                        matrix[j][i] = matrix[j][0];
                        matrix[j][0] = -temp;
                    }
                    break;
                }
            }
            if(i == n){
                return 0.0;
            }
        }

        for(i = 1; i < n; i++){
            temp = matrix[0][i] / matrix[0][0];
            for(j = 0; j < n; j++){
                matrix[j][i] -= temp*matrix[j][0];
            }
        }

        temp = matrix[0][0];
        for(i = 0; i < n - 1; i++){
            for(j = 0; j < n - 1; j++){
                matrix[i][j] = matrix[i + 1][j + 1];
            }
        }

        for(i = 0; i < n - 1; i++){
            matrix[i][0] *= temp;
        }

        n--;
    }

    if(n == 2){
        return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
    }
    else if(n == 1){
        return matrix[0][0];
    }
    else{
        return 0.0;
    }
}

//転置行列
bool CMatrix::TransposedMatrix(int n){
    m_nN = n;

    if(n > CMATRIX_MAX_N || n < 1){
        return false;
    }

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            m_nOutMatrix[i][j] = m_nInpMatrix[i][j];
        }
    }

    return TransposedMatrixSub(n, m_nOutMatrix);
}

//転置行列のサブルーチン
bool CMatrix::TransposedMatrixSub(int n, double* matrix[CMATRIX_MAX_N]){
    if(n > CMATRIX_MAX_N || n < 1){
        return false;
    }

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i > j){
                double temp = matrix[i][j];
                matrix[i][j] = matrix[j][i];
                matrix[j][i] = temp;
            }
        }
    }

    return true;
}

//逆行列

bool CMatrix::InverseMatrix(int n){
    int i, j, p, q;
    double temp;

    m_nN = n;

    if(n > CMATRIX_MAX_N || n < 2){
        return false;
    }

    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            int p_count = 0;
            for(p = 0; p < n; p++){
                if(p != i){
                    int q_count = 0;
                    for(q = 0; q < n; q++){
                        if(q != j){
                            m_nBufMatrix[p_count][q_count] = m_nInpMatrix[p][q];
                            q_count++;
                        }
                    }
                    p_count++;
                }
            }
            temp = DeterminantSub(n - 1, m_nBufMatrix);
            if((i + j) % 2 == 1){ //本とは書き方が違うがおそらく同じ意味
                temp *= -1;
            }
            m_nOutMatrix[j][i] = temp;
        }
    }

    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            m_nBufMatrix[i][j] = m_nInpMatrix[i][j];
        }
    }

    temp = DeterminantSub(n, m_nBufMatrix);
    if(temp == 0.0){
        return false;
    }

    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            m_nOutMatrix[i][j] = m_nOutMatrix[i][j] / temp;
        }
    }

    return true;
}



#ifdef MATRIX_TEST

int main(){
    CMatrix* mat = new CMatrix();

    int n = 3;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            mat->m_nInpMatrix[i][j] = i + j + 1;
        }
    }
    mat->m_nInpMatrix[2][1] = 0;

    printf("output\n");
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%.3f ",mat->m_nInpMatrix[i][j]);
        }
        printf("\n");
    }

    printf("det\n");
    printf("%f\n", mat->Determinant(n));

    printf("transpose\n");
    mat->TransposedMatrix(n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%.3f ",mat->m_nOutMatrix[i][j]);
        }
        printf("\n");
    }

    printf("inverse\n");
    mat->InverseMatrix(n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%.3f ",mat->m_nOutMatrix[i][j]);
        }
        printf("\n");
    }
}

#endif