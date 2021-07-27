#define CMATRIX_MAX_N 64

class CMatrix{
    public:
    CMatrix();
    ~CMatrix();

    double* m_nInpMatrix[CMATRIX_MAX_N];
    double* m_nOutMatrix[CMATRIX_MAX_N];
    double* m_nBufMatrix[CMATRIX_MAX_N];

    //行列の1辺の最大値
    int MaxN(){
        return CMATRIX_MAX_N;
    }

    bool InitCk();

    double Determinant(int n);
    bool TransposedMatrix(int n);
    bool InverseMatrix(int n);

    protected:
    double* m_pBuf;
    int m_nN;

    double DeterminantSub(int n, double* matrix[CMATRIX_MAX_N]);
    bool TransposedMatrixSub(int n, double* matrix[CMATRIX_MAX_N]);
};