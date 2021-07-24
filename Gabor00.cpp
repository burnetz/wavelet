#include <stdio.h>
#include <math.h>
#include "Gabor00.h"

#define _PI 3.14159265358979323846264338327950288419716939937510
#define WIN_RATE 0.01

CGabor00::CGabor00()
{
    for (int i = 0; i < 8; i++)
    {
        m_pWork[i] = NULL;
    }
}

CGabor00::~CGabor00()
{
    for (int i = 0; i < 8; i++)
    {
        if (m_pWork[i] != NULL)
        {
            delete m_pWork[i];
        }
    }
}

bool CGabor00::Prepare(double *input, int inputN)
{
    m_pInput = input;
    m_nInputN = inputN;

    return true;
}

bool CGabor00::SetParameter(double sigma, double omega, double a){
    bool ret = false;

    for(int i = 0; i < 2; i++){
        if(m_pWork[i] != NULL){
            delete m_pWork[i];
        }
        m_pWork[i] = NULL;
    }
    if(sigma <= 0.0 || omega <= 0.0 || a <= 0.0){
        return ret;
    }

    m_nGaborOffset = (int)(sqrt(-2.0*sigma*sigma*a*a*log(WIN_RATE)));
    m_nGaborN = m_nGaborOffset*2 + 1;

    int i;
    for(i = 0; i < 2; i++){
        m_pWork[i] = new double[m_nGaborN];
        if(m_pWork[i] == NULL){
            break;
        }
    }

    if(i == 2){
        m_pGaborR = m_pWork[0];
        m_pGaborI = m_pWork[1];
        ret = true;
    }

    if(ret == true){
        for(i = 0 ; i < m_nGaborN; i++){
            double x = (double)(i - m_nGaborOffset)/a;
            double d = exp(-x*x/(2.0*sigma*sigma))/sqrt(2.0*_PI*sigma*sigma);
            m_pGaborR[i] = d*cos(omega*x)/sqrt(a);
            m_pGaborI[i] = d*sin(omega*x)/sqrt(a);
        }
    }

    return ret;

}

double CGabor00::Gabor(int b){
    if(m_pWork[0] == NULL){
        return 0.0;
    }

    double dr = 0.0;
    double di = 0.0;

    double real2, imag2;
    for(int i = -(m_nGaborN - 1)/2; i<= (m_nGaborN - 1)/2; i++){
        double real0 = ReadInput(b + i);
        double imag0 = 0.0;

        double real1 = ReadGaborR(i);
        double imag1 = ReadGaborI(i);

        ComplexMultiply(real2, imag2, real0, imag0, real1, -imag1);

        dr += real2;
        di += imag2;
    }
    return dr*dr + di*di;
}

double CGabor00::ReadInput(int x){
    if(x < 0 || x >= m_nInputN){
        return 0.0;
    }
    return m_pInput[x];
}

double CGabor00::ReadGaborR(int x){
    x += m_nGaborOffset;

    if(x < 0 || x >= m_nGaborN){
        return 0.0;
    }
    return m_pGaborR[x];
}

double CGabor00::ReadGaborI(int x){
    x += m_nGaborOffset;

    if(x < 0 || x >= m_nGaborN){
        return 0.0;
    }
    return m_pGaborI[x];
}

void CGabor00::ComplexMultiply(double &rOutput, double &iOutput, double rInput0, double iInput0, double rInput1, double iInput1){
    rOutput = rInput0*rInput1 - iInput0*iInput1;
    iOutput = iInput0*rInput1 + rInput0*iInput1;
}

#ifdef GABOR_TEST

int main(int argc, char* argv[]){
    int n;

    scanf("%d", &n);

    double input[n];
    for (int i = 0; i < n; i++){
        scanf("%lf", &input[i]);
    }

    CGabor00* cgabor = new CGabor00;

    printf("%d 512\n", n);
    for(int i = 0; i < 512; i++){
        cgabor->Prepare(input, n);
        cgabor->SetParameter(2.0, _PI, (double)n/(i + 1));

        for(int j = 0; j < n; j++){
            printf("%.5f ", cgabor->Gabor(j));
        }
        printf("\n");
    }
}

#endif