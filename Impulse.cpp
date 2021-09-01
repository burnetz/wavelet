#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "RISplineMra00.h"

int main(){
	const int n = 1024;

	double input[n];

	for(int i = 0; i < n; i++){
		if((64 + i)%129 == 0){
			input[i] = 1.0;
		}
		else{
			input[i] = 0.0;
		}
	}

	printf("2 %d\n", n);
    for(int i = 0; i < n; i++){
        printf("%d %f\n", i, input[i]);
    }

	CRISplineMra00* cri = new CRISplineMra00;

	int level = 3;
    cri->Prepare(n);

	double** output;
    output = (double**)malloc(sizeof(double*) * 2);

	for(int i = 0; i < 2; i++){
        output[i] = (double*)malloc(sizeof(double) * n);
		for(int j = 0; j < n; j++){
            output[i][j] = 0;
        }
	}
    int mul = n / pow(2, level);
    cri->Mra((double**)output, input, level, mul, 0);
    //レベル3だけ残す
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < n; j++){
            if(j < n/8 || j >= n/4){
                output[i][j] = 0;
            }
        }
    }
    cri->IMra(input, (double **)output, level, mul, 0);


    for(int i = 0; i < n; i++){
        printf("%d %f\n", i, input[i]);
    }

    return 0;
}