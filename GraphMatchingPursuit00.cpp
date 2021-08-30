#include <math.h>
#include <stdio.h>
#include <float.h>
#include "GraphMatchingPursuit00.h"
#include "MatchingPursuit00.h"

#define MAXLEVEL 8

CGraphMatchingPursuit00::CGraphMatchingPursuit00()
{
	m_pBuf = NULL;
}

CGraphMatchingPursuit00::~CGraphMatchingPursuit00()
{
	if( m_pBuf != NULL )
		delete m_pBuf;
}

// 2次元マッチング追跡
// outputにはinputの(level + 1)^2倍の領域を用意すること
bool CGraphMatchingPursuit00::Pursuit( double* output, double* input, int nX, int nY,
int wName, int level, double alpha, double epsilon, double rho )
{
	int i, x, y;
	double* outBuf[ MAXLEVEL + 1 ];

	if( level < 1 || MAXLEVEL < level )
		return false;

	int max = nX;
	if( max < nY )max = nY;

	if( m_pBuf != NULL )
		delete m_pBuf;
	m_pBuf = new double[ (nX + nY)*(level + 2)*2 + nX*nY + 256 ];
	if( m_pBuf == NULL )
		return false;
	double* inp = m_pBuf;
	double* out = inp + (nX + nY)*(level + 2);
	double* box = out + (nX + nY)*(level + 2);

	CMatchingPursuit00 pursuit;

	for( y=0; y<nY*(level + 1); y++ )	// 出力バッファのクリア
		for( x=0; x<nX*(level + 1); x++ )
			output[ x + y*nX*(level + 1) ] = 0.0;

	if( pursuit.Prepare( wName, max, level ) == false )
		m_pBuf = NULL;

	// 最初のLの字の横
	for( y=0; y<nY; y++ )
	{
		// 一時的な出力バッファの配置
		for( i=0; i<=level+1; i++ )
		{
			outBuf[i] = out + (level+1-i)*nX;
		}

		// マッチング追跡を実行
		pursuit.Pursuit( outBuf, input + y*nX, nX, level, alpha, epsilon, rho, 0 );

		for( x=0; x<nX*(level+1); x++ )	// 出力
		{
			output[ x + y*nX*(level+1) ] = out[x];
		}
		for( x=0; x<nX; x++ )	// 残留成分と低域成分を加算
		{
			output[ x + y*nX*(level+1) ] += out[x + nX*(level+1)];
		}
	}

	// 最初のLの字の縦
	for( x=0; x<nX; x++ )
	{
		for( y=0; y<nY; y++ )
		{
			inp[y] = output[ x + y*nX*(level + 1) ];
		}
		for( i=1; i<=level+1; i++ )
		{
			outBuf[i] = out + (level+1-i)*nY;
		}

		pursuit.Pursuit( outBuf, inp, nY, level, alpha, epsilon, rho, 0 );

		for( y=0; y<nY*(level+1); y++ )
		{
			output[ x + y*nX*(level+1) ] = out[y];
		}
		for( y=0; y<nY; y++ )
		{
			output[ x + y*nX*(level+1) ] += out[y + nY*(level+1)];
		}
	}

	// 2回目のLの字の縦
	for( x=0; x<nX; x++ )
	{
		for( y=0; y<nY; y++ )
		{
			inp[y] = input[ x + y*nX ];
		}
		for( i=1; i<=level+1; i++ )
		{
			outBuf[i] = out + (level+1-i)*nY;
		}

		pursuit.Pursuit( outBuf, inp, nY, level, alpha, epsilon, rho, 0 );

		// 残留成分と低域成分を加算してボックスに一時保存
		for( y=0; y<nY; y++ )
		{
			box[ x + y*nX ] = out[y] + out[y + nY*(level+1)];
		}
		for( y=nY; y<nY*(level + 1); y++ )
		{
			// 最初のLの字と平均する
			output[ x + y*nX*(level+1) ] += out[y];
			output[ x + y*nX*(level+1) ] /= 2.0;
		}
	}

	// 2回目のLの字の横
	for( y=0; y<nY; y++ )
	{
		for( i=1; i<=level+1; i++ )
		{
			outBuf[i] = out + (level+1-i)*nX;
		}

		pursuit.Pursuit( outBuf, box + y*nX, nX, level, alpha, epsilon, rho, 0 );

		for( x=0; x<nX*(level+1); x++ )
		{
			output[ x + y*nX*(level+1) ] += out[x];
		}
		for( x=0; x<nX; x++ )
		{
			output[ x + y*nX*(level+1) ] += out[x + nX*(level+1)];
		}
		for( x=0; x<nX*(level+1); x++ )
		{
			output[ x + y*nX*(level+1) ] /= 2.0;
		}
	}

	// 縦への引き延ばし
	for( x=nX; x<nX*(level + 1); x++ )
	{
		for( y=0; y<nY; y++ )
		{
			inp[y] = output[ x + y*nX*(level + 1) ];
		}
		for( i=1; i<=level+1; i++ )
		{
			outBuf[i] = out + (level+1-i)*nY;
		}

		pursuit.Pursuit( outBuf, inp, nY, level, alpha, epsilon, rho, 0 );

		for( y=0; y<nY*(level+1); y++ )
		{
			output[ x + y*nX*(level+1) ] = out[y];
		}
		for( y=0; y<nY; y++ )
		{
			output[ x + y*nX*(level+1) ] += out[y + nY*(level+1)];
		}
	}

	// 横への引き伸ばし
	for( y=nY; y<nY*(level + 1); y++ )
	{
		for( i=1; i<=level+1; i++ )
		{
			outBuf[i] = out + (level+1-i)*nX;
		}

		pursuit.Pursuit( outBuf, output + y*nX*(level+1), nX, level, alpha, epsilon, rho, 0 );

		for( x=0; x<nX; x++ )
		{
			output[ x + y*nX*(level+1) ] = out[x] + out[x + nX*(level+1)];
		}
		for( x=nX; x<nX*(level+1); x++ )
		{
			// 縦の結果への加算
			output[ x + y*nX*(level+1) ] += out[x];
		}
	}

	if( m_pBuf != NULL )
		delete m_pBuf;
	m_pBuf = NULL;

	return true;
}

// 2次元マッチング追跡の逆変換
bool CGraphMatchingPursuit00::IPursuit( double* output, double* input, int nX, int nY,
int wName, int level )
{
	int i, x, y;
	double* inputBuf[ MAXLEVEL + 1 ];

	if( level < 1 || level > MAXLEVEL )
		return false;

	int max = nX;
	if( max < nY )max = nY;

	if( m_pBuf != NULL )
		delete m_pBuf;
	m_pBuf = new double[ nX*nY*(level+1) + nY*(level+2) + nY + 256 ];
	if( m_pBuf == NULL )
		return false;
	double* box = m_pBuf;
	double* inp = box + nX*nY*(level+1);
	double* out = inp + nY*(level+2);

	CMatchingPursuit00 pursuit;
	
	if( pursuit.Prepare( wName, max, level ) == false )
		return false;

	// 縦の逆変換
	for( x=0; x<nX*(level+1); x++ )
	{
		// 一時的な入力バッファにデータを転送
		for( y=0; y<nY*(level+1); y++ )
		{
			inp[y] = input[ x + y*nX*(level+1) ];
		}
		// ダミーの領域を用意
		for( y=nY*(level+1); y<nY*(level+2); y++ )
		{
			inp[y] = 0;
		}
		// 一時的な入力バッファの配置
		for( i=0; i<=level+1; i++ )
		{
			inputBuf[i] = inp + (level+1-i)*nY;
		}

		// マッチング追跡の逆変換
		pursuit.IPursuit( out, inputBuf, nY, level, 0 );

		// ボックスに一時保存
		for( y=0; y<nY; y++ )
		{
			box[ x + y*nX*(level+1) ] = out[y];
		}
	}

	// 横の逆変換
	for( y=0; y<nY; y++ )
	{
		for( x=0; x<nX*(level+1); x++ )
		{
			inp[x] = box[x + y*nX*(level+1)];
		}
		for( x=nX*(level+1); x<nX*(level+2); x++ )
		{
			inp[x] = 0;
		}
		for( i=0; i<=level+1; i++ )
		{
			inputBuf[i] = inp + (level+1-i)*nX;
		}

		pursuit.IPursuit( output + y*nX, inputBuf, nX, level, 0 );
	}

	if( m_pBuf != NULL )
		delete m_pBuf;
	m_pBuf = NULL;

	return true;
}

#ifdef GRAPH_MATCHING_TEST

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

    CGraphMatchingPursuit00 *cgmatch = new CGraphMatchingPursuit00;
    
    int level = 2;
	int outputSize = w*h*(level + 1)*(level + 1);
    double* output = new double[outputSize];
    for(int i = 0; i < outputSize; i++){
        output[i] = 0;
    }

	//3次のO splineは22
    if(cgmatch->Pursuit(output, input, w, h, 22, level, 0.99999, 0.00001, 1)){
        printf("P2\n");
        printf("%d %d\n", w*(level + 1), h*(level + 1));
        printf("255\n");

        double max = -DBL_MAX;
        double min = DBL_MAX;
        for(int i = 0; i < outputSize; i++){
            max = max < output[i] ? output[i] : max;
            min = min > output[i] ? output[i] : min;
        }

        if(min < 0){
            for(int i = 0; i < outputSize; i++){
                output[i] -= min;
            }
            max -= min;
        }
        for(int i = 0; i < outputSize; i++){
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

#ifdef GRAPH_MATCHING_SHARP_TEST

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
    }

    CGraphMatchingPursuit00 *cgmatch = new CGraphMatchingPursuit00;
    
    int level = 2;
	int outputScale = (level + 1);
	int outputSize = w*h*outputScale*outputScale;
    double* output = new double[outputSize];
    for(int i = 0; i < outputSize; i++){
        output[i] = 0;
    }

	//3次のO splineは22
    if(cgmatch->Pursuit(output, input, w, h, 22, level, 0.99999, 0.00001, 1)){  
		//level2 LL以外の全てのウェーブレット係数を2倍にする
        for(int i = 0; i < outputSize; i++){
            int x = i%(w*outputScale);
            int y = i/(w*outputScale);
            if(x >= w || y >= h){
                //printf("%d %d\n", x, y);
                //output[i] *= 2.0;
            }
        }
        
        if(cgmatch->IPursuit(input, output, w, h, 22, level)){
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

#ifdef GRAPH_MATCHING_CONTOUR_TEST

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
    }

    CGraphMatchingPursuit00 *cgmatch = new CGraphMatchingPursuit00;
    
    int level = 2;
	int outputScale = (level + 1);
	int outputSize = w*h*outputScale*outputScale;
    double* output = new double[outputSize];
    for(int i = 0; i < outputSize; i++){
        output[i] = 0;
    }

	//3次のO splineは22
    if(cgmatch->Pursuit(output, input, w, h, 22, level, 0.99999, 0.00001, 1)){  
        //level2 LL以外の全てのウェーブレット係数を2倍にする
        //level2 LLは灰色で置き換える
        for(int i = 0; i < outputSize; i++){
            int x = i%(w*outputScale);
            int y = i/(w*outputScale);
            if(x >= w || y >= h){
                //printf("%d %d\n", x, y);
                output[i] *= 2.0;
            }
            else {
                output[i] = 128;
            }
        }
        
        if(cgmatch->IPursuit(input, output, w, h, 22, level)){
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