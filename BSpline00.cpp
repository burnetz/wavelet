#include <math.h>

#include "BSpline00.h"

CBSpline00::CBSpline00()
{
	m_pBuf = NULL;

	m_nH0Begin = 0;
	m_nH0End = 0;
	m_nH0Offset = 0;

	m_nH1Begin = 0;
	m_nH1End = 0;
	m_nH1Offset = 0;

	m_nG0Begin = 0;
	m_nG0End = 0;
	m_nG0Offset = 0;

	m_nG1Begin = 0;
	m_nG1End = 0;
	m_nG1Offset = 0;

	m_nSc1Begin = 0;
	m_nSc1End = 0;
	m_nSc1Offset = 0;

	m_nSc2Begin = 0;
	m_nSc2End = 0;
	m_nSc2Offset = 0;

	m_nBetaSc2Begin = 0;
	m_nBetaSc2End = 0;
	m_nBetaSc2Offset = 0;
}

CBSpline00::~CBSpline00()
{
	if( m_pBuf != NULL )
		delete m_pBuf;
}

double CBSpline00::H0( int n )	// 分解数列H0
{
	if( m_pBuf == NULL )
		return 0.0;

	if( n<m_nH0Begin || n>m_nH0End )
		return 0.0;
	else
		return m_pH0[n + m_nH0Offset];
}

double CBSpline00::H1( int n )	// 分解数列H1
{
	if( m_pBuf == NULL )
		return 0.0;

	if( n<m_nH1Begin || n>m_nH1End )
		return 0.0;
	else
		return m_pH1[n + m_nH1Offset];
}

double CBSpline00::G0( int n )	// 再構成数列G0
{
	if( m_pBuf == NULL )
		return 0.0;

	if( n<m_nG0Begin || n>m_nG0End )
		return 0.0;
	else
		return m_pG0[n + m_nG0Offset];
}

double CBSpline00::G1( int n )	// 再構成数列G1
{
	if( m_pBuf == NULL )
		return 0.0;

	if( n<m_nG1Begin || n>m_nG1End )
		return 0.0;
	else
		return m_pG1[n+ m_nG1Offset];
}

double CBSpline00::Sc1( int n )	// オリジナルのΦ(x)
{
	if( m_pBuf == NULL )
		return 0.0;

	if( n<m_nSc1Begin || n>m_nSc1End )
		return 0.0;
	else
		return m_pSc1[n + m_nSc1Offset];
}

double CBSpline00::Sc2( int n )	// 原点を頂点とするΦ(x)
{
	if( m_pBuf == NULL )
		return 0.0;

	if( n<m_nSc2Begin || n>m_nSc2End )
		return 0.0;
	else
		return m_pSc2[n + m_nSc2Offset];
}

double CBSpline00::BetaSc2( int n ) // 原点を頂点とするΦ(x)の補間
{
	if( m_pBuf == NULL )
		return 0.0;

	if( n<m_nBetaSc2Begin || n>m_nBetaSc2End )
		return 0.0;
	else
		return m_pBetaSc2[ n + m_nBetaSc2Offset];
}


bool CBSpline00::Prepare( int m )
{
	int i;

	// バッファの用意
	if( m_pBuf != NULL )
	{
		delete m_pBuf;
	}
	m_pBuf = new double[256*7 + 256];
	if( m_pBuf != NULL )
	{
		for( i=0; i<256*7; i++ )	// 万が一はみ出して読んでもそのデータが0になるように
			m_pBuf[i] = 0.0;

		m_pH0 = m_pBuf;
		m_pH1 = m_pH0 + 256;
		m_pG0 = m_pH1 + 256;
		m_pG1 = m_pG0 + 256;
		m_pSc1 = m_pG1 + 256;
		m_pSc2 = m_pSc1 + 256;
		m_pBetaSc2 = m_pSc2 + 256;
	}
	else
	{
		return false;
	}

	switch( m )
	{

	// 3階Splineウェーブレット
	case 3:
		// 分解数列H0
		m_nH0Offset = 128;
		m_nH0Begin = -36;
		m_nH0End = 35;
		m_pH0[m_nH0Offset +  0] =  0.926791199152;
		m_pH0[m_nH0Offset +  1] =  0.048053527291;
		m_pH0[m_nH0Offset +  2] = -0.344757589865;
		m_pH0[m_nH0Offset +  3] = -0.036679011431;
		m_pH0[m_nH0Offset +  4] =  0.146104169219;
		m_pH0[m_nH0Offset +  5] =  0.016482116502;
		m_pH0[m_nH0Offset +  6] = -0.062808010228;
		m_pH0[m_nH0Offset +  7] = -0.007126488039;
		m_pH0[m_nH0Offset +  8] =  0.027039234841;
		m_pH0[m_nH0Offset +  9] =  0.003069769822;
		m_pH0[m_nH0Offset + 10] = -0.011642240637;
		m_pH0[m_nH0Offset + 11] = -0.001321822399;
		m_pH0[m_nH0Offset + 12] =  0.005012853742;
		m_pH0[m_nH0Offset + 13] =  0.000569146607;
		m_pH0[m_nH0Offset + 14] = -0.002158410916;
		m_pH0[m_nH0Offset + 15] = -0.000245060814;
		m_pH0[m_nH0Offset + 16] =  0.000929358571;
		m_pH0[m_nH0Offset + 17] =  0.000105517646;
		m_pH0[m_nH0Offset + 18] = -0.000400159019;
		m_pH0[m_nH0Offset + 19] = -0.000045434449;
		m_pH0[m_nH0Offset + 20] =  0.000172298913;
		m_pH0[m_nH0Offset + 21] =  0.000019565630;
		m_pH0[m_nH0Offset + 22] = -0.000074188370;
		m_pH0[m_nH0Offset + 23] = -0.000008430694;
		m_pH0[m_nH0Offset + 24] =  0.000031945322;
		m_pH0[m_nH0Offset + 25] =  0.000003644485;
		m_pH0[m_nH0Offset + 26] = -0.000013758672;
		m_pH0[m_nH0Offset + 27] = -0.000001602759;
		m_pH0[m_nH0Offset + 28] =  0.000005932936;
		m_pH0[m_nH0Offset + 29] =  0.000000768143;
		m_pH0[m_nH0Offset + 30] = -0.000002573944;
		m_pH0[m_nH0Offset + 31] = -0.000000515531;
		m_pH0[m_nH0Offset + 32] =  0.000001128969;
		m_pH0[m_nH0Offset + 33] =  0.000000733632;
		m_pH0[m_nH0Offset + 34] =  0.000000147537;
		m_pH0[m_nH0Offset + 35] =  0.000000005087;
		for( i=0; i<36; i++ )
			m_pH0[m_nH0Offset - 1 - i] = m_pH0[m_nH0Offset + i];

		// 分解数列H1
		m_nH1Offset = 128;
		m_nH1Begin = -34;
		m_nH1End = 33;
		m_pH1[m_nH1Offset +  0] = -0.538837275001;
		m_pH1[m_nH1Offset +  1] = -0.083989158759;
		m_pH1[m_nH1Offset +  2] =  0.253673413647;
		m_pH1[m_nH1Offset +  3] =  0.029784510337;
		m_pH1[m_nH1Offset +  4] = -0.110159128836;
		m_pH1[m_nH1Offset +  5] = -0.012549558151;
		m_pH1[m_nH1Offset +  6] =  0.047472040290;
		m_pH1[m_nH1Offset +  7] =  0.005391682418;
		m_pH1[m_nH1Offset +  8] = -0.020442024210;
		m_pH1[m_nH1Offset +  9] = -0.002321014914;
		m_pH1[m_nH1Offset + 10] =  0.008801906395;
		m_pH1[m_nH1Offset + 11] =  0.000999349767;
		m_pH1[m_nH1Offset + 12] = -0.003789887114;
		m_pH1[m_nH1Offset + 13] = -0.000430294355;
		m_pH1[m_nH1Offset + 14] =  0.001631832080;
		m_pH1[m_nH1Offset + 15] =  0.000185273939;
		m_pH1[m_nH1Offset + 16] = -0.000702626628;
		m_pH1[m_nH1Offset + 17] = -0.000079774014;
		m_pH1[m_nH1Offset + 18] =  0.000302533605;
		m_pH1[m_nH1Offset + 19] =  0.000034347852;
		m_pH1[m_nH1Offset + 20] = -0.000130263282;
		m_pH1[m_nH1Offset + 21] = -0.000014787315;
		m_pH1[m_nH1Offset + 22] =  0.000056087625;
		m_pH1[m_nH1Offset + 23] =  0.000006362354;
		m_pH1[m_nH1Offset + 24] = -0.000024148709;
		m_pH1[m_nH1Offset + 25] = -0.000002728559;
		m_pH1[m_nH1Offset + 26] =  0.000010394960;
		m_pH1[m_nH1Offset + 27] =  0.000001149506;
		m_pH1[m_nH1Offset + 28] = -0.000004469110;
		m_pH1[m_nH1Offset + 29] = -0.000000436149;
		m_pH1[m_nH1Offset + 30] =  0.000001908302;
		m_pH1[m_nH1Offset + 31] =  0.000000052652;
		m_pH1[m_nH1Offset + 32] = -0.000000774847;
		m_pH1[m_nH1Offset + 33] =  0.000000258282;
		for( i=0; i<34; i++ )
			m_pH1[m_nH1Offset - 1 - i] = - m_pH1[m_nH1Offset + i];

		// 再構成数列G0
		m_nG0Offset = 128;
		m_nG0Begin = -1;
		m_nG0End = 2;
		m_pG0[m_nG0Offset -  1] =  0.176776695297;
		m_pG0[m_nG0Offset +  0] =  0.530330085890;
		m_pG0[m_nG0Offset +  1] =  0.530330085890;
		m_pG0[m_nG0Offset +  2] =  0.176776695297;

		// 再構成数列G1
		m_nG1Offset = 128;
		m_nG1Begin = -3;
		m_nG1End = 4;
		m_pG1[m_nG0Offset -  3] =  0.003482026460;
		m_pG1[m_nG0Offset -  2] = -0.100978767343;
		m_pG1[m_nG0Offset -  1] =  0.511857889634;
		m_pG1[m_nG0Offset +  0] = -1.055054017409;
		m_pG1[m_nG0Offset +  1] =  1.055054017409;
		m_pG1[m_nG0Offset +  2] = -0.511857889634;
		m_pG1[m_nG0Offset +  3] =  0.100978767343;
		m_pG1[m_nG0Offset +  4] = -0.003482026460;

		// オリジナルのΦ(x)
		m_nSc1Offset = 128;
		m_nSc1Begin = 0;
		m_nSc1End = 1;
		m_pSc1[m_nSc1Offset + 0] = 0.674199835179;
		m_pSc1[m_nSc1Offset + 1] = 0.674199835179;

		// 原点を頂点とするΦ(x)
		m_nSc2Offset = 128;
		m_nSc2Begin = -1;
		m_nSc2End = 1;
		m_pSc2[m_nSc2Offset - 1] = 0.168549958795;
		m_pSc2[m_nSc2Offset + 0] = 1.011299752769;
		m_pSc2[m_nSc2Offset + 1] = 0.168549958795;

		// 原点を頂点とするΦ(x)による補間
		m_nBetaSc2Offset = 128;
		m_nBetaSc2Begin = -15;
		m_nBetaSc2End = 15;
		m_pBetaSc2[m_nBetaSc2Offset +  0] =  1.048808848170;
		m_pBetaSc2[m_nBetaSc2Offset +  1] = -0.179947149672;
		m_pBetaSc2[m_nBetaSc2Offset +  2] =  0.030874049863;
		m_pBetaSc2[m_nBetaSc2Offset +  3] = -0.005297149506;
		m_pBetaSc2[m_nBetaSc2Offset +  4] =  0.000908847171;
		m_pBetaSc2[m_nBetaSc2Offset +  5] = -0.000155933522;
		m_pBetaSc2[m_nBetaSc2Offset +  6] =  0.000026753963;
		m_pBetaSc2[m_nBetaSc2Offset +  7] = -0.000004590254;
		m_pBetaSc2[m_nBetaSc2Offset +  8] =  0.000000787563;
		m_pBetaSc2[m_nBetaSc2Offset +  9] = -0.000000135124;
		m_pBetaSc2[m_nBetaSc2Offset + 10] =  0.000000023184;
		m_pBetaSc2[m_nBetaSc2Offset + 11] = -0.000000003978;
		m_pBetaSc2[m_nBetaSc2Offset + 12] =  0.000000000682;
		m_pBetaSc2[m_nBetaSc2Offset + 13] = -0.000000000117;
		m_pBetaSc2[m_nBetaSc2Offset + 14] =  0.000000000020;
		m_pBetaSc2[m_nBetaSc2Offset + 15] = -0.000000000003;
		for( i=1; i<16; i++ )
			m_pBetaSc2[m_nBetaSc2Offset - i] = m_pBetaSc2[m_nBetaSc2Offset + i];

		break;

	// 4階Splineウェーブレット
	case 4:
		// 分解数列H0
		m_nH0Offset = 128;
		m_nH0Begin = -37;
		m_nH0End = 37;
		m_pH0[m_nH0Offset +  0] =  1.263123023929;
		m_pH0[m_nH0Offset +  1] =  0.566648257953;
		m_pH0[m_nH0Offset +  2] = -0.399107854109;
		m_pH0[m_nH0Offset +  3] = -0.329405166632;
		m_pH0[m_nH0Offset +  4] =  0.182551735438;
		m_pH0[m_nH0Offset +  5] =  0.178837838534;
		m_pH0[m_nH0Offset +  6] = -0.093933246127;
		m_pH0[m_nH0Offset +  7] = -0.096030209541;
		m_pH0[m_nH0Offset +  8] =  0.049817225355;
		m_pH0[m_nH0Offset +  9] =  0.051440030340;
		m_pH0[m_nH0Offset + 10] = -0.026609389196;
		m_pH0[m_nH0Offset + 11] = -0.027539380933;
		m_pH0[m_nH0Offset + 12] =  0.014236512200;
		m_pH0[m_nH0Offset + 13] =  0.014741872072;
		m_pH0[m_nH0Offset + 14] = -0.007619648555;
		m_pH0[m_nH0Offset + 15] = -0.007891151730;
		m_pH0[m_nH0Offset + 16] =  0.004078497759;
		m_pH0[m_nH0Offset + 17] =  0.004224077333;
		m_pH0[m_nH0Offset + 18] = -0.002183040415;
		m_pH0[m_nH0Offset + 19] = -0.002261234915;
		m_pH0[m_nH0Offset + 20] =  0.001168375058;
		m_pH0[m_nH0Offset + 21] =  0.001210708700;
		m_pH0[m_nH0Offset + 22] = -0.000625105074;
		m_pH0[m_nH0Offset + 23] = -0.000648654847;
		m_pH0[m_nH0Offset + 24] =  0.000334039977;
		m_pH0[m_nH0Offset + 25] =  0.000348307047;
		m_pH0[m_nH0Offset + 26] = -0.000177743582;
		m_pH0[m_nH0Offset + 27] = -0.000188485192;
		m_pH0[m_nH0Offset + 28] =  0.000093134094;
		m_pH0[m_nH0Offset + 29] =  0.000104692951;
		m_pH0[m_nH0Offset + 30] = -0.000045892022;
		m_pH0[m_nH0Offset + 31] = -0.000063001374;
		m_pH0[m_nH0Offset + 32] =  0.000015451730;
		m_pH0[m_nH0Offset + 33] =  0.000045263764;
		m_pH0[m_nH0Offset + 34] =  0.000024132083;
		m_pH0[m_nH0Offset + 35] =  0.000005324143;
		m_pH0[m_nH0Offset + 36] =  0.000000394271;
		m_pH0[m_nH0Offset + 37] =  0.000000003180;
		for( i=1; i<38; i++ )
			m_pH0[m_nH0Offset - i] = m_pH0[m_nH0Offset + i];

		// 分解数列H1
		m_nH1Offset = 128;
		m_nH1Begin = -35;
		m_nH1End = 33;
		m_pH1[m_nH1Offset -  1] =  0.614175730928;
		m_pH1[m_nH1Offset +  0] = -0.194994482794;
		m_pH1[m_nH1Offset +  1] = -0.308919674972;
		m_pH1[m_nH1Offset +  2] =  0.143937155733;
		m_pH1[m_nH1Offset +  3] =  0.162242892853;
		m_pH1[m_nH1Offset +  4] = -0.081921323130;
		m_pH1[m_nH1Offset +  5] = -0.086457328462;
		m_pH1[m_nH1Offset +  6] =  0.044448525483;
		m_pH1[m_nH1Offset +  7] =  0.046231289856;
		m_pH1[m_nH1Offset +  8] = -0.023865675342;
		m_pH1[m_nH1Offset +  9] = -0.024740861224;
		m_pH1[m_nH1Offset + 10] =  0.012783812041;
		m_pH1[m_nH1Offset + 11] =  0.013242571125;
		m_pH1[m_nH1Offset + 12] = -0.006844036508;
		m_pH1[m_nH1Offset + 13] = -0.007088378368;
		m_pH1[m_nH1Offset + 14] =  0.003663635980;
		m_pH1[m_nH1Offset + 15] =  0.003794216253;
		m_pH1[m_nH1Offset + 16] = -0.001961130910;
		m_pH1[m_nH1Offset + 17] = -0.002030887633;
		m_pH1[m_nH1Offset + 18] =  0.001049835002;
		m_pH1[m_nH1Offset + 19] =  0.001086943371;
		m_pH1[m_nH1Offset + 20] = -0.000562102051;
		m_pH1[m_nH1Offset + 21] = -0.000581537629;
		m_pH1[m_nH1Offset + 22] =  0.000301154337;
		m_pH1[m_nH1Offset + 23] =  0.000310759092;
		m_pH1[m_nH1Offset + 24] = -0.000161709823;
		m_pH1[m_nH1Offset + 25] = -0.000165359122;
		m_pH1[m_nH1Offset + 26] =  0.000087504723;
		m_pH1[m_nH1Offset + 27] =  0.000086672331;
		m_pH1[m_nH1Offset + 28] = -0.000048571352;
		m_pH1[m_nH1Offset + 29] = -0.000042932980;
		m_pH1[m_nH1Offset + 30] =  0.000028965983;
		m_pH1[m_nH1Offset + 31] =  0.000016343928;
		m_pH1[m_nH1Offset + 32] = -0.000018868339;
		m_pH1[m_nH1Offset + 33] =  0.000004717085;
		for( i=0; i<34; i++ )
			m_pH1[m_nH1Offset - 2 - i] = m_pH1[m_nH1Offset + i];

		// 再構成数列G0
		m_nG0Offset = 128;
		m_nG0Begin = -2;
		m_nG0End = 2;
		m_pG0[m_nG0Offset -  2] =  0.088388347648;
		m_pG0[m_nG0Offset -  1] =  0.353553390593;
		m_pG0[m_nG0Offset +  0] =  0.530330085890;
		m_pG0[m_nG0Offset +  1] =  0.353553390593;
		m_pG0[m_nG0Offset +  2] =  0.088388347648;

		// 再構成数列G1
		m_nG1Offset = 128;
		m_nG1Begin = -4;
		m_nG1End = 6;
		m_pG1[m_nG0Offset -  4] = -0.000059579244;
		m_pG1[m_nG0Offset -  3] =  0.007387826298;
		m_pG1[m_nG0Offset -  2] = -0.099914392758;
		m_pG1[m_nG0Offset -  1] =  0.470914347264;
		m_pG1[m_nG0Offset +  0] = -1.101143593894;
		m_pG1[m_nG0Offset +  1] =  1.445630784669;
		m_pG1[m_nG0Offset +  2] = -1.101143593894;
		m_pG1[m_nG0Offset +  3] =  0.470914347264;
		m_pG1[m_nG0Offset +  4] = -0.099914392758;
		m_pG1[m_nG0Offset +  5] =  0.007387826298;
		m_pG1[m_nG0Offset +  6] = -0.000059579244;

		// オリジナルのΦ(x)
		m_nSc1Offset = 128;
		m_nSc1Begin = -1;
		m_nSc1End = 1;
		m_pSc1[m_nSc1Offset - 1] = 0.240721868973;
		m_pSc1[m_nSc1Offset + 0] = 0.962887475890;
		m_pSc1[m_nSc1Offset + 1] = 0.240721868973;

		// 原点を頂点とするΦ(x) (オリジナルと同じ)
		m_nSc2Offset = 128;
		m_nSc2Begin = -1;
		m_nSc2End = 1;
		m_pSc2[m_nSc2Offset - 1] = 0.240721868973;
		m_pSc2[m_nSc2Offset + 0] = 0.962887475890;
		m_pSc2[m_nSc2Offset + 1] = 0.240721868973;

		// 原点を頂点とするΦ(x)による補間
		m_nBetaSc2Offset = 128;
		m_nBetaSc2Begin = -17;
		m_nBetaSc2End = 17;
		m_pBetaSc2[m_nBetaSc2Offset +  0] =  1.199206086582;
		m_pBetaSc2[m_nBetaSc2Offset +  1] = -0.321326302458;
		m_pBetaSc2[m_nBetaSc2Offset +  2] =  0.086099123251;
		m_pBetaSc2[m_nBetaSc2Offset +  3] = -0.023070190544;
		m_pBetaSc2[m_nBetaSc2Offset +  4] =  0.006181638925;
		m_pBetaSc2[m_nBetaSc2Offset +  5] = -0.001656365158;
		m_pBetaSc2[m_nBetaSc2Offset +  6] =  0.000443821706;
		m_pBetaSc2[m_nBetaSc2Offset +  7] = -0.000118921668;
		m_pBetaSc2[m_nBetaSc2Offset +  8] =  0.000031864965;
		m_pBetaSc2[m_nBetaSc2Offset +  9] = -0.000008538192;
		m_pBetaSc2[m_nBetaSc2Offset + 10] =  0.000002287802;
		m_pBetaSc2[m_nBetaSc2Offset + 11] = -0.000000613014;
		m_pBetaSc2[m_nBetaSc2Offset + 12] =  0.000000164256;
		m_pBetaSc2[m_nBetaSc2Offset + 13] = -0.000000044011;
		m_pBetaSc2[m_nBetaSc2Offset + 14] =  0.000000011789;
		m_pBetaSc2[m_nBetaSc2Offset + 15] = -0.000000003144;
		m_pBetaSc2[m_nBetaSc2Offset + 16] =  0.000000000786;
		for( i=1; i<17; i++ )
			m_pBetaSc2[m_nBetaSc2Offset - i] = m_pBetaSc2[m_nBetaSc2Offset + i];

		break;

	default:

		break;
	}

	return true;
}
