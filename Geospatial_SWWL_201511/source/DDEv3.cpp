/* Using Random sampling intervals. @Feb 18,2013 by Riqi Su*/
/* Change the mat and array storage, automatically detect the network size,
									@Feb 19, 2013 by Riqi Su*/
/* */

#include "dataioV1.h"

int Taumax=5000;				// largest tau
long Step=1000000;		// at least 5* 10^7
int DiffInt=2000;
int MinInt= 5000;

FILE *fp,*fp0,*fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7,*fp8,*fp9,*fp10, *fp11;
double Dt=5e-5;
double Speed=1;
int minDelay=40;

int **adjacency, **Amatrix, *degree, **tau_link, **tau_res;
double **weight, **k1, **k2, **k3, **k4, **xt, **xtcom, **xlagcom;
double **xtlag, **k1lag, **k2lag, **k3lag;
double **ppRKTIME, **ppWei, **ppDelay;

double RsA_G=0.2;
double RsB_G=0.2;
double RsC_G=5.7;

int N;
int t_res;   // rescaled t by Taumax length

int NeedD=500; // default measurement number;

void input(FILE *fp)
{
	int m,n,i,j;  //,j
	double temp1, temp2;
	int maxID=0, minID=1000;
	int link = 0;
	double maxDel=0, minDel=100;
	while(fscanf(fp,"%d\t%d\t%lf\t%lf\n",&m, &n, &temp1, &temp2)==4)
	{
		if ( m> maxID)	maxID=m;
		else if ( m< minID)	minID=m;

		if ( n> maxID)	maxID=n;
		else if ( n< minID)	minID=n;
	
		if ( temp2> maxDel)	maxDel=temp2;
		if ( temp2< minDel)	minDel=temp2;
		
		link++;
	}
	
	Speed= minDelay*Dt/minDel;
	Taumax= int( maxDel* Speed/Dt +20);

	printf("minD is %lf\tmaxD is %lf\nTaumax is %d\tSpeed is %lf\n", 
		minDel, maxDel, Taumax, Speed);
	
	N= maxID- minID +1;

	fseek( fp, 0,0);

	adjacency= genMatInt( N, N);
	weight= genMatDou(N, N);
	Amatrix= genMatInt(N, N);
	degree= genArrayInt(N);

	k1= genMatDou(N,3);
	k2= genMatDou(N,3);
	k3= genMatDou(N,3);
	k4= genMatDou(N,3);
	xt= genMatDou(N,3);

	xtcom= genMatDou(N,3);
	xlagcom= genMatDou(N,Taumax);
	tau_link= genMatInt(N,Taumax);
	tau_res= genMatInt(N,Taumax);

	xtlag= genMatDou( N, Taumax+10);
	k1lag= genMatDou( N, Taumax+10);
	k2lag= genMatDou( N, Taumax+10);
	k3lag= genMatDou( N, Taumax+10);

	ppWei=genMatDou( N, N);
	ppDelay= genMatDou( N,N);

	ppRKTIME=genMatDou( NeedD*3, 3*N);

	for(i=0;i<N;i++){
		degree[i]=0;
		for(j=0; j<N; j++){
			Amatrix[i][j] = 0;
			ppWei[i][j]=0;
			ppDelay[i][j]=0;
		}
	}

	link=0;
	while(fscanf(fp,"%d\t%d\t%lf\t%lf\n",&m, &n, &temp1, &temp2)==4)
	{
		adjacency[m][degree[m]] = n;
		Amatrix[m][n] = 1;

		link++;

		temp2= temp2*Speed;
		ppWei[m][n]= temp1;
		ppDelay[m][n]= temp2;
//		ppWei[n][m]= temp1;
//		ppDelay[n][m]= temp2;

		weight[m][ degree[m]]= temp1;
//		weight[n][ degree[n]]= temp1;
		tau_link[m][ degree[m]]= (int)(temp2/Dt);
//		tau_link[n][ degree[n]]= (int)(temp2/Dt);
		degree[m]++;
//		degree[n]++;
	}

	fclose(fp);
}

void initialize()    // initial angle of each oscillator
{
	int i,j;
	double temp;
	for(i=0; i<N; i++){
		for(j=0; j<3; j++){
			temp = f_RandC();
			xt[i][j] = temp;
		}
	for(j=0; j<Taumax+10; j++){
			xtlag[i][j] = xt[i][2];  // initialize xlag record matrix
			k1lag[i][j] = xt[i][2];
			k2lag[i][j] = xt[i][2];
			k3lag[i][j] = xt[i][2];
		}

	}
	for(i=0; i<N; i++){
		for(j=0; j<degree[i]; j++){
			tau_res[i][j] = 0;
		}
	}

}

double xfunction(int node, int t)
{

	int i; //, neighbor;
	double sum = 0.0;
	double value;

	for(i=0; i<degree[node]; i++){
			sum += weight[node][i] * ( xlagcom[node][i] - xtcom[node][2] );   // time-delayed coupling  ?---> 1
	}
	value = sum  - xtcom[node][1] - xtcom[node][2];
	return(value);
}

double yfunction(int node)
{
	double value ;
	value = xtcom[node][0] + RsA_G * xtcom[node][1] ;
	return(value);
}


double zfunction(int node)
{
	double value;
	value = RsB_G + xtcom[node][2] * (xtcom[node][0] - RsC_G);
	return(value);
}


void DDE_RK4( int t)
{
	int i,j;
	int neighbor;
	t_res = t % Taumax;   // position of t in the lag matrices

	for(i=0; i<N; i++){
		for(j=0; j<degree[i]; j++){
			tau_res[i][j] = (t % Taumax + Taumax - tau_link[i][j]) % Taumax;  // tau_res is position of each link in lag matrices
		}
	}

	for(i=0; i<N; i++){   // only x dimesion time delay
		for(j=0; j<3; j++){
			xtcom[i][j] = xt[i][j];
		}
		for(j=0; j<degree[i]; j++){
			neighbor = adjacency[i][j];
			xlagcom[i][j] = xtlag[neighbor][tau_res[i][j]];
		}
	}

	for(i=0; i<N; i++){
		k1[i][0] = xfunction(i,t);
		k1[i][1] = yfunction(i);
		k1[i][2] = zfunction(i);
	}

	for(i=0; i<N; i++){
		for(j=0; j<3; j++){
			xtcom[i][j] = xt[i][j] + k1[i][j] * Dt / 2.0;
		}
		for(j=0; j<degree[i]; j++){
			neighbor = adjacency[i][j];
			xlagcom[i][j] = xtlag[neighbor][tau_res[i][j]] + k1lag[neighbor][tau_res[i][j]] * Dt / 2.0;
		}
	}

	for(i=0; i<N; i++){
		k2[i][0] = xfunction(i,t);
		k2[i][1] = yfunction(i);
		k2[i][2] = zfunction(i);
	}

	for(i=0; i<N; i++){
		for(j=0; j<3; j++){
			xtcom[i][j] = xt[i][j] + k2[i][j] * Dt / 2.0;
		}
		for(j=0; j<degree[i]; j++){
			neighbor = adjacency[i][j];
			xlagcom[i][j] = xtlag[neighbor][tau_res[i][j]] + k2lag[neighbor][tau_res[i][j]] * Dt / 2.0;
		}
	}

	for(i=0; i<N; i++){
		k3[i][0] = xfunction(i,t);
		k3[i][1] = yfunction(i);
		k3[i][2] = zfunction(i);
	}

	for(i=0; i<N; i++){
		for(j=0; j<3; j++){
			xtcom[i][j] = xt[i][j] + k3[i][j] * Dt;
		}
		for(j=0; j<degree[i]; j++){
			neighbor = adjacency[i][j];
			xlagcom[i][j] = xtlag[neighbor][tau_res[i][j]] + k3lag[neighbor][tau_res[i][j]] * Dt ;
		}
	}

	for(i=0; i<N; i++){
		k4[i][0] = xfunction(i,t);
		k4[i][1] = yfunction(i);
		k4[i][2] = zfunction(i);

		for(j=0; j<3; j++){
			xt[i][j] += Dt/6.0 * (k1[i][j] + 2*k2[i][j] + 2*k3[i][j] + k4[i][j]);  // final
		}

	}

//// for lag arrays

	int dim = 2;  // 3--->1  coupling
	for(i=0; i<N; i++){
		xtlag[i][t_res] = xt[i][dim];
		k1lag[i][t_res] = k1[i][dim];
		k2lag[i][t_res] = k2[i][dim];
		k3lag[i][t_res] = k3[i][dim];

	}

}

void record(int t)
{
	int i;
	for(i=0;i<N; i++){
		ppRKTIME[t][3*i] = xt[i][0];
		ppRKTIME[t][3*i+1] = xt[i][1];
		ppRKTIME[t][3*i+2] = xt[i][2];
	}
}


int  f_Output( char *pName)
{
	long i, j;

	MATFILE matfile = MATopen( pName);

	double **Xdata= genMatDou( NeedD, N);
	double **Ydata= genMatDou( NeedD, N);
	double **Zdata= genMatDou( NeedD, N);
	for ( i = 0; i < NeedD; i += 1)
	{
		for ( j = 0; j < N; j += 1)
		{
			Xdata[i][j]= ppRKTIME[i*3+1][j*3];
			Ydata[i][j]= ppRKTIME[i*3+1][j*3 + 1];
			Zdata[i][j]= ppRKTIME[i*3+1][j*3 + 2];
		}
	}
	MATadd( matfile, "Xt", Xdata, NeedD, N);
	MATadd( matfile, "Yt", Ydata, NeedD, N);
	MATadd( matfile, "Zt", Zdata, NeedD, N);

	freeMatDou( Xdata);
	freeMatDou( Ydata);
	freeMatDou( Zdata);

	double **ppRK_H = genMatDou(1,1);
	ppRK_H[0][0] = Dt;
	MATadd( matfile, "RK_H", ppRK_H, 1,1);
	ppRK_H[0][0] = Speed;
	MATadd( matfile, "Speed", ppRK_H, 1,1);
	freeMatDou( ppRK_H);
	
	MATadd( matfile, "iniWmatrix", ppWei, N, N);
	MATadd( matfile, "DelayMat", ppDelay, N, N);

	double **dXdt= genMatDou( NeedD, N);
	double **dYdt= genMatDou( NeedD, N);
	double **dZdt= genMatDou( NeedD, N);
	double Dt2=2*Dt;
	for ( i = 0; i < NeedD; i += 1)
	{
		for ( j = 0; j < N; j += 1)
		{
			dXdt[i][j]= ( ppRKTIME[i*3 + 2][j*3] - ppRKTIME[i*3][j*3] ) / Dt2;
			dYdt[i][j]= ( ppRKTIME[i*3 + 2][j*3 + 1] - ppRKTIME[i*3][j*3 + 1] ) / Dt2;
			dZdt[i][j]= ( ppRKTIME[i*3 + 2][j*3 + 2] - ppRKTIME[i*3][j*3 + 2] ) / Dt2;
		}
	}
	MATadd( matfile, "dXdt", dXdt, NeedD, N);
	MATadd( matfile, "dYdt", dYdt, NeedD, N);
	MATadd( matfile, "dZdt", dZdt, NeedD, N);

	freeMatDou( dXdt);
	freeMatDou( dYdt);
	freeMatDou( dZdt);

	MATadd( matfile, "TimeSer", ppRKTIME, NeedD*3, N*3);

	double **ppTempD = genMatDou(1,1);
	ppTempD[0][0]= RsA_G;
	MATadd( matfile, "Rs_a", ppTempD, 1,1);
	ppTempD[0][0]= RsB_G;
	MATadd( matfile, "Rs_b", ppTempD, 1,1);
	ppTempD[0][0]= RsC_G;
	MATadd( matfile, "Rs_c", ppTempD, 1,1);

	freeMatDou( ppTempD);

	MATclose( matfile);
	return 1;
}


int main(int argc, char *argv[])
{
	int i; //,k
	int count;
	char TimeSerName[100];
	f_SRand(time(NULL));

	if( argc >= 3){
	      NeedD =atoi( argv[2]);
	}

	FILE *fp=fopen( argv[1], "r");
	input(fp);
	initialize();

	sprintf( TimeSerName, "TimeSer%s", argv[1]);
	char *pChar= strstr( TimeSerName, ".txt");
	if (pChar==NULL){	
		strcat(TimeSerName, ".mat");
	}
	else{
		pChar= strcpy(pChar, ".mat");		
	}


	long NextRec, nex;
	count = 0;
	i=0;

	for (i=0; i< 200000; i++){
		DDE_RK4(i);	
	}

	while (count< NeedD*3)	{
		NextRec= random()%DiffInt + MinInt;
		for ( nex = 0; nex < NextRec; nex += 1)	{
			DDE_RK4(i);
			i++;
		}

		for ( nex = 0; nex < 3; nex += 1)	{
			DDE_RK4(i);
			record( count);
			count++;
			i++;
		}
		//printf( "%d\t%ld\n", count/3, NextRec);
	}

	f_Output( TimeSerName);

	return 0;
}
