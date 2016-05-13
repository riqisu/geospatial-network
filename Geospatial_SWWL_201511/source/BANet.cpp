#include "dataioV1.h"

int Nodes=20;
	
double f_Tau(int i, int j, double **ppPos)
{
	double tau;
	tau=sqrt(  (ppPos[0][i]- ppPos[0][j]) *(ppPos[0][i]- ppPos[0][j]) 
			+  (ppPos[1][i]- ppPos[1][j]) *(ppPos[1][i]- ppPos[1][j]) );
	return tau;
}

int f_UpdateCumu( double *pCumu, int *pDeg, int maxID)
{
	int i, sum=0;
	
	for ( i = 0; i < maxID; i += 1){
		sum+= pDeg[i];
	}
	
	pCumu[0]=pDeg[0]*1.0/sum;
	for ( i = 1; i < maxID; i += 1){
		pCumu[i]=pDeg[i]*1.0/sum + pCumu[i-1];
	}
	
	return 1;
}

int f_findNei(double *pCumu, int maxID)
{
	double prob=f_RandO();
	int Left=0, Right=maxID, Mid;
	
	if ( prob< pCumu[0]){
		return 0;
	}
	
	Mid= ( Left+ Right)/2;
	while (Mid > Left){
		if ( prob> pCumu[Mid]){
			Left= Mid;
		}
		else	{
			Right= Mid;
		}
		Mid= ( Left+ Right)/2;
	}
	
//	int i;
//	printf("Prob=%0.4lf\tMid=%d\tp[Mid]=%0.3lf\tp[Mid-1]=%0.3lf\n",prob, Mid, pCumu[Mid], pCumu[Mid-1]);
	return (Mid+1);
}

int main (int argc, char const* argv[])
{
	int Degree=3;
	int IniNode=4;
	int GridNum=12;
	double tau, wei, speed=0.01;
	double WeiDiff, WeiMax=0.05, WeiMin=0.01;

	f_SRand( time(NULL));
	char fileName[100];
	int Stamp= (int) (f_RandO()*100);
	
	switch (argc)	{
		case 4: Stamp=atoi( argv[3]);
		case 3: Degree= atoi( argv[2]);
		case 2: Nodes= atoi( argv[1]);
				break;
		default:	break;
	}
	
	WeiDiff= WeiMax- WeiMin;
	if ( IniNode < Degree){
		IniNode= Degree;
	}
	
	sprintf( fileName, "BAN%dk%dStam%d.txt", Nodes, Degree, Stamp);
	FILE *fp= fopen( fileName, "w");
	
	sprintf( fileName, "PosN%dk%dStam%d.txt", Nodes, Degree, Stamp);
	FILE *fpos= fopen( fileName, "w");
	
	int **ppNet= genMatInt( Nodes, Nodes);
	int i, j;
	int addNode=0, Select=0;
	
	double **ppPos=genMatDou(2, Nodes);
	
	int *pDeg= genArrayInt(Nodes);
	double *pCumu=genArrayDou( Nodes);
	
	for ( i = 0; i < Nodes; i += 1) {
		pDeg[i]=0;
		for ( j = 0; j < Nodes; j += 1) {
			ppNet[i][j]=0;
		}
	}

	int AllGrid= GridNum* GridNum;	
	while ( AllGrid< 5*Nodes){ 
		// check if the Grid larger than nodes;
		GridNum++;
		AllGrid= GridNum* GridNum;
	}
	
	double xGrd=1.0/GridNum, yGrd=1.0/GridNum;
	int SelGrd=0;
	int *SelList= genArrayInt(AllGrid);
	for ( i = 0; i < AllGrid; i += 1) {
		SelList[i]=0;
	}
	
	for ( i = 0; i < Nodes; i += 1)
	{
		do{
		SelGrd= (int)( f_RandO()*AllGrid);
		}while (SelList[SelGrd]==1);
		
		SelList[SelGrd]=1;
		ppPos[0][i]= (SelGrd/GridNum)*xGrd +xGrd*f_RandO();
		ppPos[1][i]= (SelGrd%GridNum)*yGrd+yGrd*f_RandO();
		fprintf(fpos, "%lf\t%lf\n", ppPos[0][i], ppPos[1][i]);
	}
	fclose( fpos);
	
	for ( i = 0; i < IniNode; i += 1)
	{
		pDeg[i]=IniNode-1;
		for ( j = i+1; j < IniNode; j += 1)
		{
			ppNet[i][j]=1;
			ppNet[j][i]=1;
			tau= speed* f_Tau( i, j, ppPos);
			
			wei= f_RandO()*WeiDiff + WeiMin;
			fprintf(fp, "%d\t%d\t%lf\t%0.8lf\n", i, j, wei, tau);
			wei= f_RandO()*WeiDiff + WeiMin;
			fprintf(fp, "%d\t%d\t%lf\t%0.8lf\n", j, i, wei, tau);
		}
	}
	
	for ( i = IniNode; i < Nodes; i += 1)
	{
		f_UpdateCumu( pCumu, pDeg, i); // prepare the prob.
		
		addNode= Degree; // Add Degree links;
		while ( addNode>0)
		{
			Select= f_findNei(pCumu, i);
			if ( (ppNet[i][Select] !=1)&&(Select!=-1) )
			{
				ppNet[i][Select]=1;
				ppNet[Select][i]=1;
				pDeg[Select]++;
				tau= speed* f_Tau( i,  Select, ppPos);
				wei= f_RandO()*WeiDiff + WeiMin;
				fprintf(fp, "%d\t%d\t%lf\t%0.8lf\n", i, Select, wei, tau);
				wei= f_RandO()*WeiDiff + WeiMin;
				fprintf(fp, "%d\t%d\t%lf\t%0.8lf\n", Select, i, wei, tau);
				addNode--;
			}
		}
		pDeg[i]+= Degree;
	}
	
	fclose(fp);
	return 0;
}
