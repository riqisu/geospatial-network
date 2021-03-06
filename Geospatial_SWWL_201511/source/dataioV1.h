/********************************************************************\
 *                     Basic Data I/O Functions                     *
 * ---------------------------------------------------------------- *
 * OS: Unix / Linux / MacOSX / Windows                              *
 * Author: Xuan Ni                                                  *
 * Create: Thu Feb 18 22:09:08 MST 2010	                            *
 * Modify: Web Feb 13 11:02:14 MST 2013 by Riqi Su                  *
 * ---------------------------------------------------------------- *
 * CONTENTS                                                         *
 *   - Generate random variables                                    *
 *   - Generate matrix and arrays;                                  *
 *   - Write matrix to MATLAB compatible MAT-file                   *
 * More coming ...                                                  *
 * ---------------------------------------------------------------- *
 * (1) File length in bytes
 *
 *     long filelen(const char *FILENAME);
 *
 * (2) Binary file IO
 *
 *     (a) Read binary file
 *
 *         short          *readbins (const char *FILENAME,
 *                                   long       *DATALEN);
 *         unsigned short *readbinus(const char *FILENAME,
 *                                   long       *DATALEN);
 *         int            *readbini (const char *FILENAME,
 *                                   long       *DATALEN);
 *         unsigned       *readbinui(const char *FILENAME,
 *                                   long       *DATALEN);
 *         float          *readbinf (const char *FILENAME,
 *                                   long       *DATALEN);
 *         double         *readbind (const char *FILENAME,
 *                                   long       *DATALEN);
 *
 *         RETURN: pointer of specific type to array in memory,
 *                 NULL on error
 *
 *     (b) Write binary file
 *
 *         int writebins (const char     *FILENAME,
 *                        short          *DATA,
 *                        long            DATALEN);
 *         int writebinus(const char     *FILENAME,
 *                        unsigned short *DATA,
 *                        long            DATALEN);
 *         int writebini (const char     *FILENAME,
 *                        int            *DATA,
 *                        long            DATALEN);
 *         int writebinui(const char     *FILENAME,
 *                        unsigned       *DATA,
 *                        long            DATALEN);
 *         int writebinf (const char     *FILENAME,
 *                        float          *DATA,
 *                        long            DATALEN);
 *         int writebind (const char     *FILENAME,
 *                        double         *DATA,
 *                        long            DATALEN);
 *
 *         RETURN: 0 on success, -1 on error
 *
 *     (c) Other Read/Write
 *         There are two more functions, which are the called by
 *         above functions, and should not be called in most cases.
 *
 *         void *readbin(const char *FILENAME,
 *                       long       *BYTELEN);
 *         int  writebin(const char *FILENAME,
 *                       void       *DATA,
 *                       long        BYTELEN);
 *
 *     (d) Free memory return by readbin*
 *
 *         void freebin( void *DATA );
 *
 *         This function should be called for each readbin function
 *
 * (3) Text matrix file IO
 *
 *     (a) double **genmat(int ROW, int COL);
 *
 *         generate a double precision matrix pointer
 *
 *     (b) double **readmat(const char *FILENAME,
 *                          int        *ROW,
 *                          int        *COL);
 *
 *         read a matrix from a text file, the row and column of
 *         the matrix will be stored in two integers pointed by
 *         ROW and COL, respectively. This function can read matrix
 *         delimited by any occurance or combination of space,
 *         comma, and tab. MATLAB "dlmwrite" and "csvwrite"
 *         compatible.
 *
 *     (c) int writemat(const char  *FILENAME,
 *                      double     **MAT,
 *                      int          ROW,
 *                      int          COL);
 *
 *         write a matrix of ( ROW x COL ) in memory to a text file
 *
 *     (d) void freemat(double **MAT);
 *
 *         free matrix memory generated by "genmat()" or "readmat()"
 *         This function should be called to avoid memory leak
 *
 * (4) MATLAB level 5 MAT-file generator (no compression)
 *     Currently only write MAT-file is implemented
 *     The generated MAT-file can be read by any MATLAB version
 *     that support level 5 MAT-file, tested on MATLAB 2009b (Ubuntu
 *     9.10/MacOSX SnowLeopard) and MATLAB 7.0 (WinXP)
 *
 *     (a) MATFILE MATopen(const char *FILENAME);
 *
 *         open a MAT-file for write
 *
 *     (b) int MATadd(MATFILE        FILEHANDLE,
 *                    const char    *VARIABLENAME,
 *                    double       **MATRIX,
 *                    int            ROW,
 *                    int            COL );
 *
 *         add matrix to MAT-file
 *
 *     (c) void MATclose(MATFILE  FILEHANDLE);
 *
 *         close MAT-file
\********************************************************************/


#ifndef DATAIO_H
#define DATAIO_H

/* determine OS type */
#ifndef dataio_OS
#if 	defined(unix)	|| defined(__unix)	|| defined(__unix__)	||\
	defined(linux)	|| defined(__linux)	|| defined(__linux__)	||\
	defined(__MACOSX__)			|| defined(__APPLE__)	||\
	defined(sun)	|| defined(__sun)	|| \
	defined(BSD)	|| defined(__OpenBSD__)	|| defined(__FreeBSD__)
#define dataio_OS 0
#elif	defined(WIN32)	|| defined(_WIN32)	|| defined(__WIN32__) 	||\
	defined(WIN64)	|| defined(_WIN64)	|| defined(__WIN64__)	||\
	defined(_MSC_VER)
#define dataio_OS 1
#else
#define dataio_OS 2
#endif
#endif /* dataio_OS */


/* include std headers */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


/* POSIX headers for performance consideration */
#if dataio_OS==0
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/utsname.h>
#endif

#include <math.h>

/****************************************\
 *    Generate random number            *
\****************************************/
/* Period parameters */  
#define RANDNEW_N 624
#define RANDNEW_M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */


static unsigned long RANDNEW_mt[RANDNEW_N]; /* the array for the state vector  */
static int RANDNEW_mti=RANDNEW_N+1; /* RANDNEW_mti==RANDNEW_N+1 means RANDNEW_mt[RANDNEW_N] is not initialized */

void f_SRand(unsigned long s)
{
    RANDNEW_mt[0]= s & 0xffffffffUL;
    for (RANDNEW_mti=1; RANDNEW_mti<RANDNEW_N; RANDNEW_mti++) {
        RANDNEW_mt[RANDNEW_mti] = 
	    (1812433253UL * (RANDNEW_mt[RANDNEW_mti-1] ^ (RANDNEW_mt[RANDNEW_mti-1] >> 30)) + RANDNEW_mti); 
        RANDNEW_mt[RANDNEW_mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* generates a random number on [0,1]-real-interval */
double f_RandC(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (RANDNEW_mti >= RANDNEW_N) { /* generate RANDNEW_N words at one time */
        int kk;

        if (RANDNEW_mti == RANDNEW_N+1)   /* if init_genrand() has not been called, */
            f_SRand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<RANDNEW_N-RANDNEW_M;kk++) {
            y = (RANDNEW_mt[kk]&UPPER_MASK)|(RANDNEW_mt[kk+1]&LOWER_MASK);
            RANDNEW_mt[kk] = RANDNEW_mt[kk+RANDNEW_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<RANDNEW_N-1;kk++) {
            y = (RANDNEW_mt[kk]&UPPER_MASK)|(RANDNEW_mt[kk+1]&LOWER_MASK);
            RANDNEW_mt[kk] = RANDNEW_mt[kk+(RANDNEW_M-RANDNEW_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (RANDNEW_mt[RANDNEW_N-1]&UPPER_MASK)|(RANDNEW_mt[0]&LOWER_MASK);
        RANDNEW_mt[RANDNEW_N-1] = RANDNEW_mt[RANDNEW_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        RANDNEW_mti = 0;
    }
  
    y = RANDNEW_mt[RANDNEW_mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y*(1.0/4294967295.0);
}

/* generates a random number on (0,1)-real-interval */
double f_RandO(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (RANDNEW_mti >= RANDNEW_N) { /* generate RANDNEW_N words at one time */
        int kk;

        if (RANDNEW_mti == RANDNEW_N+1)   /* if f_SRand() has not been called, */
            f_SRand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<RANDNEW_N-RANDNEW_M;kk++) {
            y = (RANDNEW_mt[kk]&UPPER_MASK)|(RANDNEW_mt[kk+1]&LOWER_MASK);
            RANDNEW_mt[kk] = RANDNEW_mt[kk+RANDNEW_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<RANDNEW_N-1;kk++) {
            y = (RANDNEW_mt[kk]&UPPER_MASK)|(RANDNEW_mt[kk+1]&LOWER_MASK);
            RANDNEW_mt[kk] = RANDNEW_mt[kk+(RANDNEW_M-RANDNEW_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (RANDNEW_mt[RANDNEW_N-1]&UPPER_MASK)|(RANDNEW_mt[0]&LOWER_MASK);
        RANDNEW_mt[RANDNEW_N-1] = RANDNEW_mt[RANDNEW_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        RANDNEW_mti = 0;
    }
  
    y = RANDNEW_mt[RANDNEW_mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    
   return (((double)y) + 0.5)*(1.0/4294967296.0);
}

double f_RandN()	
{				        
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * f_RandC() - 1.0;
			x2 = 2.0 * f_RandC() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return y1;
}


/****************************************\

 *    Generate array & matrix            *
\****************************************/


int **genMatInt(int row, int col ){
	int **mat = (int **)malloc(sizeof(int*)*row);
	int *p    = (int *) malloc(sizeof(int )*row*col);
	int i;
	for(i=0;i<row;i++,p+=col){
		mat[i]=p;
	}
	return mat;
}

void freeMatInt( int **mat)
{
	free( mat[0]);
	free( mat);
}

int *genArrayInt( int col){
	int *array= ( int *)malloc( sizeof(int)* col);
	return array;
}

void freeArrayInt( int *array)
{
	free( array);
}

double **genMatDou(int row, int col ){
	double **mat = (double **)malloc(sizeof(double*)*row);
	double *p    = (double *) malloc(sizeof(double )*row*col);
	int i;
	for(i=0;i<row;i++,p+=col){
		mat[i]=p;
	}
	return mat;
}

void freeMatDou( double **mat)
{
	free( mat[0]);
	free( mat);
}

double *genArrayDou( int col)
{
	double *array= ( double *)malloc( sizeof(double)*col);
	return array;
}

void freeArrayDou( double *array)
{
	free( array);
}

/****************************************\
 *    MATLAB level 5 MAT-file write     *
\****************************************/
/* define data types used in matlab MAT-file level 5 */
#define miINT8 		1
#define miUINT8		2
#define miINT16 	3
#define miUINT16	4
#define miINT32		5
#define miUINT32 	6
#define miSINGLE 	7
#define miRESRV1	8
#define miDOUBLE  	9
#define miRESRV2	10
#define miRESRV3 	11
#define miINT64		12
#define miUINT64	13
#define miMATRIX	14
#define miCOMPRESSED	15
#define miUTF8		16
#define miUTF16		17
#define miUTF32		18
/* define array types */
#define mxCELL_CLASS	1
#define mxSTRUCT_CLASS 	2
#define mxOBJECT_CLASS 	3
#define mxCHAR_CLASS	4
#define mxSPARSE_CLASS	5
#define mxDOUBLE_CLASS 	6
#define mxSINGLE_CLASS 	7
#define mxINT8_CLASS	8
#define mxUINT8_CLASS 	9
#define mxINT16_CLASS 	10
#define mxUINT16_CLASS	11
#define mxINT32_CLASS	12
#define mxUINT32_CLASS	13

#define HEADERLEN 	128

#if dataio_OS==0
typedef int MATFILE;
#else
typedef FILE *MATFILE;
#endif

MATFILE MATopen(const char *filename){
	#if dataio_OS==0
	int fd;
	int temp;
	if((fd=open(filename,O_WRONLY|O_CREAT|O_TRUNC,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH))<0){
		printf("ERROR: open file: %s\n",filename);
		return -1;
	}
	#else
	FILE *fp;
	if(!(fp=fopen(filename,"wb"))){
		printf("ERROR: open file: %s\n",filename);
		return NULL;
	}
	#endif
	char header[125];
	short VERSION=0x0100;
	char ENDIAN1='I';
	char ENDIAN2='M';
	time_t t=time(0);
	char *tstr=ctime(&t);
	#if dataio_OS==0
	struct utsname sys;
	uname(&sys);
	sprintf(header,"MATLAB level 5 MAT-file\nPlatform:%s\nAuthor\
: Xuan Ni\nCreate: %s\n",sys.sysname,tstr);
	#else
	sprintf(header,"MATLAB level 5 MAT-file\nPlatform:%s\nAuthor\
: Xuan Ni\nCreate: %s\n","Windows",tstr);
	#endif
	int len=strlen(header);
	int i;
	for(i=len;i<124;i++) header[i]=' ';
	#if dataio_OS==0
	temp = write(fd,header,124);
	temp = write(fd,&VERSION,2);
	temp = write(fd,&ENDIAN1,1);
	temp = write(fd,&ENDIAN2,1);
	return fd;
	#else
	temp = fwrite(header,sizeof(char),124,fp);
	temp = fwrite(&VERSION,sizeof(char),2,fp);
	temp = fwrite(&ENDIAN1,sizeof(char),1,fp);
	temp = fwrite(&ENDIAN2,sizeof(char),1,fp);
	return fp;
	#endif
}

int MATadd(MATFILE file, const char *varname, double **data, int row, int col ){
	#if dataio_OS==0
	if(file<0) return -1;
	#else
	if(!file) return -1;
	#endif
	/* determine the length of varname */
	int len = strlen(varname);
	int varlen = ((len-1)/8+1)*8;
	int temp;
	/* Tag */
	int DataType=miMATRIX;
	int NumBytes=row*col*sizeof(double)+6*8+varlen;
	#if dataio_OS==0
	temp = write(file,&DataType,sizeof(int));
	temp = write(file,&NumBytes,sizeof(int));
	#else
	temp = fwrite(&DataType,sizeof(int),1,file);
	temp = fwrite(&NumBytes,sizeof(int),1,file);
	#endif
	/* array flags */
	DataType = miUINT32;
	NumBytes = 8;
	#if dataio_OS==0
	temp = write(file,&DataType,sizeof(unsigned));
	temp = write(file,&NumBytes,sizeof(int));
	#else
	temp = fwrite(&DataType,sizeof(unsigned),1,file);
	temp = fwrite(&NumBytes,sizeof(int),1,file);
	#endif
	int content=mxDOUBLE_CLASS;
	#if dataio_OS==0
	temp = write(file,&content,sizeof(int));
	#else
	temp = fwrite(&content,sizeof(int),1,file);
	#endif
	content = 0;
	#if dataio_OS==0
	temp = write(file,&content,sizeof(int));
	#else
	temp = fwrite(&content,sizeof(int),1,file);
	#endif
	/* dimensions array */
	DataType = miINT32;
	NumBytes = 8;
	#if dataio_OS==0
	temp = write(file,&DataType,sizeof(int));
	temp = write(file,&NumBytes,sizeof(int));
	temp = write(file,&row,sizeof(int));
	temp = write(file,&col,sizeof(int));
	#else
	temp = fwrite(&DataType,sizeof(int),1,file);
	temp = fwrite(&NumBytes,sizeof(int),1,file);
	temp = fwrite(&row,sizeof(int),1,file);
	temp = fwrite(&col,sizeof(int),1,file);
	#endif
	/* array name */
	DataType = miINT8;
	NumBytes = len;
	#if dataio_OS==0
	temp = write(file,&DataType,sizeof(int));
	temp = write(file,&NumBytes,sizeof(int));
	temp = write(file,varname,len);
	#else
	temp = fwrite(&DataType,sizeof(int),1,file);
	temp = fwrite(&NumBytes,sizeof(int),1,file);
	temp = fwrite(varname,sizeof(char),len,file);
	#endif
	int i,j;
	char c=0;
	for(i=len;i<varlen;i++) {
		#if dataio_OS==0
		temp = write(file,&c,1);
		#else
		temp = fwrite(&c,sizeof(char),1,file);
		#endif
	}
	/* matrix data */
	DataType = miDOUBLE;
	NumBytes = 8*row*col;
	#if dataio_OS==0
	temp = write(file,&DataType,sizeof(int));
	temp = write(file,&NumBytes,sizeof(int));
	for(j=0;j<col;j++){
		for(i=0;i<row;i++){
			temp = write(file,&data[i][j],sizeof(double));
		}
	}
	#else
	temp = fwrite(&DataType,sizeof(int),1,file);
	temp = fwrite(&NumBytes,sizeof(int),1,file);
	for(j=0;j<col;j++){
		for(i=0;i<row;i++){
			temp = fwrite(&data[i][j],sizeof(double),1,file);
		}
	}
	#endif
	return 0;
}

void MATclose(MATFILE file){
	#if dataio_OS==0
	close(file);
	#else
	fclose(file);
	#endif
}



#endif /* DATAIO_H */

