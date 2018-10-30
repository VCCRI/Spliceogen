/*  *  Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
 *  Maryland, U.S.A.  All rights reserved.

 *  sitesk.c compute a score for splice sites based on karlin's paper.
 */

#include  <stdio.h>
#include  <math.h>
#include  <stdlib.h>
#include  <string.h>
#include  <ctype.h>
#include  <limits.h>
#include  <float.h>
#include  <time.h>
#include  <assert.h>

#ifndef PATH_MAX
#define PATH_MAX 255
#endif

extern char TRAIN_DIR[PATH_MAX+1];
extern unsigned char donor_tree;
extern unsigned char acceptor_tree;
extern unsigned char use_cod_noncod;

#define  TRUE  1
#define  FALSE  0
#define  ALPHABET_SIZE  4
#define  ACCEPTOR_LEN  29                    /* Positions +44,72 in a80 */
#define  ACCEPTOR_FILE_NAME "acc"
#define  ACCEPTOR_TREE_FILE "outex"
#define  ACCEPTOR_SIGNAL_OFFSET  24          /* Start of  AG  */
extern double  ACCEPTOR_THRESHOLD;


#define  DONOR_LEN  16                        /* Positions +5,20 in d80 */
#define  DONOR_FILE_NAME "don"
#define  DONOR_TREE_FILE "outin"
#define  DONOR_SIGNAL_OFFSET  5               /* Start of  GT  */
extern double  DONOR_THRESHOLD;

#define  MARKOV_DEGREE  1
#define  MARKOV_LEN  4                     /* ALPHABET_SIZE ^ MARKOV_DEGREE */
#define  LOG_PROB_ACCEPTOR  log (1682.0 / 30191.0)
#define  LOG_PROB_NONACCEPTOR  log (28509.0 / 30191.0)
#define  LOG_PROB_DONOR  log (1682.0 / 30191.0)         /* Change this */
#define  LOG_PROB_NONDONOR  log (28509.0 / 30191.0)     /* Change this */
#define  LOW_SCORE  -99.0  /* Score if pattern does not have GT or AG signal */
#define  RETURN_TRUE_PROB  0

#define CODING_LEN 80

#ifndef  EXIT_FAILURE
  #define  EXIT_FAILURE  -1
#endif
#ifndef  EXIT_SUCCESS
  #define  EXIT_SUCCESS  0
#endif

typedef struct tree {
	int val;
   	int consens;
   	int poz;
   	int no;
   	struct tree *left;
   	struct tree *right;
   } tree;

typedef unsigned int word;

int  Is_Acceptor  (const int *, double *);
int  Is_Donor  (const int *, double *);
int  Acc  (const int *, double *, tree *t,int ind);
int  Don  (const int *, double *, tree *t,int ind);
int comp(const void *a, const void *b);
int findfile(const int * S, tree *t);
void readtree(char *line, tree *t, int start);
int find(char *line, int start);
int  Is_Cod_NonCod  (const int * , double *, int ind);

#define  Start_PosEx 56
#define  Stop_PosEx 84

#define  Start_PosIn 75
#define  Stop_PosIn 90

#define  Start_Cod 0
#define  Stop_Cod 79

#define Start_NoCod 82
#define Stop_NoCod 161


int Is_Acceptor(const int *B, double *Return_Score)
{
  FILE  * Infile;
  double Score,S1,S2;
  static tree *tacc;
  static int readtacc=FALSE;
  char line[5000];
  int i,ind;
  int T[100];
  double score1,score2,score3;
  char treefile[PATH_MAX+1];

  if(!readtacc && acceptor_tree)
    {
   
      /* read the structure of the acceptor tree */

      sprintf(treefile,"%s/%s",TRAIN_DIR,ACCEPTOR_TREE_FILE);
      Infile = fopen (treefile, "r");
      if  (Infile == NULL)
	{
	  fprintf (stderr, "ERROR:  Unable to open file %s\n", ACCEPTOR_TREE_FILE);
	  exit (EXIT_FAILURE);
	}

      tacc = (tree *) malloc(sizeof(tree));
      if (tacc == NULL) {fprintf(stderr," Memory allocation for tree failure.\n"); abort();}
      fgets(line, 5000, Infile);
      i=strlen(line);
      line[i-1]='\0';
      fclose(Infile);

      readtree(line, tacc, 0);
      readtacc=TRUE;
    }
   
  for(i=0;i<=Stop_PosEx-Start_PosEx;i++)
    T[i]=B[i+Start_PosEx];

  ind=Acc(T, &S1, tacc,0);
  if(ind==0) return(0);
  if(acceptor_tree) {
    Acc(T, &S2, tacc,1);
    score1=(S1+S2)/2;
  }
  else { score1=S1; }

  for(i=0;i<=Stop_NoCod-Start_NoCod;i++)
    T[i]=B[i+Start_NoCod];

  Is_Cod_NonCod(T,&score2,0);

  for(i=0;i<=Stop_Cod-Start_Cod;i++)
    T[i]=B[i+Start_Cod];

  Is_Cod_NonCod(T,&score3,1);


  Score=score1+score2+score3;

  *Return_Score=Score;
  
  
  return Score >= ACCEPTOR_THRESHOLD;
	  
      
}  

int Is_Donor(const int *B, double *Return_Score)
{
  FILE  * Infile;
  double Score,S1,S2;
  static tree *tdon;
  static int readtdon=FALSE;
  char line[5000];
  int ind,i;
  int T[100];
  double score1,score2,score3;
  char treefile[PATH_MAX+1];

  if(!readtdon && donor_tree)
    {
   
      /* read the structure of the donor tree */

      sprintf(treefile,"%s/%s",TRAIN_DIR,DONOR_TREE_FILE);
      Infile = fopen (treefile, "r");
      if  (Infile == NULL)
	{
	  fprintf (stderr, "ERROR:  Unable to open file %s\n", DONOR_TREE_FILE);
	  exit (EXIT_FAILURE);
	}

      tdon = (tree *) malloc(sizeof(tree));
      if (tdon == NULL) {fprintf(stderr,"Memory allocation for tree failure.\n"); abort();}
      fgets(line, 5000, Infile);
      i=strlen(line);
      line[i-1]='\0';
      fclose(Infile);

      readtree(line, tdon, 0);
      readtdon=TRUE;
    }

  for(i=0;i<=Stop_PosIn-Start_PosIn;i++)
    T[i]=B[i+Start_PosIn];

  ind=Don(T, &S1, tdon,0);
  if(ind==0) return(0);
  if(donor_tree) {
    Don(T, &S2, tdon,1);
    score1=(S1+S2)/2;
  }
  else {
    score1=S1;
  }

  //  if(score1<=THR_DON) score1=-99;

  for(i=0;i<=Stop_Cod-Start_Cod;i++)
    T[i]=B[i+Start_Cod];

  Is_Cod_NonCod(T,&score2,2);

  //  if(score2<=THR_DON_EX) score2=-99;
  
  for(i=0;i<=Stop_NoCod-Start_NoCod;i++)
    T[i]=B[i+Start_NoCod];

  Is_Cod_NonCod(T,&score3,3);

  //  if(score3<=THR_DON_IN) score3=-99;

  Score=score1+score2+score3;

  *Return_Score=Score;
  
  return Score >= DONOR_THRESHOLD;
	  
      
}  


    
void readtree(char *line, tree *t, int start)
{
 int len;
 int i,n;
 char part[10];
 len=strlen(line);

 i=start;
 while((line[i]=='(')||(line[i]==' ')) i++;
 n=i;
 while(line[i]!=' ')
 {
	part[i-n]=line[i];
	i++;
 }
 part[i-n]='\0';
 t->val=atoi(part);

 i++;
 n=i;
 while(line[i]!=' ')
 { 
	part[i-n]=line[i];
	i++;
 }
 part[i-n]='\0';
 t->consens=atoi(part);

 i++;
 n=i;
 while(line[i]!=' ')
 { 
	part[i-n]=line[i];
	i++;
 }
 part[i-n]='\0';
 t->poz=atoi(part);

 i++;
 n=i;
 while(line[i]!=' ')
 { 
	part[i-n]=line[i];
	i++;
 }
 part[i-n]='\0';
 t->no=atoi(part);

 t->left=NULL;
 t->right=NULL;

 i+=2;n=i;
 if(line[i]=='(') 
 	{
 		i=find(line,i+1);
		t->left = (tree *) malloc(sizeof(tree));
   		if (t->left == NULL) {fprintf(stderr,"Memory allocation for tree failure.\n"); abort();}
        readtree(line,t->left,n);
     }
	
 i+=2;n=i;
 if(line[i]=='(') 
 	{
 		i=find(line,i+1);
		t->right = (tree *) malloc(sizeof(tree));
   		if (t->right == NULL) {fprintf(stderr,"Memory allocation for tree failure.\n"); abort();}
        readtree(line,t->right,n);
     }
}

int find(char *line, int start)
{
 int stop,i;

 i=start;

 while(line[i]!=')')
 	if(line[i]=='(') i=find(line,i+1);
 	else i++;
 stop=i+1;
 return(stop);
}
 	

int comp(const void *a, const void *b)
{ 
  if(*(double *)a > *(double *)b) return(1);
  else if (*(double *)a==*(double *)b) return(0);
  else return(0);

}  
  

int findfile(const int * S, tree *t)
{
	int val, cons, poz;
	val=t->val;

	cons=t->consens;
	if( cons !=-1)
	{ 
		poz=t->poz;
	    if(S[poz]==cons)
	    	val=findfile(S,t->left);
	    else val=findfile(S, t->right);
	}

	return(val);
}


int  Acc  (const int * S, double * Return_Score, tree *t,int ind)

/* Evaluate string  S [0 .. (ACCEPTOR_LEN -1)] and
*  return  TRUE  or  FALSE  as to whether it is a likely acceptor
*  site.  Also set  Return_Score  to the probability that it is an acceptor
*  site. */

  {
   FILE  * Infile;
   static float  Positive_Table[300][ACCEPTOR_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static float  Negative_Table[300][ACCEPTOR_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static int Tables_Loaded[300]={FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE};
   double  Positive_Sum, Negative_Sum, Score;
   char accname[PATH_MAX+1];
#if  RETURN_TRUE_PROB
   double  X, Y;
#endif
   int  i, j, k, Ct, Sub, no;

/* see which acceptor you should use */

if(ind) 
  {
	no=findfile(S,t);
	sprintf(accname,"%s/%s%d",TRAIN_DIR,ACCEPTOR_FILE_NAME,no);

  }
else 
  {
    strcpy(accname,TRAIN_DIR);
    strcat(accname,"/acc1.mar");
    no=0;
  }

   if  (! Tables_Loaded[no])
       {
        Infile = fopen (accname, "r");
        if  (Infile == NULL)
            {
             fprintf (stderr, "ERROR:  Unable to open acceptor file \"%s\"\n",
                        accname);
             exit (EXIT_FAILURE);
            }

        for  (i = MARKOV_DEGREE - 1;  i < ACCEPTOR_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Positive_Table [no][i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading acceptor file \"%s\"\n", 
                                ACCEPTOR_FILE_NAME);
                    exit (EXIT_FAILURE);
                   }
              }

        for  (i = MARKOV_DEGREE - 1;  i < ACCEPTOR_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Negative_Table [no][i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading acceptor file \"%s\"\n", 
                                ACCEPTOR_FILE_NAME);
                    exit (EXIT_FAILURE);
                   }
              }

        fclose (Infile);

        Tables_Loaded[no]  = TRUE;
       }

   if  (S [ACCEPTOR_SIGNAL_OFFSET] != 0
           || S [ACCEPTOR_SIGNAL_OFFSET + 1] != 2)    /* AG */
       {
        * Return_Score = LOW_SCORE;
        return  FALSE;
       }

   Sub = 0;
   for  (i = 0;  i < MARKOV_DEGREE;  i ++)
     Sub = ALPHABET_SIZE * Sub + S [i];

   Positive_Sum = Positive_Table [no][MARKOV_DEGREE - 1] [0] [Sub];
   Negative_Sum = Negative_Table [no][MARKOV_DEGREE - 1] [0] [Sub];

   for  (i = MARKOV_DEGREE;  i < ACCEPTOR_LEN;  i ++)
     {
      j = S [i];
      Positive_Sum += Positive_Table [no] [i] [j] [Sub];
      Negative_Sum += Negative_Table [no] [i] [j] [Sub];
      Sub = ALPHABET_SIZE * (Sub % (MARKOV_LEN / ALPHABET_SIZE)) + j;
     }
  


   Score = Positive_Sum - Negative_Sum;

#if  RETURN_TRUE_PROB
   X = exp (Positive_Sum + LOG_PROB_ACCEPTOR);
   Y = exp (Negative_Sum + LOG_PROB_NONACCEPTOR);
   * Return_Score = log (X / (X + Y));
#else
   * Return_Score = Score;
#endif

   return(1);
  }



int  Don  (const int * S, double * Return_Score, tree *t,int ind)

/* Evaluate string  S [0 .. (DONOR_LEN -1)] and
*  return  TRUE  or  FALSE  as to whether it is a likely donor
*  site.  Also set  Return_Score  to the probability that it is an donor
*  site. */
{
   FILE  * Infile;
   static float  Positive_Table [300][DONOR_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static float  Negative_Table [300][DONOR_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static int Tables_Loaded[300]={FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE};
   double  Positive_Sum, Negative_Sum, Score;
   char donname[PATH_MAX+1];
   int no;

#if  RETURN_TRUE_PROB
   double  X, Y;
#endif
   int  i, j, k, Ct, Sub;

   /* see which donor file you should use */
if(ind)
   { no=findfile(S,t);
    sprintf(donname,"%s/%s%d",TRAIN_DIR,DONOR_FILE_NAME,no);
   }
else 
{
  strcpy(donname,TRAIN_DIR);
  strcat(donname,"/don1.mar");
  no=0;
}

   if  (! Tables_Loaded[no] )
       {

        Infile = fopen (donname, "r");
        if  (Infile == NULL)
            {
             fprintf (stderr, "ERROR:  Unable to open donor file \"%s\"\n",
                        donname);
             exit (EXIT_FAILURE);
            }

        for  (i = MARKOV_DEGREE - 1;  i < DONOR_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Positive_Table [no] [i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading donor file \"%s\"\n", 
                                DONOR_FILE_NAME);
                    exit (EXIT_FAILURE);
                   }
              }

        for  (i = MARKOV_DEGREE - 1;  i < DONOR_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Negative_Table [no] [i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading donor file \"%s\"\n", 
                                DONOR_FILE_NAME);
                    exit (EXIT_FAILURE);
                   }
              }

        fclose (Infile);

        Tables_Loaded [no] = TRUE;
       }


   if  (S [DONOR_SIGNAL_OFFSET] != 2
           || S [DONOR_SIGNAL_OFFSET + 1] != 3)    /* GT */
       {
        * Return_Score = LOW_SCORE;
        return  FALSE;
       }

   Sub = 0;
   for  (i = 0;  i < MARKOV_DEGREE;  i ++)
     Sub = ALPHABET_SIZE * Sub + S [i];

   Positive_Sum = Positive_Table [no] [MARKOV_DEGREE - 1] [0] [Sub];
   Negative_Sum = Negative_Table [no] [MARKOV_DEGREE - 1] [0] [Sub];

   for  (i = MARKOV_DEGREE;  i < DONOR_LEN;  i ++)
     {
      j = S [i];
      Positive_Sum += Positive_Table [no] [i] [j] [Sub];
      Negative_Sum += Negative_Table [no] [i] [j] [Sub];
      Sub = ALPHABET_SIZE * (Sub % (MARKOV_LEN / ALPHABET_SIZE)) + j;
     }
 
   Score = Positive_Sum - Negative_Sum;

#if  RETURN_TRUE_PROB
   X = exp (Positive_Sum + LOG_PROB_DONOR);
   Y = exp (Negative_Sum + LOG_PROB_NONDONOR);
   * Return_Score = log (X / (X + Y));
#else
   * Return_Score = Score;
#endif

	if(Score==-99) printf("look one\n");

   return(1);
  }

  
int  Is_Cod_NonCod  (const int * S, double * Return_Score, int ind)

/* Evaluate string  S [0 .. (CODING_LEN -1)] and
*  return  TRUE  or  FALSE  as to whether it is a likely donor
*  site.  Also set  Return_Score  to the probability that it is an donor
*  site. */

  {
   FILE  * Infile;
   static float  Positive_Table [4][CODING_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static float  Negative_Table [4][CODING_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static  int  Tables_Loaded[4] = {FALSE,FALSE,FALSE,FALSE};
   double  Positive_Sum, Negative_Sum, Score, Threshold;
   char filename[PATH_MAX+1];
   int no;


#if  RETURN_TRUE_PROB
   double  X, Y;
#endif
   int  i, j, k, Ct, Sub;

   no=ind;

   switch (no) {
   case 0: // case of exon in acceptor
     strcpy(filename,TRAIN_DIR);
     strcat(filename,"/score_ex.acc");
     break;
   case 1: // case of intron in acceptor
     strcpy(filename,TRAIN_DIR);
     strcat(filename,"/score_in.acc");
     break;
   case 2: // case of exon in donor
     strcpy(filename,TRAIN_DIR);
     strcat(filename,"/score_ex.don");
     break;
   case 3: // case of intron in donor
     strcpy(filename,TRAIN_DIR);
     strcat(filename,"/score_in.don");
     break;
   }

   if  (! Tables_Loaded[no] )
       {
        Infile = fopen (filename, "r");
        if  (Infile == NULL)
            {
             fprintf (stderr, "ERROR:  Unable to open donor file \"%s\"\n",
                        filename);
             exit (EXIT_FAILURE);
            }

        for  (i = MARKOV_DEGREE - 1;  i < CODING_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Positive_Table [no] [i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading donor file \"%s\"\n", 
                                DONOR_FILE_NAME);
                    exit (EXIT_FAILURE);
                   }
              }

        for  (i = MARKOV_DEGREE - 1;  i < CODING_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Negative_Table [no] [i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading donor file \"%s\"\n", 
                                DONOR_FILE_NAME);
                    exit (EXIT_FAILURE);
                   }
              }

        fclose (Infile);

        Tables_Loaded [no] = TRUE;
       }

   Sub = 0;
   for  (i = 0;  i < MARKOV_DEGREE;  i ++)
     Sub = ALPHABET_SIZE * Sub + S [i];

   Positive_Sum = Positive_Table [no] [MARKOV_DEGREE - 1] [0] [Sub];
   Negative_Sum = Negative_Table [no] [MARKOV_DEGREE - 1] [0] [Sub];

   for  (i = MARKOV_DEGREE;  i < CODING_LEN;  i ++)
     {
      j = S [i];
      Positive_Sum += Positive_Table [no] [i] [j] [Sub];
      Negative_Sum += Negative_Table [no] [i] [j] [Sub];
      Sub = ALPHABET_SIZE * (Sub % (MARKOV_LEN / ALPHABET_SIZE)) + j;
     }
 


   Score = Positive_Sum - Negative_Sum;

#if  RETURN_TRUE_PROB
   X = exp (Positive_Sum + LOG_PROB_DONOR);
   Y = exp (Negative_Sum + LOG_PROB_NONDONOR);
   * Return_Score = log (X / (X + Y));
#else
   * Return_Score = Score;
#endif

   return (1);
  }




