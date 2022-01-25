#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define MAXN 4097
#define INF 9999

clock_t start_s, stop_s;

long long int count_fw=0, count_fwi=0,count_nfw=0;
int A[MAXN][MAXN], Ai[MAXN][MAXN], An[MAXN][MAXN];

void fw(int N)
{
    int i,j,k;

    printf("Floyd-Warshall Algorithm\n");
    start_s=clock();
       for (k=0;k<N;k++)
         for (i=0;i<N;i++)
            for (j=0;j<N;j++)
            {
                count_fw++;
                if ((A[i][k]<INF) && (A[k][j]<INF)&& ((A[i][k]+A[k][j])<A[i][j]))
                   A[i][j]=A[i][k]+A[k][j];
            }
    stop_s=clock();
    printf("\nExecution time: %f\n",(double)(stop_s - start_s));
}

void fwi(int N)
{
    int i,j,k;

    printf("Floyd-Warshall Algorithm improved\n");
    start_s=clock();
    for (k=0;k<N;k++)
         for (i=0;i<N;i++)
            if (Ai[i][k]<INF)
               for (j=0;j<N;j++)
                  if (Ai[k][j]<INF)
                 {
                     count_fwi++;
                     if ((Ai[i][k]<INF) && (Ai[k][j]<INF)&& ((Ai[i][k]+Ai[k][j])<Ai[i][j]))
                        Ai[i][j]=Ai[i][k]+Ai[k][j];
                  }
    stop_s=clock();
    printf("\nExecution time: %f\n",(double)(stop_s - start_s));
}

void new_fw(int N)
{
    static int inc[MAXN],outc[MAXN];
    static int inlist[MAXN][MAXN], outlist[MAXN][MAXN];
    int i,j,k,kk;
    static int select_k[MAXN],mininxout,mink;

    printf("New Floyd Warshall Algorithm\n");

    for (i=0;i<N;i++)
        inc[i]=0,outc[i]=0,select_k[i]=0;
    start_s=clock();
    for (i=0;i<N;i++)
       for (j=0;j<N;j++)
       {
         count_nfw++;  // initial calculations 
         if ((An[i][j]!=0) && (An[i][j]<INF))
         { inc[j]++;outc[i]++;inlist[j][inc[j]-1]=i;outlist[i][outc[i]-1]=j;}
       }


    for (kk=0;kk<N;kk++)
    {
        mink=-1;
        mininxout=2*N*N;
        for (k=0;k<N;k++)
        {
           count_nfw++;  // for each k: find min
           if ((select_k[k]==0) &&(inc[k]*outc[k]<mininxout))
           {
              mink=k;
              mininxout=inc[k]*outc[k];
           }
       }
       k=mink;
       select_k[k]=1;
       for (i=0;i<inc[k];i++)
            for (j=0;j<outc[k];j++)
            {
                count_nfw++;
                if ((An[inlist[k][i]][k]+An[k][outlist[k][j]])<An[inlist[k][i]][outlist[k][j]])
                {
                   if (An[inlist[k][i]][outlist[k][j]]==INF) //NO NEGATIVE CYCLE ASSUMED
                   {
                      outc[inlist[k][i]]++; outlist[inlist[k][i]][outc[inlist[k][i]]-1]=outlist[k][j];
                      inc[outlist[k][j]]++; inlist[outlist[k][j]][inc[outlist[k][j]]-1]=inlist[k][i];
                   }
                   An[inlist[k][i]][outlist[k][j]]=An[inlist[k][i]][k]+An[k][outlist[k][j]];
               }
           }
    }
    stop_s=clock();
    printf("\nExecution time: %f\n",(double)(stop_s - start_s));
}


int main()
{
    int N, i, j;

    /* input from user */
    /***
    scanf("%d",&N);
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
           scanf("%d",&A[i][j]),An[i][j]=Ai[i][j]=A[i][j];
    ***/
    

     /* random input */

    N=1024; 
    int lgN=log(N)/log(2);
    printf("Clocks per sec %d\n",CLOCKS_PER_SEC);
    int TESTS[]={1,2,4,8,16,4*lgN,8*lgN,4*N/lgN,2*N};
    for (int x=0;x<9;x++) for (int y=0;y<5;y++)
    {
       count_fw= count_fwi= count_nfw=0; // correct counts

       for (i=0;i<N;i++)
          for (j=0;j<N;j++)
             if (i==j) A[i][j]=Ai[i][j]=An[i][j]=0;
             else A[i][j]=Ai[i][j]=An[i][j]=INF;

       int c=0;
       for (i=0;i<N;i++)
           for (j=0;j<N;j++)
               if ((i!=j)&&(rand()%4*N>=(4*N-TESTS[x]))) A[i][j]=Ai[i][j]=An[i][j]=1,c++;
       printf("%d\n",c);

       fwi(N);

       new_fw(N);

       for (i=0;i<N;i++)
           for (j=0;j<N;j++)
               if (An[i][j] != Ai[i][j])
                   printf("%d %d\n",i,j);

       printf("FW %ld FWi %ld NFW %ld\n",count_fw,count_fwi,count_nfw);
    }
    return 0;
}
/*
5
0 6 9999 5 9999
2 0 3 -1 2
-2 9999 0 2 9999
-1 1 2 0 -1
1 9999 9999 9999 0
*/
