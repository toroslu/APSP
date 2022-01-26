#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define MAXN 4097
#define INF 9999

clock_t start_s, stop_s;

int G[MAXN][MAXN] , distance[MAXN] ,  parent[MAXN], h[MAXN] ;
int APSPD[MAXN][MAXN], APSPD2[MAXN][MAXN],APSPJ[MAXN][MAXN];

struct neighbour { int v; int w; struct neighbour* n; };
struct vertex { int odegree; struct neighbour *first, *last;};
struct vertex Graph [MAXN];

long long int count_dijkstra=0, count_bf=0, count_johnson=0, count_fw=0, count_fwi;

// don't use index 0
// heap values, positions in the position array, position of keys of index 1 to N in heap
int heap_val[MAXN], heap_pos[MAXN], position[MAXN];

struct min_vertex {int val, id;};

void swap(int *a, int *b){int t=*a; *a=*b; *b=t;}

// decrease the value at postion i to a new value

void decreasekey_heap(int i, int new_value)
{
    int k1;
    int k=position[i]; // position in a heap
    heap_val[k]=new_value;
    while (k>1)
    {
        count_dijkstra++;
        count_johnson++;
        k1=k / 2; //parent
        if (heap_val[k1]>heap_val[k])
        {
            swap(&position[heap_pos[k1]],&position[heap_pos[k]]);
            swap(&heap_val[k1],&heap_val[k]);
            swap(&heap_pos[k1],&heap_pos[k]);
        }
        else break;
        k=k1;
    }
}

void addkey_heap(int i, int new_value, int *heap_size)
{
    (*heap_size)++;
    int k=*heap_size;
    position[i]=k;
    heap_val[k]=new_value;
    heap_pos[k]=i;
    decreasekey_heap(i,new_value);
}

struct min_vertex deletemin_heap(int *heap_size)
{
    int k1, k2, k=1;
    struct min_vertex r;
    r.val=heap_val[k];
    r.id=heap_pos[k];
    heap_val[k]=heap_val[*heap_size];
    heap_pos[k]=heap_pos[*heap_size];
    position[heap_pos[k]]=k;


    (*heap_size)--;
    while (k<*heap_size)
    {
        count_dijkstra++;
        count_johnson++;
        k1=2*k; k2=2*k+1; //children
        if (k1>*heap_size) break;
        if (k2<=*heap_size)
            if (heap_val[k2]<heap_val[k1]) k1=k2;
        if (heap_val[k1]<heap_val[k])
        {
            swap(&position[heap_pos[k1]],&position[heap_pos[k]]);
            swap(&heap_val[k1],&heap_val[k]);
            swap(&heap_pos[k1],&heap_pos[k]);
        }
        else break;
        k=k1;
    }
    return r;
}

void apspd( int V, int S)
{
    int i ;
    for ( i = 1 ; i <= V ; i++ )
        APSPD[S][i]=distance[i];
}

void dijkstra ( int S , int V )
{
    int u;
    struct min_vertex r;
    int N=V;
    while( N >0 )
    {
        r=deletemin_heap(&N);
        u=r.id;
        distance[u]=r.val;
        for (struct neighbour *t=Graph[u].first; t!=NULL;t=t->n)
        {
                if ((distance[u]!=INF)&&( distance[u] + t->w < distance[t->v]))
                    decreasekey_heap(t->v,distance[u] + t->w ), distance[t->v]=distance[u]+t->w, parent[t->v] = u ;
        }
    }
    apspd(V,S);
}

void  main_Dijkstra(int V, int S)
{
    int i;
    int N=0;

    for ( i = 1 ; i <= V ; i++ )
        addkey_heap(i,INF,&N), distance[i] = INF , parent[i] = -1 ;

    decreasekey_heap(S,0);
    distance[S] = 0 ;
    dijkstra (S, V) ;

}

int Bellman_Ford0(int V,  int S)
{
    int i,u,flag=1; //distance[MAXN],parent[MAXN],flag=1;
    for(i=0;i<=V;i++)
        distance[i] = INF , parent[i] = -1 ;

    distance[S]=0 ;
    for(i=0;i<=V-1;i++)
    {
        for (u=0;u<=V;u++) //,printf("%d\n",count_johnson))
            for (struct neighbour *t=Graph[u].first; t!=NULL;t=t->n)
            {
                count_johnson++;
                count_bf++;
                if((distance[u]!=INF)&&( distance[u] + t->w < distance[t->v]))
                   distance[t->v]=distance[u]+t->w, parent[t->v] = u ;
            }
    }

        for (u=0;u<=V;u++)
            for (struct neighbour *t=Graph[u].first; t!=NULL;t=t->n)
            {
                count_johnson++;
                count_bf++;
                if((distance[u]!=INF)&&( distance[u] + t->w < distance[t->v]))
                   flag=0;
            }

        return flag;
}
void main_Bellman_Ford0(int V)  // add vertex 0, find SP from vertex 0
{
    struct neighbour *t;
    int i;
    Graph[0].odegree=V;
    Graph[0].first=(struct neighbour *)malloc(sizeof (struct neighbour));
    for ( i=1,t=Graph[0].first;i<=V;i++)
    {
        t->v=i;
        t->w=0;
        if (i<V) {t->n=(struct neighbour *)malloc(sizeof (struct neighbour)); t=t->n;}
        else t=NULL;
    }

    if(Bellman_Ford0(V,0))
        ;  // No negative weight cycle
    else printf("\nNegative weight cycle exists\n");
}

int A[MAXN][MAXN], Ai[MAXN][MAXN];

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
    printf("Execution time: %f\n",(double)(stop_s - start_s));
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
    printf("Execution time: %f\n",(double)(stop_s - start_s));
}


void recalculate_weights(int V)
{
    int u;
    for (u=0;u<=V;u++)
        for (struct neighbour *t=Graph[u].first; t!=NULL;t=t->n)
        {
            count_johnson++;
            t->w=t->w+h[u]-h[t->v];
        }
}

void johnson(int V)
{
    int i,j,S;
    printf("Johnson's Algorithm\n");
    main_Bellman_Ford0(V);

    for (i=0;i<=V;i++)
    {
        h[i]=distance[i];
    }
    recalculate_weights(V);
    for (S=1;S<=V;S++)
        main_Dijkstra(V, S);

    //    APSP Johnson
    for (i=1;i<=V;i++)
        for (j=1;j<=V;j++)
            if (APSPD[i][j] !=INF) APSPJ[i][j]=APSPD[i][j]+h[j]-h[i];
            else APSPJ[i][j]=APSPD[i][j];
}

int main()
{
    int N, i, j, S;

    /* input from user */
    /***
    scanf("%d",&N);
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
           scanf("%d",&A[i][j]),G[i][j]=Ai[i][j]=A[i][j];
    ***/
    

     /* random input */
     /***/
    N=1024; 
    int lgN=log(N)/log(2);
    printf("Clocks per sec %d\n",CLOCKS_PER_SEC);
    int TESTS[]={1,2,4,8,16,4*lgN,8*lgN,4*N/lgN,2*N};
    for (int x=0;x<9;x++) for (int y=0;y<5;y++)
    {
       count_fw= count_fwi= count_dijkstra=count_dijkstra=count_bf=0; // correct counts
       for (i=0;i<N;i++)
          for (j=0;j<N;j++)
              if (i==j) G[i][j]=Ai[i][j]=A[i][j]=0;
              else G[i][j]=Ai[i][j]=A[i][j]=INF;
       int c=0;
       for (i=0;i<N;i++)
           for (j=0;j<N;j++)
               if ((i!=j)&&(rand()%(4*N)>=(4*N-TESTS[x]))) A[i][j]=Ai[i][j]=G[i][j]=1,c++;
       printf("\nNumber of Edges: %d\n",c);
       /***/
       fwi(N);

       for ( i=1 ; i<=N ; i++ )
       {
          Graph[i].odegree=0;
          Graph[i].first=NULL;
          Graph[i].last=NULL;
          for( j = 1 ; j <= N ; j++ )
          {
             if ((G[i-1][j-1]!=INF) && (G[i-1][j-1] != 0)) 
             {   // to be able to handle negative wights
                 Graph[i].odegree++;
                 if (Graph[i].odegree==1)
                 {
                     Graph[i].first=(struct neighbour *)malloc(sizeof (struct neighbour));
                     Graph[i].last=Graph[i].first;
                 }
                 else
                 {
                     Graph[i].last->n=(struct neighbour *)malloc(sizeof (struct neighbour));
                     Graph[i].last=Graph[i].last->n;
                 }
                 Graph[i].last->v=j;
                 Graph[i].last->w=G[i-1][j-1];
                 Graph[i].last->n=NULL;
              }
          }
      }

    start_s=clock();
    johnson(N);  
    stop_s=clock();
    printf("Johnson's Execution time: %f\n",(double)(stop_s - start_s));

    for (i=1;i<=N;i++)
        for (j=1;j<=N;j++)
            if (Ai[i-1][j-1] != APSPJ[i][j])
              printf("%d %d\n",i,j);
   
    printf("BF %ld Dij %ld John %ld FWi %ld \n",count_bf,count_dijkstra,count_johnson,count_fwi);
  }
  return 0;
}