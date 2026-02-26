#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
/*******************************************************************/
/* HP-modellen på kvadratiskt gitter                               */
/* Kedjerepresentation:                                            */
/*       i=0,...,N-1            Index for kedjans pkter            */
/*       seq[i]                 +1 for H och 0 for P               */
/*       x[i],y[i]              Cartesiska koordinater             */ 
/*       0 <= x[i],y[i] < L                                        */
/*                                                                 */
/*       occ[x][y]              =i+1 om (x,y) ockuperad av kedje-  */
/*                              pkt i, =0 om (x,y) ej ockuperad    */
/*******************************************************************/                
#define NIT  100000000
#define NTHERM 10000000
#define IRT (NIT/100000)

long seed;
int N;
int L;
int N1;
int N2;
int *seq;
int *x;
int *y;
int *xf, *yf;
int **occ;
int dx[4];                                /* navigator             */
int dy[4];                                /* navigator             */
int nhh;                                  /* # HH-kontakter        */
double beta;
char rt[100];

void init(void);
int one_bead(double beta);
int two_bead(double beta);
int pivot(double beta);
int energy(int x[], int y[]); 
double ran3(long *ip);          

int main()
{
  extern int *x,*y,*seq,nhh,**occ;
  extern double beta;
  int i;
  long it;
  FILE *fp;
  
  init();
  
  fp=fopen(rt,"w");
  
  printf("Starting simulation...\n");
  
  for (it=1;it<=NIT;it++) {
    
    for (i=0;i<N1;i++) one_bead(beta);
    for (i=0;i<N2;i++) two_bead(beta);
    nhh=energy(x,y); 
    pivot(beta);
    
    if (it%IRT==0) fprintf(fp,"%li\t%d\n",it,-nhh);
    
  }
  
  printf("Simulation terminated successfully!\n");
  
  fclose(fp);
  
  for (i=0;i<L;i++) free(occ[i]);
  free(seq);
  free(x); free(y);
  free(xf); free(yf);
  free(occ);
  
  return 1;
}

void init(void) 
{
  extern int *x,*y,**occ,dx[4],dy[4],nhh,*seq;
  extern long seed;
  int i,j;
  double temp;
  char sigma[500];
  
  printf("HP model on square lattice \n");
  printf("-------------------------------------------\n");
  printf("Please enter sequence: ");
  scanf("%s", sigma);
  printf("Please enter temperature: ");
  scanf("%lf",&temp);
  
  N=strlen(sigma); N1=N-1; N2=N-2; L=2*N+3;
  
  seq=malloc(sizeof(int)*N);
  x=malloc(sizeof(int)*N);
  y=malloc(sizeof(int)*N);
  xf = malloc(sizeof(int)*N);
  yf = malloc(sizeof(int)*N);
  
  occ=malloc(sizeof(int*)*L);
  for (i=0;i<L;i++) occ[i]=malloc(sizeof(int)*L);
  
  sprintf(rt,"output_%4.2lf.dat",temp);
  
  /*** define navigators ***/ 
  
  dx[0]=1; dy[0]=0;
  dx[1]=0; dy[1]=1;
  dx[2]=-1; dy[2]=0;
  dx[3]=0; dy[3]=-1;
  
  for (i=0;i<L;i++) {for (j=0;j<L;j++) occ[i][j]=0;}
  
  for (i=0;i<N;i++) { 
    x[i]=N+1+i;
    y[i]=N;
    occ[x[i]][y[i]]=i+1;
  } 
  
  /*** read sequence ***/   
  
  for (j=0;j<N;j++) {
    if (sigma[j]=='H' || sigma[j]=='h') seq[j]=1;
    else if (sigma[j]=='P' || sigma[j]=='p') seq[j]=0;
    else {printf("Invalid sequence!\n"); exit(-1);}
  }
  nhh=0;
  
  beta=1./temp;
  
  seed=time(NULL);
  
}        
/**************************************************************************/
void neighbours(int x, int y, int xt[], int yt[])
{
  extern int dx[4], dy[4];
  int i;
  
  for (i=0; i<4; i++) {
    xt[i]=x+dx[i];
    yt[i]=y+dy[i];
  }
}
/**************************************************************************/
int one_bead(double beta)
{
  extern long seed;
  extern int *x,*y,dx[4],dy[4],**occ,*seq;
  int i,bead,nfree,xt[4],yt[4],xt2[4],yt2[4],lx[3],ly[3],
  xf,yf,deltaE,tmp,ax,ay,bx,by;
  
  bead=(N-1)*ran3(&seed)+1;
  
  if (bead==N-1) {
    /*** END-BEAD ***/
    /*** enumerate possible moves ***/
    nfree=0;
    neighbours(x[N-2],y[N-2],xt,yt);
    for (i=0;i<4;i++) {
      if (occ[(lx[nfree]=xt[i])][(ly[nfree]=yt[i])]==0) {
	nfree++;
      }
    }
    if (nfree==0) return 0;
    /*** select move ***/
    i=nfree*ran3(&seed); 
    xf=lx[i]; yf=ly[i];
    /*** energy change ***/
    deltaE=0;
    if (seq[N-1]==1) {
      neighbours(x[N-1],y[N-1],xt,yt);
      neighbours(xf,yf,xt2,yt2);
      for (i=0;i<4;i++) {
	tmp=occ[xt[i]][yt[i]]-1;
	if (tmp>=0 && tmp!=N-2) deltaE+=seq[tmp];
	tmp=occ[xt2[i]][yt2[i]]-1;
	if (tmp>=0 && tmp!=N-2) deltaE-=seq[tmp];
      }
    }
    /*** accept/reject ***/
    if (ran3(&seed)<exp(-beta*deltaE)) {
      occ[x[N-1]][y[N-1]]=0;
      x[N-1]=xf; y[N-1]=yf;
      occ[x[N-1]][y[N-1]]=N; 
      return 1;
    } 
    else
      return 0;
  }
  else {
    /*** INTERNAL BEAD ***/
    ax=(x[bead]-x[bead-1]); ay=(y[bead]-y[bead-1]);
    bx=(x[bead+1]-x[bead]); by=(y[bead+1]-y[bead]);
    if (ax*bx+ay*by!=0) return 0; 
    xf=(x[bead-1]+bx); yf=(y[bead-1]+by); 
    /*** check self-avoidance ***/
    if (occ[xf][yf]>0) return 0;
    /*** energy change ***/
    deltaE=0;
    if (seq[bead]==1) {
      neighbours(x[bead],y[bead],xt,yt);
      neighbours(xf,yf,xt2,yt2);
      for (i=0;i<4;i++) {
	tmp=occ[xt[i]][yt[i]]-1;
	if (tmp>=0 && tmp!=bead-1 && tmp!=bead+1) deltaE+=seq[tmp];
	tmp=occ[xt2[i]][yt2[i]]-1;
	if (tmp>=0 && tmp!=bead-1 && tmp!=bead+1) deltaE-=seq[tmp];
      }
    }
    /*** accept/reject ***/
    if (ran3(&seed)<exp(-beta*deltaE)) {
      occ[x[bead]][y[bead]]=0;
      x[bead]=xf; y[bead]=yf;
      occ[x[bead]][y[bead]]=bead+1;
      return 1;
    } 
    else
      return 0;
  }  
}
/**************************************************************************/
int two_bead(double beta)
{
  extern long seed;
  extern int *x,*y,dx[4],dy[4],**occ,*seq;
  int i,j,bead1,bead2,nfree,xt1[4],yt1[4],xt2[4],yt2[4],lx1[9],ly1[9],
  lx2[9],ly2[9],xf1,yf1,xf2,yf2,deltaE,d2,ax,ay,bx,by,cx,cy,tmp;
  
  bead1=(N-2)*ran3(&seed)+1; 
  bead2=bead1+1; 
  
  if (bead1==N-2) {
    /*** END-BEADS ***/
    /*** enumerate possible moves ***/
    nfree=0;
    neighbours(x[N-3],y[N-3],xt1,yt1);
    for (i=0;i<4;i++) {
      if (xt1[i]==x[N-4] && yt1[i]==y[N-4]) continue;
      if (occ[xt1[i]][yt1[i]]>0 && occ[xt1[i]][yt1[i]]<N-1) continue;
      neighbours(xt1[i],yt1[i],xt2,yt2);
      for (j=0;j<4;j++) {
	if (xt2[j]==x[N-3] && yt2[j]==y[N-3]) continue;
	if (occ[xt2[j]][yt2[j]]>0 && occ[xt2[j]][yt2[j]]<N-1) continue;
	if (xt1[i]==x[N-2] && yt1[i]==y[N-2] && 
	  xt2[j]==x[N-1] && yt2[j]==y[N-1])continue; 
	lx1[nfree]=xt1[i]; ly1[nfree]=yt1[i]; 
	lx2[nfree]=xt2[j]; ly2[nfree]=yt2[j]; 
	nfree++; 
      }
    }
    if (nfree==0) return 0;
    /*** select move ***/
    i=nfree*ran3(&seed); 
    xf1=lx1[i]; yf1=ly1[i];
    xf2=lx2[i]; yf2=ly2[i];
    /*** energy change ****/
    deltaE=0;
    if (seq[N-1]==1) {
      neighbours(x[N-1],y[N-1],xt1,yt1);
      neighbours(xf2,yf2,xt2,yt2);
      for (i=0;i<4;i++) {
	tmp=occ[xt1[i]][yt1[i]]-1;
	if (tmp>=0 && tmp!=N-2) deltaE+=seq[tmp];
	if (xt2[i]==xf1 && yt2[i]==yf1) continue;
	tmp=occ[xt2[i]][yt2[i]]-1;
	if (tmp>=0 && tmp<N-2) deltaE-=seq[tmp];
      }
    }
    if (seq[N-2]==1) {
      neighbours(x[N-2],y[N-2],xt1,yt1);
      neighbours(xf1,yf1,xt2,yt2);
      for (i=0;i<4;i++) {
	tmp=occ[xt1[i]][yt1[i]]-1;
	if (tmp>=0 && tmp!=N-3 && tmp!=N-1) deltaE+=seq[tmp];
	if (xt1[i]==xf2 && yt1[i]==yf2) continue;
	tmp=occ[xt2[i]][yt2[i]]-1;
	if (tmp>=0 && tmp<N-3) deltaE-=seq[tmp];
      }
    }
    /*** accept/reject ***/
    if (ran3(&seed)<exp(-beta*deltaE)) {
      occ[x[N-2]][y[N-2]]=0;
      occ[x[N-1]][y[N-1]]=0;
      x[N-2]=xf1; y[N-2]=yf1;
      x[N-1]=xf2; y[N-1]=yf2;
      occ[x[N-2]][y[N-2]]=N-1;
      occ[x[N-1]][y[N-1]]=N;
      return 1;
    } 
    else
      return 0;
  }
  else {
    /*** INTERNAL BEADS ***/
    /*** d2=1: crankshaft    d2=5: L-flip ***/
    d2=((x[bead2+1]-x[bead1-1])*(x[bead2+1]-x[bead1-1])+
    (y[bead2+1]-y[bead1-1])*(y[bead2+1]-y[bead1-1]))%L;
    ax=(x[bead1]-x[bead1-1]); ay=(y[bead1]-y[bead1-1]);
    cx=(x[bead2+1]-x[bead2]); cy=(y[bead2+1]-y[bead2]);
    
    if (d2!=1 && (ax*cx+ay*cy)!=0) return 0;
    if (d2==1) {
      xf1=(x[bead1-1]-ax); yf1=(y[bead1-1]-ay);
    }
    else {
      xf1=(x[bead1-1]+cx); yf1=(y[bead1-1]+cy);
    } 
    bx=(x[bead2]-x[bead1]); by=(y[bead2]-y[bead1]);
    if (xf1>=L) xf1-=L; if (xf1<0) xf1+=L;
    if (yf1>=L) yf1-=L; if (yf1<0) yf1+=L;
    xf2=(xf1+bx); yf2=(yf1+by);
    if (xf2>=L) xf2-=L; if (xf2<0) xf2+=L;
    if (yf2>=L) yf2-=L; if (yf2<0) yf2+=L;
    /*** check self-avoidance ***/ 
    if (occ[xf1][yf1]>0 || occ[xf2][yf2]>0) return 0;
    /*** energy change ***/
    deltaE=0;
    if (seq[bead1]==1) {
      neighbours(x[bead1],y[bead1],xt1,yt1);
      neighbours(xf1,yf1,xt2,yt2);
      for (i=0;i<4;i++) {
	tmp=occ[xt1[i]][yt1[i]]-1;
	if (tmp>=0 && tmp!=bead1-1 && tmp!=bead2) deltaE+=seq[tmp];
	tmp=occ[xt2[i]][yt2[i]]-1;
	if (tmp>=0 && tmp!=bead1-1 && tmp!=bead2) deltaE-=seq[tmp];
      }
    }
    if (seq[bead2]==1) {
      neighbours(x[bead2],y[bead2],xt1,yt1);
      neighbours(xf2,yf2,xt2,yt2);
      for (i=0;i<4;i++) {
	tmp=occ[xt1[i]][yt1[i]]-1;
	if (tmp>=0 && tmp!=bead1 && tmp!=bead2+1) deltaE+=seq[tmp];
	tmp=occ[xt2[i]][yt2[i]]-1;
	if (tmp>=0 && tmp!=bead2+1 && tmp!=bead1) deltaE-=seq[tmp];
      }
    }
    /*** accept/reject ***/
    if (ran3(&seed)<exp(-beta*deltaE)) {
      occ[x[bead1]][y[bead1]]=0;
      occ[x[bead2]][y[bead2]]=0;
      x[bead1]=xf1; y[bead1]=yf1;
      x[bead2]=xf2; y[bead2]=yf2;
      occ[x[bead1]][y[bead1]]=bead1+1;
      occ[x[bead2]][y[bead2]]=bead2+1;
      return 1;
    } 
    else
      return 0;
  }
}    
/**************************************************************************/
int pivot(double beta)
{
  extern long seed;
  extern int *x,*y,**occ,*seq,nhh;
  extern int *xf, *yf;
  int i,j,k,symm_op,tmp,nhhf;
  
  
  /*** generate trial configuration ***/
  
  k=(N-2)*ran3(&seed)+1;                        /* pivotpkt               */
  
  symm_op=7*ran3(&seed);                        /* typ av forandring      */
  
  for (i=k+1;i<N;i++) {                         /* koord. rel. pivotpkten */
    xf[i]=x[i]-x[k]; yf[i]=y[i]-y[k];           
  }
  
  switch (symm_op) {
    case 0:                                       /* rot. +90               */
      for (i=k+1;i<N;i++) {
	tmp=xf[i];
	xf[i]=-yf[i];
	yf[i]=tmp;
      }
      break;
    case 1:                                       /* rot. -90               */
      for (i=k+1;i<N;i++) {
	tmp=xf[i];
	xf[i]=yf[i];
	yf[i]=-tmp;
      }
      break;
    case 2:                                       /* rot. 180               */
      for (i=k+1;i<N;i++) {    
	xf[i]=-xf[i];
	yf[i]=-yf[i];
      }
      break;
    case 3:                                       /* spegling (1,0)         */ 
      for (i=k+1;i<N;i++) {
	xf[i]=xf[i];
	yf[i]=-yf[i];
      }
      break;
    case 4:                                       /* spegling (0,1)         */
      for (i=k+1;i<N;i++) {
	xf[i]=-xf[i];
	yf[i]=yf[i];
      }
      break;
    case 5:                                       /* spegling (1,1)         */
      for (i=k+1;i<N;i++) {
	tmp=xf[i];
	xf[i]=yf[i];
	yf[i]=tmp;
      }
      break;
    case 6:                                       /* spegling (1,-1)        */
      for (i=k+1;i<N;i++) {
	tmp=xf[i];
	xf[i]=-yf[i];
	yf[i]=-tmp;
      }
      break;
  }
  
  for (i=k+1;i<N;i++) {xf[i]+=x[k]; yf[i]+=y[k];}/* satt samman          */
    for (i=0;i<=k;i++) {xf[i]=x[i]; yf[i]=y[i];}
    
    /*** nya konfigurationen klar --- OK? ***/
    
    for (i=k+1;i<N;i++) occ[x[i]][y[i]]=0;
    
    for (i=k+1;i<N;i++) {
      if (occ[xf[i]][yf[i]]>0) {
        for (j=k+1;j<i;j++) occ[xf[j]][yf[j]]=0;
        for (j=k+1;j<N;j++) occ[x[j]][y[j]]=j+1;
        return 0;
      }
      else
	occ[xf[i]][yf[i]]=i+1; 
    }
    
    /*** acceptera? ***/
    
    nhhf=energy(xf,yf);
    //printf("vcvlou %d %d %lf \n",nhh, nhhf, beta);
    
    if (ran3(&seed)<exp(-beta*(-nhhf+nhh))) {
      for (i=k+1;i<N;i++) { 
        x[i]=xf[i];
        y[i]=yf[i];
      }
      
      nhh=nhhf;
      return 1;
    }
    else {
      for (i=k+1;i<N;i++) occ[xf[i]][yf[i]]=0;
      for (i=k+1;i<N;i++) occ[x[i]][y[i]]=i+1;
      
      return 0;
    }
    
}
/*******************************************************************/
int energy(int x[],int y[]) 
{
  extern int **occ, dx[4], dy[4];
  int i,j,nhh, xt, yt;
  
  nhh=0;
  
  for (i=0; i<N; i++) {
    if (seq[i]!=1) continue;
    for (j=0; j<2; j++) {
      xt=x[i]+dx[j]; yt=y[i]+dy[j];
      if (occ[xt][yt]-1>=i-1 && occ[xt][yt]-1<=i+1) continue;
      if (seq[occ[xt][yt]-1]==1) nhh++;
    }
  }
  
  return nhh;
}
/*******************************************************************/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(long *seed)
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  double ret_val;
  
  if (*seed < 0 || iff == 0) {
    iff=1;
    mj=MSEED-(*seed < 0 ? -*seed : *seed);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
	ma[i] -= ma[1+(i+30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
      inext=0;
    inextp=31;
    *seed=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  ret_val = mj*FAC;
  if (mj == 0) ret_val = FAC;
  return ret_val;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
/******************************************************************/
