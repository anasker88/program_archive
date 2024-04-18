/********************************************************************
* ファイル名  ：mat3d.c                                             *
* タイトル    ：線形システム解法プログラム                          *
* 製作者      ：平野　博之                                          *
* 所属        ：岡山理科大学 工学部 応用化学科                      *
* 製作日      ：2003.12.25                                          *
* 言語        ：C                                                   *
*********************************************************************
*  注意：Ｃ言語では配列の添字が０から始まるので、FORTRAN のB(NE)を  *
*  そのままb[NE]と宣言すると、b[0]からb[NE-1]を意味する．ここでは   *
*  できる限りFORTRANプログラムに近くなるように、b[NE+1]として添字   *
*  を1からNEまで変化させてある．                                    *
********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NX 5/* NX:x方向格子数 */
#define NY 5/* NY:y方向格子数 */
#define NZ 5/* NZ:z方向格子数 */
#define NE 125/* NE:全格子数=NX*NY */
#define NXY 25 /* NXY=NX*NY */
#define NX1 6 /* NX1=NX+1 */
#define NY1 6 /* NY1=NY+1 */
#define NZ1 6 /* NZ1=NZ+1 */
#define NE1 126 /* NE1=NE+1 */
#define NXY2 51 /* NXY2=NXY*2+1 */

void dlu (double [NE1][NE1], double [NE1], double [NE1], int ne);
void dfgp (double [NE1][NE1],double [NE1],double [NE1],int ne);
void dfg (double [NE1][NE1], double [NE1], double [NE1], int ne);
void dgjp (double [NE1][NE1],double [NE1],double [NE1],int ne);
void dgbnd (double [NE1], double [NE1], double [NE1],
            double [NE1], double [NE1], double [NE1], double [NE1],
            double [NE1], double [NE1], int ny, int ne, int nxy);
void hjacob (double [NE1][NE1], double [NE1], double [NE1], int ne);
void hsor (double [NE1][NE1], double [NE1], double [NE1], int ne);
void sorbnd (double [NE1], double [NE1], double [NE1], double [NE1], 
             double [NE1], double [NE1], double [NE1], double [NE1], double [NE1],
             int ny, int nxy, int ne);
void lblbnd (double [NE1] ,double [NE1], double [NE1],
             double [NE1], double [NE1], double [NE1], double [NE1],
             double [NE1], double [NE1],
             int nx, int ny, int nz, int nxy, int ne);
void kcr (double [NE1][NE1], double [NE1], double [NE1], int ne);
void profmv (double [NE1][NE1], double [NE1], double [NE1] ,int ne);
double provv (double [NE1], double [NE1], int ne);
void kcrbnd (double [NE1], double [NE1], double [NE1], double [NE1],
             double [NE1], double [NE1], double [NE1],
             double [NE1], double [NE1],
             int ny, int nxy, int ne);
void probmv (double [NE1], double [NE1], double [NE1], double [NE1],
             double [NE1], double [NE1], double [NE1],
             double [NE1], double [NE1],
             int ny, int nxy, int ne);
void kbicg (double [NE1][NE1], double [NE1], double [NE1], int ne);
void kbibnd (double [NE1], double [NE1], double [NE1], double [NE1], double [NE1],
             double [NE1], double [NE1], double [NE1], double [NE1],
             int ny, int nxy, int ne);
void main()
{
     double a[NE1][NE1], b[NE1], x1[NE1], x3[NX1][NY1][NZ1];
/* バンドマトリックス用配列の定義 */
     double at1[NE1],at2[NE1],at3[NE1],at4[NE1],at5[NE1],at6[NE1],at7[NE1];
/* 上記離散化方程式中のΔ=1を設定 */
     double delta=1.0;
     int i,ii,ix,iy,iz,nx,ny,nz,ne,nxy;
     int ij;
     double b1,b2,b3,b4,b5,b6,b7;
     char s[250];
     nx=NX; ny=NY; nz=NZ; ne=NE; nxy=NXY;
     next_cal:{};/* 解法の変更毎の戻り点 */
/* A,Bを設定するためにあらかじめx2の値を格子番号と同じになるようにする */
     i = 1;
     for (iz=1;iz<=nz;iz++){
       for (ix=1;ix<=nx;ix++){
         for (iy=1;iy<=ny;iy++){
           x3[ix][iy][iz]=i; i=i+1;
         }
       }
     }
     for (i=1;i<=ne;i++){/* X3以外の配列のゼロクリア */
       at1[i]=0.0;at2[i]=0.0;at3[i]=0.0;at4[i]=0.0;at5[i]=0.0;at6[i]=0.0;at7[i]=0.0;
       b[i]=0.0; x1[i]=0.0;
       for (ii=1;ii<=ne;ii++){
         a[i][ii] = 0.0;
       }
     }
/* X_{i,j,k}=X_{m=(k-1)*(NX*NY)+(i-1)*NY+j}=mとなるようにマトリックスA,Bを求める */
     i=1;
     for (iz=1;iz<=NZ;iz++){
       for (ix=1;ix<=nx;ix++){
         for (iy=1;iy<=ny;iy++){
         /*f(i,j,k)は常に計算領域内*/
         at3[i] = -2.0/(delta*delta) - 2.0/(delta*delta) - 2.0/(delta*delta);
         a[i][i] = at3[i];
         b3 = x3[ix][iy][iz]*(-2.0/(delta*delta)-2.0/(delta*delta)-2.0/(delta*delta));
         /* f(i-1,j,k)が計算領域外（左側）の境界条件の処理 */
         if (ix==1){
           at1[i]=0.0;/* f(0,j,k)=0として処理する */
           b1=0.0;/* AT1*f(0,j,k)をB1とする */
         }
         /* f(i-1,j,k)が計算領域内の場合 */
         else {
           at1[i] = 1.0/(delta*delta) - 1.0/(2.0*delta);
           a[i][i-ny] = at1[i];
           b1 = x3[ix-1][iy][iz] * ( 1.0/(delta*delta) - 1.0/(2.0*delta) );
         }
         /* f(i+1,j,k)が計算領域外（右側）の境界条件の処理 */
         if (ix==nx){
           at5[i] = 0.0;/* f(NX,j,k)=0として処理 */
           b5 = 0.0;/* AT5*f(NX,j,k)=B5とする */
         }
         /* f(i+1,j,k)が計算領域内の場合 */
         else {
           at5[i] = 1.0/(delta*delta) + 1.0/(2.0*delta);
           a[i][i+ny] = at5[i];
           b5 = x3[ix+1][iy][iz]*( 1.0/(delta*delta) + 1.0/(2.0*delta) );
         }
         /* f(i,j-1,k)が計算領域外（下側）の境界条件の処理 */
         if (iy==1){
           at2[i] = 0.0;/* f(i,0,k)=0として処理 */
           b2 = 0.0;/* AT2*f(i,0,k)=B2とする */
         }
         /* f(i,j-1,k)が計算領域内の場合 */
         else {
           at2[i] = 1.0/(delta*delta) - 1.0/(2.0*delta);
           a[i][i-1] = at2[i];
           b2 = x3[ix][iy-1][iz]*( 1.0/(delta*delta) - 1.0/(2.0*delta) );
         }
         /* f(i,j+1,k)が計算領域外（上側）の境界条件の処理 */
         if (iy==ny){
           at4[i] = 0.0;/* f(i,NY,k)=0として処理 */
           b4 = 0.0;/* AT4*f(i,NY,k)=B4とする */
         }
         /* f(i,j+1,k)が計算領域内の場合 */
         else {
           at4[i] = 1.0/(delta*delta) + 1.0/(2.0*delta);
           a[i][i+1] = at4[i];
           b4 = x3[ix][iy+1][iz]*( 1.0/(delta*delta) + 1.0/(2.0*delta) );
         }
         /* f(i,j,k-1)が計算領域外（後面）の境界条件の処理 */
         if (iz==1){
           at6[i] = 0.0;/* f(i,j,0)=0として処理 */
           b6 = 0.0;/* AT6*f(i,j,0)=B6とする */
         }
         /* f(i,j,k-1)が計算領域内の場合 */
         else{
           at6[i] = 1.0/(delta*delta) - 1.0/(2.0*delta);
           a[i][i-nxy] = at6[i];
           b6=x3[ix][iy][iz-1]*( 1.0/(delta*delta)-1.0/(2.0*delta) );
         }
         /* f(i,j,k+1)が計算領域外（前面）の境界条件の処理 */
         if (iz==nz){
           at7[i] = 0.0;/* f(i,j,NZ)=0として処理 */
           b7 = 0.0;/* AT7*f(i,j,NZ)=B7とする */
         }
         /* f(i,j,k+1)が計算領域内の場合 */
         else{
           at7[i] = 1.0/(delta*delta) + 1.0/(2.0*delta);
           a[i][i+nxy] = at7[i];
           b7=x3[ix][iy][iz+1]*( 1.0/(delta*delta) + 1.0/(2.0*delta) );
         }
         /* Bの計算 */
         b[i] = b1+b2+b3+b4+b5+b6+b7; i=i+1;
         }
       }
     }
     printf("  1 : LU Decomposition \n"); /* 各サブルーチンを用いた計算 */
     printf("  2 : Gauss Elimination \n");
     printf("  3 : Gauss Elimination without Pivotting \n");
     printf("  4 : Gauss Jordan \n");
     printf("  5 : Band Gauss \n");
     printf("  6 : Jacobi \n");
     printf("  7 : point-SOR \n");
     printf("  8 : point-SOR for Band Matrix \n");
     printf("  9 : line-SOR for Band Matrix \n");
     printf(" 10 : Conjugate Residual \n");
     printf(" 11 : Conjugate Residual for Band Matrix \n");
     printf(" 12 : Bi-CGSTAB \n");
     printf(" 13 : Bi-CGSTAB for Band Matrix \n");
     printf(" 14 : END \n");
     printf(" Input Number ---> ");
     gets(s); ij=atoi(s);
     switch(ij){
       case  1 : dlu (a,b,x1,ne); break;
       case  2 : dfgp (a,b,x1,ne); break;
       case  3 : dfg (a,b,x1,ne); break;
       case  4 : dgjp (a,b,x1,ne); break;
       case  5 : dgbnd (at1,at2,at3,at4,at5,at6,at7,b,x1,ny,nxy,ne); break;
       case  6 : hjacob (a,b,x1,ne); break;
       case  7 : hsor (a,b,x1,ne); break;
       case  8 : sorbnd (at1,at2,at3,at4,at5,at6,at7,b,x1,ny,nxy,ne); break;
       case  9 : lblbnd (at1,at2,at3,at4,at5,at6,at7,b,x1,nx,ny,nz,nxy,ne); break;
       case 10 : kcr (a,b,x1,ne); break;
       case 11 : kcrbnd (at1,at2,at3,at4,at5,at6,at7,b,x1,ny,nxy,ne); break;
       case 12 : kbicg (a,b,x1,ne); break;
       case 13 : kbibnd (at1,at2,at3,at4,at5,at6,at7,b,x1,ny,nxy,ne); break;
     }
     if (ij != 14){
       goto next_cal;
     }
}
void dlu (double a[NE1][NE1], double b[NE1], double x[NE1], int ne)
{
     double y[NE1];
     int ipv[NE1];
     int i,j,k,l,ipvex;
     double eps,apv,t;
     for (i=1;i<=ne;i++){/* 配列のゼロクリア */
       y[i] = 0.0;
     }
     for (i=1;i<=ne;i++){/* 部分軸選択用配列の初期設定 */
       ipv[i] = i;
     }
     eps = 1.0e-50;/* 部分軸選択による特異性の判定値 */
     for (k=1;k<=ne;k++){/* 部分軸選択 */
       l = k; apv = fabs( a[ipv[l]][k] );
       for (i=k+1;i<=ne;i++){
         if ( fabs( a[ipv[i]][k] ) > apv ){/* 部分選択 */
           l = i; apv = fabs(a[ipv[l]][k]);
         }
       }
       if (l != k ){/* 部分軸選択を行った方がよいと判断したら入れ替えを行う */
         ipvex = ipv[k]; ipv[k] = ipv[l]; ipv[l] = ipvex;
       }
       if ( fabs( a[ipv[k]][k] )<=eps ){/*部分軸選択を行っても計算不可能なとき:行列は特異(singular)*/
         printf("  Matrix is singular at k= %d ",k);
         return;
       }
       a[ipv[k]][k] = 1.0 / a[ipv[k]][k];/* U を求める */
       for (i=k+1;i<=ne;i++){
         a[ipv[i]][k] = a[ipv[i]][k] * a[ipv[k]][k];
         for (j=k+1;j<=ne;j++){
           a[ipv[i]][j] = a[ipv[i]][j] - a[ipv[i]][k]*a[ipv[k]][j];
         }
       }
     }
     y[1] = b[ipv[1]];/* 第1段階 LY = B を解く */
     for (i=2;i<=ne;i++){
       t = b[ipv[i]];
       for (j=1;j<=i-1;j++){
         t = t - a[ipv[i]][j] * y[j];
       }
       y[i] = t;
     }
     x[ne] = y[ne] * a[ipv[ne]][ne];/* 第2段階 UX = Y を解く */
     for (i=ne-1;i>=1;i--){
       t = y[i];
       for (j=i+1;j<=ne;j++){
         t = t - a[ipv[i]][j] * x[j];
       }
       x[i] = t * a[ipv[i]][i];
     }
     for (i=1;i<=ne;i++){/* 結果出力 */
       printf(" X_1_( %d ) = %15.6lE \n", i, x[i]);
     }
}
void dfgp (double a[NE1][NE1],double b[NE1],double x[NE1],int ne)
{
     int ipv[NE1];
     int i,j,k,l,ipvex;
     double eps,apv,t;
     eps = 1.0e-50;/* 部分軸選択による特異性の判定値 */
     for(i=1;i<=ne;i++){/* 部分軸選択用配列の初期設定 */
       ipv[i] = i;
     }
     for (k=1;k<=ne;k++){
       l = k; apv = fabs( a[ipv[l]][k] );
       for (i=k+1;i<=ne;i++){/* 部分軸選択 */
         if ( fabs( a[ipv[i]][k] ) > apv ){
           l=i; apv = fabs( a[ipv[l]][k] );
         }
       }
       if ( l != k ){/* 部分軸選択を行った方がよいと判断したら入れ替えを行う */
         ipvex = ipv[k]; ipv[k] = ipv[l]; ipv[l] = ipvex;
       }
       if ( fabs( a[ipv[k]][k] )<=eps ){/*部分軸選択を行っても計算不可能なとき:行列は特異*/
         printf(" Matrix is singular at k= %d ",k );
         return;
       }
       for (i=k+1;i<=ne;i++){/* 前進消去 */
         a[ipv[k]][i] = a[ipv[k]][i]/a[ipv[k]][k];
       }
       b[ipv[k]] = b[ipv[k]]/a[ipv[k]][k];
       for (i=k+1;i<=ne;i++){
         for (j=k+1;j<=ne;j++){
           a[ipv[i]][j] = a[ipv[i]][j] - a[ipv[i]][k]*a[ipv[k]][j];
         }
         b[ipv[i]] = b[ipv[i]]-a[ipv[i]][k]*b[ipv[k]];
       }
     }
     x[ne] = b[ipv[ne]];/* 後退代入*/
     for (i=ne-1;i>=1;i--){
       t = b[ipv[i]];
       for (j=i+1;j<=ne;j++){
         t = t - a[ipv[i]][j] * x[j];
       }
       x[i] = t;
     }
     for (i=1;i<=ne;i++){/* 結果出力 */
       printf(" X_2_( %d ) = %15.6lE \n", i, x[i]);
     }
}
void dfg (double a[NE1][NE1], double b[NE1], double x[NE1], int ne)
{
     int i,j,k;
     double eps,t;
     eps = 1.0e-50;/* 特異性の判定値 */
     for (k=1;k<=ne;k++){/* 前進消去 */
       if ( fabs( a[k][k] )<=eps ){/* 除算を行う係数が小さく計算不可能なとき:行列は特異 */
         printf(" Matrix is singular at k= %d ",k );
         return;
       }
       for (i=k+1;i<=ne;i++){
         a[k][i] = a[k][i]/a[k][k];
       }
       b[k] = b[k]/a[k][k];
       for (i=k+1;i<=ne;i++){
         for (j=k+1;j<=ne;j++){
           a[i][j] = a[i][j] - a[i][k]*a[k][j];
         }
         b[i] = b[i] - a[i][k]*b[k];
       }
     }
     x[ne] = b[ne];/* 後退代入 */
     for (i=ne-1;i>=1;i--){
       t = b[i];
       for (j=i+1;j<=ne;j++){
         t = t - a[i][j] * x[j];
       }
       x[i] = t;
     }
     for (i=1;i<=ne;i++){/* 結果出力 */
       printf(" X_3_( %d ) = %15.6lE \n", i, x[i]);
     }
}
void dgjp (double a[NE1][NE1],double b[NE1],double x[NE1],int ne)
{
     int ipv[NE1];
     int i,j,k,l,ipvex;
     double eps,aipv;
     eps = 1.0e-75;/* 部分軸選択による特異性の判定値 */
     for (i=1;i<=ne;i++){/* 部分軸選択用配列の初期設定 */
       ipv[i] = i;
     }
     for (k=1;k<=ne;k++){
       l = k; aipv = fabs( a[ipv[l]][k] );
       for (i=k+1;i<=ne;i++){/* 部分軸選択 */
         if ( fabs( a[ipv[i]][k] ) > aipv ){
           l = i; aipv = fabs( a[ipv[l]][k] );
         }
       }
       if ( l != k ){/* 部分軸選択を行った方がよいと判断したら入れ替えを行う */
         ipvex = ipv[k]; ipv[k] = ipv[l] ; ipv[l] = ipvex;
       }
       if ( fabs( a[ipv[k]][k] )<=eps ){/*部分軸選択を行っても計算不可能なとき*/
         printf(" Matrix is singular at k= %d ",k );
         return;
       }
       for (i=k+1;i<=ne;i++){/* 消去過程 */
         a[ipv[k]][i] = a[ipv[k]][i]/a[ipv[k]][k];
       }
       b[ipv[k]] = b[ipv[k]]/a[ipv[k]][k];
       for (i=1;i<=ne;i++){
         if ( i != k ){
           for (j=k+1;j<=ne;j++){
             a[ipv[i]][j] = a[ipv[i]][j]-a[ipv[i]][k]*a[ipv[k]][j];
           }
           b[ipv[i]] = b[ipv[i]]-a[ipv[i]][k]*b[ipv[k]];
         }
       }
     }
     for (i=1;i<=ne;i++){/* 上の消去過程が終了するとIX=Bとなり解はすぐ求まる */
       x[i] = b[ipv[i]];
     }
     for (i=1;i<=ne;i++){/* 結果出力 */
       printf(" X_4_( %d ) = %15.6lE \n", i, x[i]);
     }
}
void dgbnd (double a1[NE1], double a2[NE1], double a3[NE1],
            double a4[NE1], double a5[NE1], double a6[NE1], double a7[NE1],
            double b[NE1], double x[NE1],
            int ny, int nxy, int ne)
{
     double a[NXY2][NE1];
     int ine,i,j,n,k;
     double aa,s;
     for (ine=1;ine<=ne;ine++){/* マトリックスAのゼロクリア */
       for (i=-nxy;i<=nxy;i++){
         a[i+nxy][ine] = 0.0;
       }
     }
     for (ine=1;ine<=ne;ine++){/* A1からA7までをaに格納する */
       a[ -ny+nxy][ine] = a1[ine];
       a[  -1+nxy][ine] = a2[ine];
       a[   0+nxy][ine] = a3[ine];
       a[   1+nxy][ine] = a4[ine];
       a[  ny+nxy][ine] = a5[ine];
       a[-nxy+nxy][ine] = a6[ine];
       a[ nxy+nxy][ine] = a7[ine];
     }
     for (i=1;i<=ne-1;i++){/* 前進消去 */
       if ( i <= ne-nxy ){
         for (j=1;j<=nxy;j++){
           aa = a[-j+nxy][i+j]/a[0+nxy][i];
           b[i+j] = b[i+j] - b[i]*aa;
           n=1;
           for (k=-j+1;k<=nxy-j;k++){
             a[k+nxy][i+j] = a[k+nxy][i+j]-a[n+nxy][i]*aa;
             n=n+1;
           }
         }
       }
       else{
         for (j=1;j<=ne-i;j++){
           aa = a[-j+nxy][i+j]/a[0+nxy][i];
           b[i+j] = b[i+j] - b[i]*aa;
           n=1;
           for (k=-j+1;k<=ne-i-j;k++){
             a[k+nxy][i+j] = a[k+nxy][i+j]-a[n+nxy][i]*aa;
             n=n+1;
           }
         }
       }
     }
     /* 係数行列の特異性を判定 */
     if ( fabs( a[0+nxy][ne] ) <= 1.0e-20 ){
        printf("  Matrix is singular : |A(0,NE)| < 1E-20 ");
     }
     x[ne] = b[ne] / a[0+nxy][ne];/* 後退代入 */
     for (i=ne-1;i>=1;i--){
       s = 0.0;
       if ( i >  ne-nxy ){
         for (n=1;n<=ne-i;n++){
           s = s + a[n+nxy][i] * x[i+n];
         }
         x[i] = ( b[i]-s ) / a[0+nxy][i];
       }
       else{
         for (n=1;n<=nxy;n++){
           s = s + a[n+nxy][i] * x[i+n];
         }
         x[i] = ( b[i]-s ) / a[0+nxy][i];
       }
     }
     for (ine=1;ine<=ne;ine++){/* 結果出力 */
         printf(" X_5_( %d ) = %15.6lE \n", ine, x[ine]);
     }
}
void hjacob (double a[NE1][NE1], double b[NE1], double xn[NE1], int ne)
{
     double x[NE1];
     int nitr,i,itr,j;
     double eitr,bnorm,rnorm,sum,zansa;
     nitr=200; eitr=1.0e-9;/*最大繰り返し回数(nitr)と収束判定値(eitr)の設定*/
     for (i=1;i<=ne;i++){/* 配列のゼロクリア */
       x[i] = 0.0;
     }
     bnorm = 0.0;/* Bの2乗ノルムの計算 */
     for (i=1;i<=ne;i++){
       bnorm = bnorm + b[i]*b[i];
     }
     itr = 0;/* ITR : 繰り返し回数のカウンタ */
     while (++itr>0){/* 反復回数あるいはzansaが条件を満たすまで */
       for (i=1;i<=ne;i++){/* 先に得られた値を旧値に設定し直す */
         x[i] = xn[i];
       }
       for (i=1;i<=ne;i++){/* 新値の計算 */
         sum = 0.0;
         for (j=1;j<=ne;j++){
          if (i != j){
            sum = sum + a[i][j]*x[j];
          }
          xn[i] = ( b[i]-sum )/a[i][i];
        }
       }
       rnorm = 0.0;/*収束判定のためのノルムの計算*/
       for (i=1;i<=ne;i++){/* RNORM : 残差ベクトル(R=AX-B)の2乗ノルム */
         sum = 0.0;
         for (j=1;j<=ne;j++){
           sum = sum + a[i][j]*xn[j];
         }
         rnorm = rnorm + ( b[i]-sum )*( b[i]-sum );
       }
       zansa = sqrt( rnorm/bnorm );/* 残差の計算 ZANSA = || R || / || B || */
       if (zansa <= eitr){/* 収束判定 */
         printf (" Converged : Total ITR =  %d \n",itr);
         goto HYOUJI;
       }
       if (itr > nitr){
         printf ("   Not converged ! \n");
         goto HYOUJI;
       }
     }
     HYOUJI:{/* 結果出力 */
       for (i=1;i<=ne;i++){
         printf(" X_6_( %d ) = %15.6lE \n", i, xn[i]);
       }
     }
}
void hsor (double a[NE1][NE1], double b[NE1], double x[NE1], int ne)
{
     double eitr,omg,bnorm,rnorm,xold,sum,xnew,zansa;
     int nitr,i,itr,j;
     char s[250];
     nitr=200; eitr=1.0e-9;/*最大繰り返し回数(nitr)と収束判定値(eitr)の設定*/
     /* 緩和係数の設定 ( OMG=1.0D0ならGauss-Seidel法 ) */
     printf(" Input omg ( if Gauss-Seidel, omg=1.0e0 ) ---> ");
     gets(s); omg=atof(s);
     bnorm = 0.0;/* Bの2乗ノルムの計算 */
     for (i=1;i<=ne;i++){
       bnorm = bnorm + b[i]*b[i];
     }
     itr = 0;/* ITR : 繰り返し回数のカウンタ */
     while (++itr>0){
       for (i=1;i<=ne;i++){
         xold = x[i]; sum = 0.0;
         for (j=1;j<=ne;j++){
           if ( i != j ){
             sum = sum + a[i][j]*x[j];
           }
         }
         xnew = ( b[i]-sum )/a[i][i];
         x[i] = xold + omg * ( xnew - xold );
       }
       rnorm = 0.0;/* 収束判定のためのノルムの計算 */
       for (i=1;i<=ne;i++){/* RNORM : 残差ベクトル(R=AX-B)の2乗ノルム */
         sum = 0.0;
         for (j=1;j<=ne;j++){
           sum = sum + a[i][j]*x[j];
         }
         rnorm = rnorm + ( b[i]-sum )*( b[i]-sum );
       }
       zansa = sqrt( rnorm / bnorm );/* 残差の計算 ZANSA = || R || / || B || */
       if (zansa <= eitr){/* 収束判定 */
         printf (" Converged : Total ITR =  %d \n",itr);
         goto HYOUJI;
       }
       if (itr > nitr){
         printf ("   Not converged ! \n");
         goto HYOUJI;
       }
     }
     HYOUJI:{/* 結果出力 */
       for (i=1;i<=ne;i++){
         printf(" X_7_( %d ) = %15.6lE \n", i, x[i]);
       }
     }
}
void sorbnd (double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1], 
             double a5[NE1], double a6[NE1], double a7[NE1], double b[NE1], double x[NE1],
             int ny, int nxy, int ne)
{
     int nitr,i,j;
     double eitr,omg,bnorm,rnorm,xold,sum,xnew,zansa;
     char s[250];
     nitr = 200; eitr = 1.0e-9;/*最大繰り返し回数(nitr)と収束判定値(eitr)の設定*/
     /* 緩和係数の設定 ( OMG=1.0D0ならGauss-Seidel法 ) */
     printf(" Input omg ( if Gauss-Seidel, omg=1.0e0 ) ---> ");
     gets(s); omg=atof(s);
     bnorm = 0.0;/* Bの2乗ノルムの計算 */
     for (i=1;i<=ne;i++){
       bnorm = bnorm + b[i]*b[i];
     }
     for (j=1;j<=nitr;j++){
       i=1;/* a3,a4,a5,a7 の範囲 */
         xold = x[i];
         sum = a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + omg * ( xnew - xold );
       for (i=2;i<=ny;i++) {/* a2,a3,a4,a5,a7の範囲 */
         xold = x[i];
         sum = a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + omg * ( xnew - xold );
       }
       for (i=ny+1;i<=nxy;i++) {/* a1-a5,a7 の範囲 */
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + omg * ( xnew - xold );
       }
       for (i=nxy+1;i<=ne-nxy;i++) {/* a6,a1-a5,a7 の範囲 */
         xold = x[i];
         sum = a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1]
              +a5[i]*x[i+ny]+a7[i]*x[i+nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + omg * ( xnew - xold );
       }
       for (i=ne-nxy+1;i<=ne-ny;i++) {/* a6, a1-a5 の範囲 */
         xold = x[i];
         sum = a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1]
              +a5[i]*x[i+ny];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + omg * ( xnew - xold );
       }
       for (i=ne-ny+1;i<=ne-1;i++) {/* a6, a1-a4 の範囲 */
         xold = x[i];
         sum = a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + omg * ( xnew - xold );
       }
       i=ne;/* a6, a1-a3 の範囲 */
         xold = x[i];
         sum=a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + omg * ( xnew - xold );
       rnorm = 0.0;/* RNORM : 残差ベクトル(R=AX-B)の2乗ノルム */
       for (i=1;i<=ne;i++){
         sum = 0.0;
         if (i==1){
           sum=a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];
         }
         if (i>=2) { if (i<=ny) {
           sum=a2[i]*x[i-1]+a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];}
         }
         if (i>=ny+1) { if (i<=nxy) {
           sum=a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]
              +a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];}
         }
         if (i>=nxy+1) { if (i<=ne-nxy) {
           sum=a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]
              +a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];}
         }
         if (i>=ne-nxy+1) { if (i<=ne-ny) {
           sum=a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]
              +a4[i]*x[i+1]+a5[i]*x[i+ny];}
         }
         if (i>=ne-ny+1) { if (i<=ne-1) {
           sum=a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]
              +a4[i]*x[i+1];}
         }
         if (i==ne) {
           sum=a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i];
         }
         rnorm= rnorm + (b[i]-sum)*(b[i]-sum);
       }
       zansa = sqrt(rnorm/bnorm);/* 残差の計算 ZANSA = || R || / || B || */
       if ( zansa <= eitr ){/* 収束判定 */
         printf (" Converged : Total ITR =  %d \n",j);
         goto label_700;
       }
     }
     printf (" Not converged ! ");/* NITRまで計算しても収束せず */
     label_700:{}/* 収束と判定されたときの分岐点 */
     for (i=1;i<=ne;i++){/* 結果出力 */
       printf(" X_8_( %d ) = %15.6lE \n", i, x[i]);
     }
}
void lblbnd (double at1[NE1] ,double at2[NE1], double at3[NE1],
             double at4[NE1], double at5[NE1], double at6[NE1], double at7[NE1],
             double bx[NE1], double xn[NE1],
             int nx, int ny, int nz, int nxy, int ne)
{
     double x1[NE1],xo[NE1],xold[NE1],a[NE1],b[NE1],c[NE1],d[NE1],u[NE1],y[NE1];
     double rnorm,bnorm,eitr,omg,sum,zansa;
     int k, inx, iy, ix, iz, iny, j, i, inz, inzz, nitr;
     char s[250];
     nitr = 200; eitr = 1.0e-9;/*最大繰り返し回数(nitr)と収束判定値(eitr)の設定*/
     /* 緩和係数の設定 ( OMG=1.0D0ならGauss-Seidel法 ) */
     printf(" Input omg ( if Gauss-Seidel, omg=1.0e0 ) ---> ");
     gets(s); omg=atof(s);
     bnorm = 0.0;/* Bの2乗ノルムの計算 */
     for (i=1;i<=ne;i++){
       bnorm = bnorm + bx[i]*bx[i];
     }
     for (k=1;k<=nitr;k++){
       inx = 1;/* x 軸方向への掃引 : トーマス法による*/
       for (iz=1;iz<=nz;iz++){
         for (iy=1;iy<=ny;iy++){
           for (ix=1;ix<=nx;ix++){
           inz = iy + (ix-1)*ny + (iz-1)*nxy;
           /* トーマス法のための係数A,B,C,Dの設定 : XNは最新のX */
           a[inx]=at1[inz]; b[inx]=at3[inz]; c[inx]=at5[inz]; d[inx]=bx[inz];
           if (inz-1 >= 1){
             d[inx]=d[inx]-at2[inz]*xn[inz-1];
           }
           if (inz+1 <= ne){
             d[inx]=d[inx]-at4[inz]*xn[inz+1];
           }
           if (inz-nxy >= 1){
             d[inx]=d[inx]-at6[inz]*xn[inz-nxy];
           }
           if (inz+nxy <= ne){
             d[inx]=d[inx]-at7[inz]*xn[inz+nxy];
           }
           xold[inz]=xn[inz];/* x方向への掃引の前の値をXOLDに保存 */
           xo[inz]=xn[inz];/* トーマス法で答えを求める前のXNをXOに保存 */
           inx=inx+1;
           }
         }
       }
       u[1]=c[1]/b[1];/* Ly=b を解く */
       for (j=2;j<=ne-1;j++){
         u[j]=c[j]/( b[j]-a[j]*u[j-1] );
       }
       y[1]=d[1]/b[1];
       for (j=2;j<=ne;j++){
         y[j]=( d[j]-a[j]*y[j-1] ) / ( b[j]-a[j]*u[j-1] );
       }
       x1[ne] = y[ne];/* Ux=y を解く */
       for (j=ne-1;j>=1;j--){
         x1[j] = y[j] - u[j]*x1[j+1];
       }
       inx = 1;
       for (iz=1;iz<=nz;iz++){
         for (iy=1;iy<=ny;iy++){
           for (ix=1;ix<=nx;ix++){
           inz=iy+(ix-1)*ny+(iz-1)*nxy;
           /* 得られたX1と反復前のXOにより最新のXNを緩和 */
           xn[inz]=(1.0-omg)*xo[inz]+omg*x1[inx];
           inx=inx+1;
           }
         }
       }
       iny = 1;/* y 軸方向への掃引 : トーマス法による */
       for (iz=1;iz<=nz;iz++){
         for (ix=1;ix<=nx;ix++){
           for (iy=1;iy<=ny;iy++){
           inz = iy + (ix-1)*nx + (iz-1)*nxy;
           /* トーマス法のための係数A,B,C,Dの設定 : XNは最新のX */
           a[iny]=at2[inz]; b[iny]=at3[inz]; c[iny]=at4[inz]; d[iny]=bx[inz];
           if (inz-ny >= 1){
             d[iny]=d[iny]-at1[inz]*xn[inz-ny];
           }
           if (inz+ny <= ne){
             d[iny]=d[iny]-at5[inz]*xn[inz+ny];
           }
           if (inz-nxy >= 1){
             d[iny]=d[iny]-at6[inz]*xn[inz-nxy];
           }
           if (inz+nxy <= ne){
             d[iny]=d[iny]-at7[inz]*xn[inz+nxy];
           }
           xo[inz]=xn[inz];/* トーマス法で答えを求める前のXNをXOに保存 */
           iny=iny+1;
           }
         }
       }
       u[1] = c[1] / b[1];/* Ly=b を解く */
       for (j=2;j<=ne-1;j++){
         u[j] = c[j]/( b[j]-a[j]*u[j-1] );
       }
       y[1] = d[1] / b[1];
       for (j=2;j<=ne;j++){
         y[j] = ( d[j]-a[j]*y[j-1] ) / ( b[j]-a[j]*u[j-1] );
       }
       x1[ne] = y[ne];/* Ux=y を解く */
       for (j=ne-1;j>=1;j--){
         x1[j] = y[j] - u[j]*x1[j+1];
       }
       iny = 1;
       for (iz=1;iz<=nz;iz++){
         for (ix=1;ix<=nx;ix++){
           for (iy=1;iy<=ny;iy++){
           inz = iy + (ix-1)*ny + (iz-1)*nxy;
           /* 得られたX1と反復前のXOにより最新のXNを緩和 */
           xn[inz] = (1.0-omg)*xo[inz]+omg*x1[iny];
           iny = iny + 1;
           }
         }
       }
       inzz = 1;/* z 軸方向への掃引 : トーマス法による */
       for (iy=1;iy<=ny;iy++){
         for (ix=1;ix<=nx;ix++){
           for (iz=1;iz<=nz;iz++){
           inz = iy + (ix-1)*nx + (iz-1)*nxy;
           /* トーマス法のための係数A,B,C,Dの設定 : XNは最新のX */
           a[inzz]=at6[inz]; b[inzz]=at3[inz]; c[inzz]=at7[inz]; d[inzz]=bx[inz];
           if (inz-ny >= 1){
             d[inzz]=d[inzz]-at1[inz]*xn[inz-ny];
           }
           if (inz+ny <= ne){
             d[inzz]=d[inzz]-at5[inz]*xn[inz+ny];
           }
           if (inz-1 >= 1){
             d[inzz]=d[inzz]-at2[inz]*xn[inz-1];
           }
           if (inz+1 <= ne){
             d[inzz]=d[inzz]-at4[inz]*xn[inz+1];
           }
           xo[inz]=xn[inz];/* トーマス法で答えを求める前のXNをXOに保存 */
           inzz=inzz+1;
           }
         }
       }
       u[1] = c[1] / b[1];/* Ly=b を解く */
       for (j=2;j<=ne-1;j++){
         u[j] = c[j]/( b[j]-a[j]*u[j-1] );
       }
       y[1] = d[1] / b[1];
       for (j=2;j<=ne;j++){
         y[j] = ( d[j]-a[j]*y[j-1] ) / ( b[j]-a[j]*u[j-1] );
       }
       x1[ne] = y[ne];/* Ux=y を解く */
       for (j=ne-1;j>=1;j--){
         x1[j] = y[j] - u[j]*x1[j+1];
       }
       inzz = 1;
       for (iy=1;iy<=ny;iy++){
         for (ix=1;ix<=nx;ix++){
           for (iz=1;iz<=nz;iz++){
           inz = iy + (ix-1)*ny + (iz-1)*nxy;
           /* 得られたX1と反復前のXOにより最新のXNを緩和 */
           xn[inz] = (1.0-omg)*xo[inz]+omg*x1[inzz];
           inzz = inzz + 1;
           }
         }
       }
       rnorm = 0.0;/* RNORM : 残差ベクトル(R=AX-B)の2乗ノルム */
       for (i=1;i<=ne;i++){
         sum = 0.0;
         if (i==1){
           sum=at3[i]*xn[i]+at4[i]*xn[i+1]+at5[i]*xn[i+ny]+at7[i]*xn[i+nxy];
         }
         if (i>=2) { if (i<=ny) {
           sum=at2[i]*xn[i-1]+at3[i]*xn[i]+at4[i]*xn[i+1]
              +at5[i]*xn[i+ny]+at7[i]*xn[i+nxy];}
         }
         if (i>=ny+1) { if (i<=nxy) {
           sum=at1[i]*xn[i-ny]+at2[i]*xn[i-1]+at3[i]*xn[i]
              +at4[i]*xn[i+1]+at5[i]*xn[i+ny]+at7[i]*xn[i+nxy];}
         }
         if (i>=nxy+1) { if (i<=ne-nxy) {
           sum=at6[i]*xn[i-nxy]+at1[i]*xn[i-ny]+at2[i]*xn[i-1]+at3[i]*xn[i]
              +at4[i]*xn[i+1]+at5[i]*xn[i+ny]+at7[i]*xn[i+nxy];}
         }
         if (i>=ne-nxy+1) { if (i<=ne-ny) {
           sum=at6[i]*xn[i-nxy]+at1[i]*xn[i-ny]+at2[i]*xn[i-1]+at3[i]*xn[i]
              +at4[i]*xn[i+1]+at5[i]*xn[i+ny];}
         }
         if (i>=ne-ny+1) { if (i<=ne-1) {
           sum=at6[i]*xn[i-nxy]+at1[i]*xn[i-ny]+at2[i]*xn[i-1]+at3[i]*xn[i]
              +at4[i]*xn[i+1];}
         }
         if (i==ne) {
           sum=at6[i]*xn[i-nxy]+at1[i]*xn[i-ny]+at2[i]*xn[i-1]+at3[i]*xn[i];
         }
         rnorm= rnorm + (bx[i]-sum)*(bx[i]-sum);
       }
       zansa = sqrt(rnorm/bnorm);/* 残差の計算 ZANSA = || R || / || B || */
       if ( zansa <= eitr ){/* 収束判定 */
         printf (" Converged : Total ITR =  %d \n",k);
         goto label_700;
       }
     }
     printf (" Not converged ! ");/* NITRまで計算しても収束せず */
     label_700:{}/* 収束と判定されたときの分岐点 */
     for (i=1;i<=ne;i++){/* 結果出力 */
       printf(" X_9_( %d ) = %15.6lE \n", i, xn[i]);
     }
}
void kcr (double a[NE1][NE1], double b[NE1], double x[NE1], int ne)
{
     double r[NE1], p[NE1], ap[NE1], ar[NE1];/* 作業用配列 */
     double eitr,bnorm,rap,apap,alp,rnorm,sum,zansa,arap,beta;
     int nitr,j,i,itr;
     nitr = 200;/* 最大繰り返し回数の設定 */
/* 収束判定条件(変数 zansa がこの値以下になれば nitr 以下でも収束と判定する) */
     eitr = 1.0e-9;
     for (j=1;j<=ne;j++){/* 配列のゼロクリア */
       r[j] = 0.0; p[j] = 0.0; ap[j] = 0.0; ar[j] = 0.0;
     }
     profmv (a,x,r,ne);/* 配列 r に AX を代入 */
     bnorm = 0.0;/* bのノルムの計算 */
     for (i=1;i<=ne;i++){
       bnorm = bnorm + b[i]*b[i];
     }
     for (i=1;i<=ne;i++){/*r_{0}とp_{0}(初期値)の設定*/
       r[i] = b[i] - r[i]; p[i] = r[i];/* r_{0}=b-AX, p_{0}=r_{0} */
     }
     profmv (a,p,ap,ne);/* apに A p_{0} を代入(以降のAPは式(A)で求める)*/
     for (itr=1; itr<=nitr; itr++){/* 繰り返し計算 */
       rap = provv (r,ap,ne);/* ( r_{k}, A p_{k} )の計算 => RAP */
       apap = provv (ap,ap,ne);/* ( A p_{k}, A p_{k} )の計算 => APAP */
       alp = rap / apap; rnorm = 0.0;/*α_{k}の計算*/
       for (i=1;i<=ne;i++){
         x[i]=x[i]+alp*p[i]; r[i]=r[i]-alp*ap[i];
       }/* x_{k+1}=x_{k}+α_{k}p_{k},r_{k+1}=r_{k}-α_{k}Ap_{k} */
       rnorm = 0.0;/* RNORM : 残差ベクトル(R=AX-B)の2乗ノルム */
       for (i=1;i<=ne;i++){
         sum = 0.0;
         for (j=1;j<=ne;j++){
           sum = sum + a[i][j]*x[j];
         }
         rnorm = rnorm + ( b[i]-sum )*( b[i]-sum );
       }
       zansa = sqrt( rnorm / bnorm );/* 残差の計算 ZANSA = || R || / || B || */
       if (zansa <= eitr){/* 収束判定 */
         printf (" Converged : Total ITR =  %d \n",itr);
         goto label_700;
       }
       profmv (a,r,ar,ne);
       arap = provv( ar,ap,ne );/* ( A r_{k+1}, A p_{k} )の計算 => ARAP */
       beta = - arap / apap;/* β_{k}= - ( A r_{k+1}, A p_{k} ) / ( A p_{k}, A p_{k} ) */
       for (i=1;i<=ne;i++){
         p[i]=r[i]+beta*p[i]; ap[i]=ar[i]+beta*ap[i];
       }/*p_{k+1}=r_{k+1}+β_{k}p_{k},Ap_{k+1}=Ar_{k+1}+β_{k}A p_{k}<----- 式(A) */
     }
     printf ("   Not converged ! \n");/* nitr まで計算しても収束せず */
     label_700:{}/* 収束とみなされたときの分岐 */
     for (i=1;i<=ne;i++){/* 結果出力 */
       printf(" X_10_( %d ) = %15.6lE \n", i, x[i]);
     }
}
void profmv (double a[NE1][NE1], double b[NE1], double c[NE1] ,int ne)
{
     int i,j;
     double s;
     for (i=1;i<=ne;i++){
       s = 0.0;
       for (j=1;j<=ne;j++){
         s = s + a[i][j]*b[j];
       }
       c[i] = s;
     }
}
double provv (double a[NE1], double b[NE1], int ne)
{
      int i;
      double c;

      c = 0.0;
      for (i=1;i<=ne;i++){
        c = c + a[i]*b[i];
      }
      return(c);
}
void kcrbnd (double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1],
             double a5[NE1], double a6[NE1], double a7[NE1],
             double b[NE1], double x[NE1],
             int ny, int nxy, int ne)
{
     double r[NE1], p[NE1], ap[NE1], ar[NE1];/* 作業用配列 */
     double eitr,bnorm,rap,apap,alp,rnorm,sum,zansa,arap,beta;
     int nitr,j,i,itr;
     nitr = 200;/* 最大繰り返し回数の設定 */
/* 収束判定条件(変数 zansa がこの値以下になれば nitr 以下でも収束と判定する) */
     eitr = 1.0e-9;
     for (j=1;j<=ne;j++){/* 配列のゼロクリア */
       r[j] = 0.0; p[j] = 0.0; ap[j] = 0.0; ar[j] = 0.0;
     }
     probmv (a1,a2,a3,a4,a5,a6,a7,x,r,ny,nxy,ne);/* 配列 r に AX を代入 */
     bnorm = 0.0;/* bのノルムの計算 */
     for (i=1;i<=ne;i++){
       bnorm = bnorm + b[i]*b[i];
     }
     for (i=1;i<=ne;i++){/*r_{0}とp_{0}(初期値)の設定*/
       r[i] = b[i] - r[i]; p[i] = r[i];/* r_{0}=b-AX, p_{0}=r_{0} */
     }
     probmv (a1,a2,a3,a4,a5,a6,a7,p,ap,ny,nxy,ne);/* apに A p_{0} を代入(以降のAPは式(A)で求める)*/
     for (itr=1; itr<=nitr; itr++){/* 繰り返し計算 */
       rap = provv (r,ap,ne);/* ( r_{k}, A p_{k} )の計算 => RAP */
       apap = provv (ap,ap,ne);/* ( A p_{k}, A p_{k} )の計算 => APAP */
       alp = rap / apap; rnorm = 0.0;/*α_{k}の計算*/
       for (i=1;i<=ne;i++){
         x[i]=x[i]+alp*p[i]; r[i]=r[i]-alp*ap[i];
       }/* x_{k+1}=x_{k}+α_{k}p_{k},r_{k+1}=r_{k}-α_{k}Ap_{k} */
       rnorm = 0.0;/* RNORM : 残差ベクトル(R=AX-B)の2乗ノルム */
       for (i=1;i<=ne;i++){
         sum = 0.0;
         if (i==1){
           sum=a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];
         }
         if (i>=2) { if (i<=ny) {
           sum=a2[i]*x[i-1]+a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];}
         }
         if (i>=ny+1) { if (i<=nxy) {
           sum=a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]
              +a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];}
         }
         if (i>=nxy+1) { if (i<=ne-nxy) {
           sum=a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]
              +a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];}
         }
         if (i>=ne-nxy+1) { if (i<=ne-ny) {
           sum=a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]
              +a4[i]*x[i+1]+a5[i]*x[i+ny];}
         }
         if (i>=ne-ny+1) { if (i<=ne-1) {
           sum=a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]
              +a4[i]*x[i+1];}
         }
         if (i==ne) {
           sum=a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i];
         }
         rnorm= rnorm + (b[i]-sum)*(b[i]-sum);
       }
       zansa = sqrt( rnorm / bnorm );/* 残差の計算 ZANSA = || R || / || B || */
       if (zansa <= eitr){/* 収束判定 */
         printf (" Converged : Total ITR =  %d \n",itr);
         goto label_700;
       }
       probmv (a1,a2,a3,a4,a5,a6,a7,r,ar,ny,nxy,ne);
       arap = provv( ar,ap,ne );/* ( A r_{k+1}, A p_{k} )の計算 => ARAP */
       beta = - arap / apap;/* β_{k}= - ( A r_{k+1}, A p_{k} ) / ( A p_{k}, A p_{k} ) */
       for (i=1;i<=ne;i++){
         p[i]=r[i]+beta*p[i]; ap[i]=ar[i]+beta*ap[i];
       }/*p_{k+1}=r_{k+1}+β_{k}p_{k},Ap_{k+1}=Ar_{k+1}+β_{k}A p_{k}<----- 式(A) */
     }
     printf ("   Not converged ! \n");/* nitr まで計算しても収束せず */
     label_700:{}/* 収束とみなされたときの分岐 */
     for (i=1;i<=ne;i++){/* 結果出力 */
       printf(" X_11_( %d ) = %15.6lE \n", i, x[i]);
     }
}
void probmv (double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1],
             double a5[NE1], double a6[NE1], double a7[NE1],
             double b[NE1], double c[NE1],
             int ny, int nxy, int ne)
{
     int i;

     i=1; c[i]=a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny]+a7[i]*b[i+nxy];

     for (i=2;i<=ny;i++){
       c[i]=a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny]+a7[i]*b[i+nxy];
     }

     for (i=ny+1;i<=nxy;i++){
       c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny]+a7[i]*b[i+nxy];
     }

     for (i=nxy+1;i<=ne-nxy;i++){
       c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny]
            +a6[i]*b[i-nxy]+a7[i]*b[i+nxy];
     }

     for (i=ne-nxy+1;i<=ne-ny;i++){
       c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny]
           +a6[i]*b[i-nxy];
     }
     for (i=ne-ny+1;i<=ne-1;i++){
       c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1]+a6[i]*b[i-nxy];
     }

     i=ne; c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i]+a6[i]*b[i-nxy];

}
void kbicg (double a[NE1][NE1], double b[NE1], double x[NE1], int ne)
{
     double r[NE1],ap[NE1],at[NE1],p[NE1],s[NE1],t[NE1];/* 作業用配列 */
     double eitr,bnorm,sr1,sap,alpha,att,atat,xi,rnorm,sum,zansa,sr2,beta;
     int nitr,i,j,m;
     nitr = 200;/* 最大繰り返し回数の設定 */
     /* 収束判定条件(変数 ZANSA がこの値以下になれば NITR 以下でも収束と判定する) */
     eitr = 1.0e-9;
     bnorm = 0.0e0;/* Bのノルムの計算 */
     for (i=1;i<=ne;i++){
       bnorm = bnorm + b[i]*b[i];
     }
     for (j=1;j<=ne;j++){/* 配列のゼロクリア */
       r[j]=0.0; ap[j]=0.0; at[j]=0.0; p[j]=0.0; s[j]=0.0; t[j]=0.0;
     }
     profmv (a,x,r,ne);/* R に AX を代入 */
     for (i=1;i<=ne;i++){/* r_{0}とp_{0}(初期値)，そして s=r_{0} の設定 */
       r[i]=b[i]-r[i]; p[i]=r[i]; s[i]=r[i];
     }
     for (j=1;j<=nitr;j++){/* 繰り返し計算 */
       sr1 = provv (s,r,NE);/* ( s, r_{k} ) の計算 => SR1 */
       profmv (a,p,ap,ne);/* A p_{k} の計算 => AP(NE) */
       sap = provv (s,ap,ne);/* ( s, A p_{k} ) の計算 => SAP */
       alpha = sr1/sap;/* α_{k} = ( s, r_{k} ) / ( s, A p_{k} ) */
       for (i=1;i<=ne;i++){
         t[i] = r[i] - alpha*ap[i];/* t_{k} = r_{k} - α_{k} A p_{k} */
       }
       profmv (a,t,at,ne);/* A t_{k} の計算 => AT(NE) */
       att = provv (at,t,ne);/* ( A t_{k}, t_{k} ) の計算 => ATT */
       atat = provv (at,at,ne);/* ( A t_{k}, A t_{k} ) の計算 => ATAT */
       xi = att/atat;/* ξ_{k} = ( A t_{k}, t_{k} ) / ( A t_{k}, A t_{k} ) */
       for (i=1;i<=ne;i++){
         x[i]=x[i]+alpha*p[i]+xi*t[i];/* x_{k+1} = x_{k} + α_{k} p_{k} + ξ_{k} t_{k} */
         r[i]=t[i]-xi*at[i];/* r_{k+1} = t_{k} - ξ_{k} A t_{k} */
       }
       rnorm=0.0;/* RNORM : 残差ベクトル(R=AX-B)の2乗ノルム */
       for (i=1;i<=ne;i++){
         sum = 0.0;
         for (m=1;m<=ne;m++){
           sum = sum + a[i][m]*x[m];
         }
         rnorm = rnorm + (b[i]-sum)*(b[i]-sum);
       }
       zansa = sqrt(rnorm/bnorm);/* 残差の計算 ZANSA = || R || / || B || */
       if (zansa <= eitr){/*       ZANSA が EITR 以下なら収束とみなして 900 へ */
         printf (" Converged : Total ITR =  %d \n",j);
         goto label_900;
        }
        /* 収束せずの場合 : β_{k}と p_{k+1} を求めて繰り返し計算 */
        sr2 = provv (s,r,ne);/* ( s, r_{k+1} ) の計算 => SR2 */
        /* β_{k} = ( α_{k} / ξ_{k} ) * ( s, r_{k+1} )/( s, r_{k} ) */
        beta = (alpha / xi) * (sr2 / sr1);
       for (i=1;i<=ne;i++){
         /* p_{k+1} = r_{k+1} + β_{k}( p_{k}-ξ_{k} A p_{k} ) */
         p[i] = r[i] + beta * ( p[i] - xi*ap[i] );
       }
     }
     printf (" Not converged ! \n");/* NITR まで計算しても収束せず */
     label_900:{}/* 収束と判定されたときの分岐点 */
     for (i=1;i<=ne;i++){/* 結果出力 */
       printf(" X_12_( %d ) = %15.6lE \n", i, x[i]);
     }
}
void kbibnd (double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1], double a5[NE1],
             double a6[NE1], double a7[NE1], double b[NE1], double x[NE1],
             int ny, int nxy, int ne)
{
     double r[NE1],ap[NE1],at[NE1],p[NE1],s[NE1],t[NE1];/* 作業用配列 */
     double eitr,bnorm,sr1,sap,alpha,att,atat,xi,rnorm,sum,zansa,sr2,beta;
     int nitr,i,j;
     nitr = 200;/* 最大繰り返し回数の設定 */
     /* 収束判定条件(変数 ZANSA がこの値以下になれば NITR 以下でも収束と判定する) */
     eitr = 1.0e-9;
     bnorm = 0.0e0;/* Bのノルムの計算 */
     for (i=1;i<=ne;i++){
       bnorm = bnorm + b[i]*b[i];
     }
     for (j=1;j<=ne;j++){/* 配列のゼロクリア */
       r[j]=0.0; ap[j]=0.0; at[j]=0.0; p[j]=0.0; s[j]=0.0; t[j]=0.0;
     }
     probmv (a1,a2,a3,a4,a5,a6,a7,x,r,ny,nxy,ne);/* R に AX を代入 */
     for (i=1;i<=ne;i++){/* r_{0}とp_{0}(初期値)，そして s=r_{0} の設定 */
       r[i]=b[i]-r[i]; p[i]=r[i]; s[i]=r[i];
     }
     for (j=1;j<=nitr;j++){/* 繰り返し計算 */
       sr1 = provv (s,r,NE);/* ( s, r_{k} ) の計算 => SR1 */
       probmv (a1,a2,a3,a4,a5,a6,a7,p,ap,ny,nxy,ne);/* A p_{k} の計算 => AP(NE) */
       sap = provv (s,ap,ne);/* ( s, A p_{k} ) の計算 => SAP */
       alpha = sr1/sap;/* α_{k} = ( s, r_{k} ) / ( s, A p_{k} ) */
       for (i=1;i<=ne;i++){
         t[i] = r[i] - alpha*ap[i];/* t_{k} = r_{k} - α_{k} A p_{k} */
       }
       probmv (a1,a2,a3,a4,a5,a6,a7,t,at,ny,nxy,ne);/* A t_{k} の計算 => AT(NE) */
       att = provv (at,t,ne);/* ( A t_{k}, t_{k} ) の計算 => ATT */
       atat = provv (at,at,ne);/* ( A t_{k}, A t_{k} ) の計算 => ATAT */
       xi = att/atat;/* ξ_{k} = ( A t_{k}, t_{k} ) / ( A t_{k}, A t_{k} ) */
       for (i=1;i<=ne;i++){
         x[i]=x[i]+alpha*p[i]+xi*t[i];/* x_{k+1} = x_{k} + α_{k} p_{k} + ξ_{k} t_{k} */
         r[i]=t[i]-xi*at[i];/* r_{k+1} = t_{k} - ξ_{k} A t_{k} */
       }
       rnorm=0.0;/* RNORM : 残差ベクトル(R=AX-B)の2乗ノルム */
       for (i=1;i<=ne;i++){
         sum = 0.0;
         if (i==1){
           sum=a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];
         }
         if (i>=2) { if (i<=ny) {
           sum=a2[i]*x[i-1]+a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];}
         }
         if (i>=ny+1) { if (i<=nxy) {
           sum=a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]
              +a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];}
         }
         if (i>=nxy+1) { if (i<=ne-nxy) {
           sum=a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]
              +a4[i]*x[i+1]+a5[i]*x[i+ny]+a7[i]*x[i+nxy];}
         }
         if (i>=ne-nxy+1) { if (i<=ne-ny) {
           sum=a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]
              +a4[i]*x[i+1]+a5[i]*x[i+ny];}
         }
         if (i>=ne-ny+1) { if (i<=ne-1) {
           sum=a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]
              +a4[i]*x[i+1];}
         }
         if (i==ne) {
           sum=a6[i]*x[i-nxy]+a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i];
         }
         rnorm= rnorm + (b[i]-sum)*(b[i]-sum);
       }
       zansa = sqrt(rnorm/bnorm);/* 残差の計算 ZANSA = || R || / || B || */
       if (zansa <= eitr){/*       ZANSA が EITR 以下なら収束とみなして 900 へ */
         printf (" Converged : Total ITR =  %d \n",j);
         goto label_900;
        }
        /* 収束せずの場合 : β_{k}と p_{k+1} を求めて繰り返し計算 */
        sr2 = provv (s,r,ne);/* ( s, r_{k+1} ) の計算 => SR2 */
        /* β_{k} = ( α_{k} / ξ_{k} ) * ( s, r_{k+1} )/( s, r_{k} ) */
        beta = (alpha / xi) * (sr2 / sr1);
       for (i=1;i<=ne;i++){
         /* p_{k+1} = r_{k+1} + β_{k}( p_{k}-ξ_{k} A p_{k} ) */
         p[i] = r[i] + beta * ( p[i] - xi*ap[i] );
       }
     }
     printf (" Not converged ! \n");/* NITR まで計算しても収束せず */
     label_900:{}/* 収束と判定されたときの分岐点 */
     for (i=1;i<=ne;i++){/* 結果出力 */
       printf(" X_13_( %d ) = %15.6lE \n", i, x[i]);
     }
}
