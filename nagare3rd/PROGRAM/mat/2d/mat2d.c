/********************************************************************
* ファイル名  ：mat2d.c                                             *
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
#define NX 5 /* NX:x方向格子数 */
#define NY 5 /* NY:y方向格子数 */
#define NE 25 /* NE:全格子数=NX*NY */
#define NX1 6 /* NX1=NX+1 */
#define NY1 6 /* NY1=NY+1 */
#define NE1 26 /* NE1=NE+1 */
#define NY2 11 /* NY2=NY*2+1 */
void dlu (double [NE1][NE1], double [NE1], double [NE1], int ne);
void dfgp (double [NE1][NE1], double [NE1], double [NE1], int ne);
void dfg (double [NE1][NE1], double [NE1], double [NE1], int ne);
void dgjp (double [NE1][NE1], double [NE1], double [NE1], int ne);
void dgbnd (double [NE1], double [NE1], double [NE1],
            double [NE1], double [NE1], double [NE1],
            double [NE1], int ny, int ne);
void hjacob (double [NE1][NE1], double [NE1], double [NE1], int ne);
void hgs (double [NE1][NE1], double [NE1], double [NX1][NY1], int ne);
void hsor (double [NE1][NE1], double [NE1], double [NE1], int ne);
void sorbnd (double [NE1], double [NE1], double [NE1], double [NE1], 
             double [NE1], double [NE1], double [NE1],
             int ny, int ne);
void lblbnd (double [NE1] ,double [NE1], double [NE1],
             double [NE1], double [NE1], double [NE1], double [NE1],
             int nx, int ny, int ne);
void kcr (double [NE1][NE1], double [NE1], double [NE], int ne);
void profmv (double [NE1][NE1], double [NE1], double [NE1] ,int ne);
double provv (double[NE1], double[NE1], int ne);
void kcrbnd (double [NE1], double [NE1], double [NE1],
             double [NE1], double [NE1], double [NE1],
             double [NE1], int ny, int ne);
void probmv (double [NE1], double [NE1], double [NE1], double [NE1],
             double [NE1], double [NE1], double [NE1],
             int ny, int ne);
void kbicg (double [NE1][NE1], double [NE1], double [NE1], int ne);
void kbibnd  (double [NE1], double [NE1], double [NE1],
              double [NE1], double [NE1], double [NE1],
              double [NE1], int ny, int ne);
void main()
{
     double a[NE1][NE1], b[NE1], x1[NE1], x2[NX1][NY1];
/* バンドマトリックス用配列の定義 */
     double at1[NE1],at2[NE1],at3[NE1],at4[NE1],at5[NE1];
     double delta=1.0;/* 上記離散化方程式中のΔ=1を設定 */
     int i,ii,ix,iy,nx,ny,ne;
     int ij;
     double b1,b2,b3,b4,b5;
     char s[250];
     nx=NX; ny=NY; ne=NE;
     next_cal:{};/* 解法の変更毎の戻り点 */
/* A,Bを設定するためにあらかじめx2の値を格子番号と同じになるようにする */
     i = 1;
     for (ix=1;ix<=nx;ix++){
       for (iy=1;iy<=ny;iy++){
         x2[ix][iy]=i; i=i+1;
       }
     }
     for (i=1;i<=ne;i++){/* x2以外の配列のゼロクリア */
       at1[i]=0.0; at2[i]=0.0; at3[i]=0.0; at4[i]=0.0; at5[i]=0.0;
       b[i]=0.0; x1[i]=0;
       for (ii=1;ii<=ne;ii++){
         a[i][ii] = 0.0;
       }
     }
/* X_{i,j}=k=(i-1)*NY+jとなるようにA,Bを求める */
     i=1;
     for (ix=1;ix<=nx;ix++){
       for (iy=1;iy<=ny;iy++){
         if (ix==1){/*f(i-1,j)が計算領域外（左側）の境界条件の処理*/
           at1[i]=0.0;/*f(0,j)=0として処理する*/
           b1=0.0;/*AT1*f(i-1,j)をb1とする*/
         }
         else {/*f(i-1,j)が計算領域内の場合*/
           at1[i] = 1.0/(delta*delta) - 1.0/(2.0*delta);
           a[i][i-ny] = at1[i];
           b1 = x2[ix-1][iy] * ( 1.0/(delta*delta) - 1.0/(2.0*delta) );
         }
         at3[i] = -2.0/(delta*delta)-2.0/(delta*delta);/*f(i,j)は常に計算領域内*/
         a[i][i] = at3[i];/*AT3*f(i,j)をb3とする*/
         b3 = x2[ix][iy] * ( -2.0/(delta*delta) - 2.0/(delta*delta) );
         if (ix==nx){/*f(i+1,j)が計算領域外（右側）の境界条件の処理*/
           at5[i] = 0.0;/*f(i+1,j)=0として処理*/
           b5 = 0.0;/*AT5*f(i+1,j)=b5とする*/
         }
         else {/*f(i+1,j)が計算領域内の場合*/
           at5[i] = 1.0/(delta*delta) + 1.0/(2.0*delta);
           a[i][i+ny] = at5[i];
           b5 = x2[ix+1][iy]*( 1.0/(delta*delta) + 1.0/(2.0*delta) );
         }
         if (iy==1){/*f(i,j-1)が計算領域外（下側）の境界条件の処理*/
           at2[i] = 0.0;/*f(i,j-1)=0として処理*/
           b2 = 0.0;/*AT2*f(i,j-1)=b2とする*/
         }
         else {/*f(i,j-1)が計算領域内の場合*/
           at2[i] = 1.0/(delta*delta) - 1.0/(2.0*delta);
           a[i][i-1] = at2[i];
           b2 = x2[ix][iy-1]*( 1.0/(delta*delta) - 1.0/(2.0*delta) );
         }
         if (iy==ny){/*f(i,j+1)が計算領域外（上側）の境界条件の処理*/
           at4[i] = 0.0; b4 = 0.0;/*f(i,j+1)=0として処理*/
         }/*AT4*f(i,j+1)=b4とする*/
         else {/*f(i,j+1)が計算領域内の場合*/
           at4[i] = 1.0/(delta*delta) + 1.0/(2.0*delta);
           a[i][i+1] = at4[i];
           b4 = x2[ix][iy+1]*( 1.0/(delta*delta) + 1.0/(2.0*delta) );
         }
         b[i] = b1+b2+b3+b4+b5; i=i+1;/*bの計算*/
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
       case  5 : dgbnd (at1,at2,at3,at4,at5,b,x1,ny,ne); break;
       case  6 : hjacob (a,b,x1,ne); break;
       case  7 : hsor (a,b,x1,ne); break;
       case  8 : sorbnd (at1,at2,at3,at4,at5,b,x1,ny,ne); break;
       case  9 : lblbnd (at1,at2,at3,at4,at5,b,x1,nx,ny,ne); break;
       case 10 : kcr (a,b,x1,ne); break;
       case 11 : kcrbnd (at1,at2,at3,at4,at5,b,x1,ny,ne); break;
       case 12 : kbicg (a,b,x1,ne); break;
       case 13 : kbibnd (at1,at2,at3,at4,at5,b,x1,ny,ne); break;
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
void dgbnd (double at1[NE1], double at2[NE1], double at3[NE1],
            double at4[NE1], double at5[NE1], double b[NE1],
            double x[NE1], int ny, int ne)
{
     double at[NY2][NE1];
     int ine,i,j,n,k;
     double aa,s;
     for (ine=1;ine<=ne;ine++){/* マトリックスAのゼロクリア */
       for (i=-ny;i<=ny;i++){
         at[i+ny][ine] = 0.0;
       }
     }
     for (ine=1;ine<=ne;ine++){/* AT1からAT5を at に格納 */
       at[-ny+ny][ine] = at1[ine];
       at[ -1+ny][ine] = at2[ine];
       at[  0+ny][ine] = at3[ine];
       at[  1+ny][ine] = at4[ine];
       at[ ny+ny][ine] = at5[ine];
     }
     for (i=1;i<=ne-1;i++){/* 前進消去 */
       if ( i <= ne-ny ){
         for (j=1;j<=ny;j++){
           aa = at[-j+ny][i+j]/at[0+ny][i];
           b[i+j] = b[i+j] - b[i]*aa;
           n=1;
           for (k=-j+1;k<=ny-j;k++){
             at[k+ny][i+j] = at[k+ny][i+j]-at[n+ny][i]*aa;
             n=n+1;
           }
         }
       }
       else{
         for (j=1;j<=ne-i;j++){
           aa = at[-j+ny][i+j]/at[0+ny][i];
           b[i+j] = b[i+j] - b[i]*aa;
           n=1;
           for (k=-j+1;k<=ne-i-j;k++){
             at[k+ny][i+j] = at[k+ny][i+j]-at[n+ny][i]*aa;
             n=n+1;
           }
         }
       }
     }
     /* 係数行列の特異性を判定 */
     if ( fabs( at[0+ny][ne] ) <= 1.0e-20 ){
        printf("  Matrix is singular : |A(0,NE)| < 1E-20 ");
     }
     x[ne] = b[ne] / at[0+ny][ne];/* 後退代入 */
     for (i=ne-1;i>=1;i--){
       s = 0.0;
       if ( i >  ne-ny ){
         for (n=1;n<=ne-i;n++){
           s = s + at[n+ny][i] * x[i+n];
         }
         x[i] = ( b[i]-s ) / at[0+ny][i];
       }
       else{
         for (n=1;n<=ny;n++){
           s = s + at[n+ny][i] * x[i+n];
         }
         x[i] = ( b[i]-s ) / at[0+ny][i];
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
             double a5[NE1], double b[NE1], double x[NE1],
             int ny, int ne)
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
       i=1;/* a3,a4,a5 の範囲 */
         xold = x[i];
         sum = a4[i]*x[i+1]+a5[i]*x[i+ny];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + omg * ( xnew - xold );
       for (i=2;i<=ny;i++){/* a2,a3,a4,a5の範囲 */
         xold = x[i];
         sum = a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + omg * ( xnew - xold );
       }
       for (i=ny+1;i<=ne-ny;i++){/* a1-a5 の範囲 */
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + omg * ( xnew - xold );
       }
       for (i=ne-ny+1;i<=ne-1;i++){/* a1-a4 の範囲 */
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + omg * ( xnew - xold );
       }
       i=ne;/* a1-a3 の範囲 */
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + omg * ( xnew - xold );
       rnorm=0.0;/* RNORM : 残差ベクトル(R=AX-B)の2乗ノルム */
       for (i=1;i<=ne;i++){
         sum = 0.0;
         if (i==1){
           sum=a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny];
         }
         if (2<=i){ if(i<=ny){/* if ( (2<=i)*(i<=ny) ){}も可．以降同様 */
           sum=a2[i]*x[i-1]+a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny];}
         }
         if (ny+1 <= i){ if (i<=ne-ny){
           sum=a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny];}
         }
         if (ne-ny+1 <= i){ if (i<=ne-1){
           sum=a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]+a4[i]*x[i+1];}
         }
         if (i==ne){
           sum=a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i];
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
             double at4[NE1], double at5[NE1], double bx[NE1], double xn[NE1],
             int nx, int ny, int ne)
{
     double x1[NE1], xo[NE1], a[NE1], b[NE1], c[NE1], d[NE1], u[NE1], y[NE1];
     double eitr, omg, bnorm, rnorm, sum, zansa;
     int nitr, i, k, inx, iy, ix, iny, j;
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
       for (iy=1;iy<=ny;iy++){
         for (ix=1;ix<=nx;ix++){
           iny = iy + (ix-1)*ny;
           /* トーマス法のための係数A,B,C,Dの設定 : XNは最新のX */
           a[inx]=at1[iny]; b[inx]=at3[iny]; c[inx]=at5[iny]; d[inx]=bx[iny];
           if (iny-1 >= 1){
             d[inx]=d[inx]-at2[iny]*xn[iny-1];
           }
           if (iny+1 <= ne){
             d[inx]=d[inx]-at4[iny]*xn[iny+1];
           }
           xo[iny]=xn[iny];/* トーマス法で答えを求める前のXNをXOに保存 */
           inx=inx+1;
         }
       }
       u[1]=c[1]/b[1];/* Ly=b を解く*/
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
       for (iy=1;iy<=ny;iy++){
         for (ix=1;ix<=nx;ix++){
           iny=iy+(ix-1)*ny;
           /* 得られたX1と反復前のXOにより最新のXNを緩和 */
           xn[iny]=(1.0-omg)*xo[iny]+omg*x1[inx];
           inx=inx+1;
         }
       }
       iny = 1;/* y 軸方向への掃引 : トーマス法による */
       for (ix=1;ix<=nx;ix++){
         for (iy=1;iy<=ny;iy++){
           iny = ix + (iy-1)*nx;
           /* トーマス法のための係数A,B,C,Dの設定 : XNは最新のX */
           a[iny]=at2[iny]; b[iny]=at3[iny]; c[iny]=at4[iny]; d[iny]=bx[iny];
           if (iny-ny >= 1){
             d[iny]=d[iny]-at1[iny]*xn[iny-ny];
           }
           if (iny+ny <= ne){
             d[iny]=d[iny]-at5[iny]*xn[iny+ny];
           }
           xo[iny]=xn[iny];/* トーマス法で答えを求める前のXNをXOに保存 */
           iny=iny+1;
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
       for (ix=1;ix<=nx;ix++){
         for (iy=1;iy<=ny;iy++){
           /* 得られたX1と反復前のXOにより最新のXNを緩和 */
           xn[iny] = (1.0-omg)*xo[iny]+omg*x1[iny];
           iny = iny + 1;
         }
       }
       rnorm= 0.0;/* RNORM : 残差ベクトル(R=AX-B)の2乗ノルム */
       for (i=1;i<=ne;i++){/*  : サブルーチンSORBNDの場合分けを参照 */
         sum = 0.0;
         if (i==1){
           sum=at3[i]*xn[i]+at4[i]*xn[i+1]+at5[i]*xn[i+ny];
         }
         if (2<=i){ if(i<=ny){
           sum=at2[i]*xn[i-1]+at3[i]*xn[i]+at4[i]*xn[i+1]+at5[i]*xn[i+ny];}
         }
         if (ny+1 <= i){ if (i<=ne-ny){
           sum=at1[i]*xn[i-ny]+at2[i]*xn[i-1]+at3[i]*xn[i]+at4[i]*xn[i+1]+at5[i]*xn[i+ny];}
         }
         if (ne-ny+1 <= i){ if (i<=ne-1){
           sum=at1[i]*xn[i-ny]+at2[i]*xn[i-1]+at3[i]*xn[i]+at4[i]*xn[i+1];}
         }
         if (i==ne){
           sum=at1[i]*xn[i-ny]+at2[i]*xn[i-1]+at3[i]*xn[i];
         }
         rnorm= rnorm + (bx[i]-sum)*(bx[i]-sum);
       }
       zansa = sqrt(rnorm/bnorm);/* 残差の計算 ZANSA = || R || / || B || */
       if ( zansa <= eitr ){/* 収束判定 */
         printf (" Converged : Total ITR =  %d \n",k);
         goto label_900;
       }
     }
     printf (" Not converged ! ");/* NITRまで計算しても収束せず */
     label_900:{}/* 収束と判定されたときの分岐点 */
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
void kcrbnd (double a1[NE1], double a2[NE1], double a3[NE1],
             double a4[NE1], double a5[NE1], double b[NE1],
             double x[NE1], int ny, int ne)
{
     double r[NE1], p[NE1], ap[NE1], ar[NE1];/* 作業用配列 */
     double eitr,bnorm,rap,apap,alp,sum,zansa,arap,beta,rnorm;
     int nitr,i,itr;
     nitr = 200;/* 最大繰り返し回数の設定 */
/* 収束判定条件(変数 zansa がこの値以下になれば nitr 以下でも収束と判定する) */
     eitr = 1.0e-9;
/* 配列のゼロクリア */
     for (i=1;i<=ne;i++){
       r[i]=0.0; p[i]=0.0; ap[i]=0.0; ar[i]=0.0;
     }
     probmv (a1,a2,a3,a4,a5,x,r,ny,ne);/* 配列 r に AX を代入 */
     bnorm = 0.0;/* bのノルムの計算 */
     for (i=1;i<=ne;i++){
       bnorm = bnorm + b[i]*b[i];
     }
     for (i=1;i<=ne;i++){/*r_{0}とp_{0}(初期値)の設定*/
       r[i] = b[i] - r[i]; p[i] = r[i];/* r_{0}=b-AX, p_{0}=r_{0} */
     }
     probmv (a1,a2,a3,a4,a5,p,ap,ny,ne);/* apに A p_{0} を代入(以降のAPは式(A)で求める)*/
     itr=0;
     for (itr=1; itr<=nitr; itr++){/* 繰り返し計算 */
       rap = provv (r,ap,ne);/* ( r_{k}, A p_{k} )の計算 => RAP */
       apap = provv (ap,ap,ne);/* ( A p_{k}, A p_{k} )の計算 => APAP */
       alp = rap / apap; rnorm = 0.0;/*α_{k}の計算*/
       for (i=1;i<=ne;i++){
         x[i]=x[i]+alp*p[i]; r[i]=r[i]-alp*ap[i];
       }/* x_{k+1}=x_{k}+α_{k}p_{k},r_{k+1}=r_{k}-α_{k}Ap_{k} */
       rnorm= 0.0;/* RNORM : 残差ベクトル(R=AX-B)の2乗ノルム */
       for (i=1;i<=ne;i++){/*  : サブルーチンSORBNDの場合分けを参照 */
         sum = 0.0;
         if (i==1){
           sum=a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny];
         }
         if (2<=i){ if(i<=ny){
           sum=a2[i]*x[i-1]+a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny];}
         }
         if (ny+1 <= i){ if (i<=ne-ny){
           sum=a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny];}
         }
         if (ne-ny+1 <= i){ if (i<=ne-1){
           sum=a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]+a4[i]*x[i+1];}
         }
         if (i==ne){
           sum=a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i];
         }
         rnorm= rnorm + (b[i]-sum)*(b[i]-sum);
       }
       zansa = sqrt( rnorm / bnorm );/* 残差の計算 ZANSA = || R || / || B || */
       if (zansa <= eitr){/* 収束判定 */
         printf (" Converged : Total ITR =  %d \n",itr);
         goto label_700;
       }
       probmv (a1,a2,a3,a4,a5,r,ar,ny,ne);
       arap = provv( ar,ap,ne );/* ( A r_{k+1}, A p_{k} )の計算 => ARAP */
       beta = - arap / apap;/* β_{k}= - ( A r_{k+1}, A p_{k} ) / ( A p_{k}, A p_{k} ) */
       for (i=1;i<=ne;i++){
         p[i]=r[i]+beta*p[i]; ap[i]=ar[i]+beta*ap[i];
       }/*p_{k+1}=r_{k+1}+β_{k}p_{k},Ap_{k+1}=Ar_{k+1}+β_{k}A p_{k}<----- 式(A) */
     }
     printf ("   Not converged ! \n");/* nitrまで計算しても収束せず */
     label_700:{}/* 収束とみなされたときの分岐点 */
     for (i=1;i<=ne;i++){
       printf(" X_11_( %d ) = %15.6lE \n", i, x[i]);
     }
}
void probmv (double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1],
            double a5[NE1], double b[NE1], double c[NE1],
            int ny, int ne)
{
     int i;

     i=1; c[i]=a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny];

     for (i=2;i<=ny;i++){
       c[i]=a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny];
     }

     for (i=ny+1;i<=ne-ny;i++){
       c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny];
     }

     for (i=ne-ny+1;i<=ne-1;i++){
       c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1];
     }

     i=ne; c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i];

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
void kbibnd  (double a1[NE1], double a2[NE1], double a3[NE1],
              double a4[NE1], double a5[NE1], double b[NE1],
              double x[NE1], int ny, int ne)
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
     probmv (a1,a2,a3,a4,a5,x,r,ny,ne);/* R に AX を代入 */
     for (i=1;i<=ne;i++){/* r_{0}とp_{0}(初期値)，そして s=r_{0} の設定 */
       r[i]=b[i]-r[i]; p[i]=r[i]; s[i]=r[i];
     }
     for (j=1;j<=nitr;j++){/* 繰り返し計算 */
       sr1 = provv (s,r,NE);/* ( s, r_{k} ) の計算 => SR1 */
       probmv (a1,a2,a3,a4,a5,p,ap,ny,ne);/* A p_{k} の計算 => AP(NE) */
       sap = provv (s,ap,ne);/* ( s, A p_{k} ) の計算 => SAP */
       alpha = sr1/sap;/* α_{k} = ( s, r_{k} ) / ( s, A p_{k} ) */
       for (i=1;i<=ne;i++){
         t[i] = r[i] - alpha*ap[i];/* t_{k} = r_{k} - α_{k} A p_{k} */
       }
       probmv (a1,a2,a3,a4,a5,t,at,ny,ne);/* A t_{k} の計算 => AT(NE) */
       att = provv (at,t,ne);/* ( A t_{k}, t_{k} ) の計算 => ATT */
       atat = provv (at,at,ne);/* ( A t_{k}, A t_{k} ) の計算 => ATAT */
       xi = att/atat;/* ξ_{k} = ( A t_{k}, t_{k} ) / ( A t_{k}, A t_{k} ) */
       for (i=1;i<=ne;i++){
         x[i]=x[i]+alpha*p[i]+xi*t[i];/* x_{k+1} = x_{k} + α_{k} p_{k} + ξ_{k} t_{k} */
         r[i]=t[i]-xi*at[i];/* r_{k+1} = t_{k} - ξ_{k} A t_{k} */
       }
       rnorm=0.0;/* RNORM : 残差ベクトル(R=AX-B)の2乗ノルム */
       for (i=1;i<=ne;i++){/*  : サブルーチンSORBNDの場合分けを参照 */
         sum = 0.0;
         if (i==1){
           sum=a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny];
         }
         if (2<=i){ if(i<=ny){
           sum=a2[i]*x[i-1]+a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny];}
         }
         if (ny+1 <= i){ if (i<=ne-ny){
           sum=a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]+a4[i]*x[i+1]+a5[i]*x[i+ny];}
         }
         if (ne-ny+1 <= i){ if (i<=ne-1){
           sum=a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i]+a4[i]*x[i+1];}
         }
         if (i==ne){
           sum=a1[i]*x[i-ny]+a2[i]*x[i-1]+a3[i]*x[i];
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
