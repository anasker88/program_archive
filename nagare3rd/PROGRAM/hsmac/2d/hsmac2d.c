/********************************************************************
* ファイル名  ：hsmac2d.c                                           *
* タイトル    ：HSMAC法による2次元熱流動解析プログラム              *
* 製作者      ：平野　博之                                          *
* 所属        ：岡山理科大学 工学部 応用化学科                      *
* 製作日      ：2003.12.25                                          *
* 言語        ：C                                                   *
*********************************************************************
* プログラムを実行すると，"in2d.mac"(必ず小文字）を読みに行く．     *
* 詳細はFORTRANプログラムのコメントを参照．                         *
*                                                                   *
*  注意：Ｃ言語では配列の添字が０から始まるので，FORTRAN のX(N)を   *
*  そのままX[N]と宣言すると，X[0]からX[N-1]を意味する．ここでは     *
*  できる限りFORTRANプログラムに近くなるように，X[N+1]として添字    *
*  を1からNまで変化させてある．                                     *
********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define NX 20 /* NX:x方向格子数 */
#define NY 20 /* NY:y方向格子数 */
#define NX1 21 /* NX1=NX+1 */
#define NX2 22 /* NX2=NX+2 */
#define NY1 21 /* NY1=NY+1 */
#define NY2 22 /* NY2=NY+2 */
void cinit();
void adv();
void calvel();
void press();
void caltem();
void velbnd();
void tbnd();
void prout();
void tecplt( char *fname_tec);
double DX, DY, DT;
double VIS, ALP, BUO;
double RE, PR, GR, TIME, OMG, EPSP;
int  ICYCLE, ITR, IFLG, IRELP, METHOD;
double DMAX;
int ITYPE;
char fname[10][80];
char *fname10;
FILE *in_10, *out_11, *out_12, *out_13, *out_14, *out_21;
FILE *in_15, *in_16, *in_17, *in_18;
double UO[NX1][NY2], UN[NX1][NY2], VO[NX2][NY1], VN[NX2][NY1];
double PO[NX2][NY2], TO[NX2][NY2], TN[NX2][NY2];
void main()
{
     int i,j,NITR,NCYCLE;
     double DLX,DLY;
     char buff[80];
     fname10="in2d.mac";
     strncpy(fname[0],fname10,8);
     in_10=fopen(fname[0],"rt"); /*パラメータファイルのオープン*/
     for (i=1;i<=9;i++){/*出力ファイル名の読み込み*/
       fgets(buff, sizeof buff,in_10);
       for (j=0;buff[j]!=0;j++);
       strncpy(fname[i], buff, j-1);
     }
     out_11=fopen(fname[1],"wb");/*Uの計算結果出力用ファイルオープン*/
     out_12=fopen(fname[2],"wb");/*Vの計算結果出力用ファイルオープン*/
     out_13=fopen(fname[3],"wb");/*Pの計算結果出力用ファイルオープン*/
     out_14=fopen(fname[4],"wb");/*Tの計算結果出力用ファイルオープン*/
     fgets(buff, sizeof buff,in_10);/*in2d.mac中のコメント行(10行目)のスキップ*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(11行目)*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(12行目)*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(13行目)*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(14行目)*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(15行目)*/
     fscanf(in_10," %d %d %d %d ",&ITYPE,&ICYCLE,&NITR,&NCYCLE);
     fgets(buff, sizeof buff,in_10);/*in2d.mac中のコメント行(17行目)のスキップ*/
     fgets(buff, sizeof buff,in_10);/*in2d.mac中のコメント行(18行目)のスキップ*/
     fscanf(in_10," %lg %lg ",&EPSP,&OMG);
     fgets(buff, sizeof buff,in_10);/*in2d.mac中のコメント行(20行目)のスキップ*/
     fgets(buff, sizeof buff,in_10);/*in2d.mac中のコメント行(21行目)のスキップ*/
     fscanf(in_10," %lg %lg %lg %lg ",&DT,&RE,&PR,&GR);
     fgets(buff, sizeof buff,in_10);/*in2d.mac中のコメント行(23行目)のスキップ*/
     fgets(buff, sizeof buff,in_10);/*in2d.mac中のコメント行(24行目)のスキップ*/
     fscanf(in_10," %lg %lg %d %d",&DLX,&DLY,&IRELP,&METHOD);
     printf (" ITYPE = %d ICYCLE= %d NITR= %d NCYCLE= %d \n",
             ITYPE,ICYCLE,NITR,NCYCLE);
     printf(" EPSP= %12.3lE  OMG = %12.3lE \n",EPSP,OMG);
     printf(" DT = %12.3lE  RE = %12.3lE  PR = %12.3lE  GR = %12.3lE \n",DT,RE,PR,GR);
     printf(" DLX = %12.3lE  DLY = %12.3lE  IRELP = %d  METHOD = %d \n",DLX,DLY,IRELP,METHOD);
     if ( ICYCLE != 0 ){/*継続の計算の場合*/
       in_15=fopen(fname[5],"rb");/*Uデータファイルのオープン*/
       in_16=fopen(fname[6],"rb");/*Vデータファイルのオープン*/
       in_17=fopen(fname[7],"rb");/*Pデータファイルのオープン*/
       in_18=fopen(fname[8],"rb");/*Tデータファイルのオープン*/
     }
     DX=DLX/(double)NX; DY=DLY/(double)NY;/*x,y方向の格子幅DX,DY*/
     VIS=PR; ALP=1.0; BUO=( GR * ( PR*PR ) );
      /*VIS:運動方程式中の拡散項の係数(ここではPr)*/
      /*ALP:エネルギー方程式中の拡散項の係数(ここでは1)*/
      /*BUO:浮力項の係数(ここでは Gr * Pr**2)*/
     if ( ITYPE==1 ){/* 等温場なら浮力項の係数をゼロに設定 */
       BUO = 0.0;
     }
     cinit();/*初期値の設定*/
     label_700:{};/*時間進行のための戻り点*/
     adv();/*時間進行*/
     calvel();/*速度場の計算*/
     ITR=1;
     label_710:{};/*圧力反復のための戻り点*/
     IFLG=0;
     press();/*仮の圧力場の計算*/
     if ( IFLG==0 ){/*Newton-Raphson法による圧力場の計算が収束したとき*/
       if ( ITYPE==2 ){/*非等温場計算の場合*/
         caltem();/*温度場を計算*/
       }
     }
     if ( IFLG==1 ){/*圧力場の計算が収束していないとき*/
       if ( ITR<NITR ){/*圧力計算の反復回数があらかじめ設定された最大値NITRより小さいとき*/
         ITR = ITR + 1;/*さらに反復を繰り返す*/
         goto label_710;
       }
       else{/*圧力計算の反復回数がNITR以上のとき発散とみなして計算終了*/
         printf(" calculation has diverged \n");
         prout();
         goto label_900;
       }
     }
     if ( ICYCLE < NCYCLE ){/*時間進行サイクル(ICYCLE)がNCYCLEより小さい時*/
       goto label_700;
     }
     else{/*時間進行サイクルがNCYCLEになったら->計算終了*/
       prout();
     }
     label_900:{};
     tecplt(fname[9]);/*Tecplot用データの出力*/
     fclose (in_10); fclose (out_11); fclose (out_12);
     fclose (out_13); fclose (out_14);
}
/*初期設定*/
void cinit()
{
     int ix,iy;
     if ( ICYCLE == 0 ){/*新規計算の場合*/
       for (ix=0;ix<=NX;ix++){/*Uの初期値設定*/
         for (iy=0;iy<=NY+1;iy++){
           UN[ix][iy] = 0.0;
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*Vの初期値設定*/
         for (iy=0;iy<=NY;iy++){
           VN[ix][iy] = 0.0;
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*Pの初期値設定*/
         for (iy=0;iy<=NY+1;iy++){
           PO[ix][iy] = 0.0;
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*Tの初期値設定(領域内は高温(+0.5)と低温(-0.5)の中間温度)*/
         for (iy=0;iy<=NY+1;iy++){/*（注意）浮力項の計算で温度の配列を使用して*/
           TN[ix][iy] = 0.0;/*いるので等温場でもT=0として初期条件だけは設定する*/
         }/*必要がある．ゼロ以外の値を入れると浮力項が計算される可能性があるので注意．*/
       }
       for (iy=0;iy<=NY+1;iy++){/*Tの境界：右側壁（冷却）T=-0.5*/
         TN[NX+1][iy] = 2.0 * ( -0.5 ) - TN[NX][iy];
       }
       for (iy=0;iy<=NY+1;iy++){/*Tの境界：左側壁（加熱）T=+0.5*/
         TN[0][iy] = 2.0 * ( +0.5 ) - TN[1][iy];
       }
       for (ix=1;ix<=NX;ix++){/*Tの境界：上面（断熱）*/
         TN[ix][NY+1] = TN[ix][NY];
       }
       for (ix=1;ix<=NX;ix++){/*Tの境界：下面（断熱）*/
         TN[ix][0] = TN[ix][1];
       }
     }
     else{/*継続計算（すでにある計算結果からスタート）の場合*/
       /*Uデータファイルからの読み込み[Unit No.=15]*/
       fread(UN, sizeof(double), NX1*NY2, in_15);
       /*Vデータファイルからの読み込み[Unit No.=16]*/
       fread(VN, sizeof(double), NX2*NY1, in_16);
       /*Pデータファイルからの読み込み[Unit No.=17]*/
       fread(PO, sizeof(double), NX2*NY2, in_17);
       /*Tデータファイルからの読み込み[Unit No.=18]*/
       fread(TN, sizeof(double), NX2*NY2, in_18);
       fclose (in_15); fclose (in_16); fclose (in_17); fclose (in_18);
     }
}
/*時間進行*/
void adv()
{
     int ix,iy;
     TIME = DT*(double)ICYCLE; ICYCLE = ICYCLE + 1;
     if ( (ICYCLE%100) ==0 ){/*ICYCLEを100回毎に表示*/
       printf ("  CYC = %8d \n",ICYCLE);
     }
     for (ix=0;ix<=NX;ix++){/* UN -> UO : 必要なら入れ替える前にUNとUOから変動量を求める*/
       for (iy=0;iy<=NY+1;iy++){/*UN:前の時間ステップにおいて最終的に得られた値，圧力補正の度に更新される*/
         UO[ix][iy] = UN[ix][iy];/*UO:新しい時間ステップでの初期値．UNを保存．*/
       }
     }
     for (ix=0;ix<=NX+1;ix++){/* VN -> VO : 必要なら入れ替える前にVNとVOから変動量を求める*/
       for (iy=0;iy<=NY;iy++){/*VN:前の時間ステップにおいて最終的に得られた値，圧力補正の度に更新される*/
         VO[ix][iy] = VN[ix][iy];/*VO:新しい時間ステップでの初期値．VNを保存．*/
       }
     }
     for (ix=0;ix<=NX+1;ix++){/* TN -> TO : 必要なら入れ替える前にTNとTOから変動量を求める*/
       for (iy=0;iy<=NY+1;iy++){/*TN:前の時間ステップでの計算値*/
         TO[ix][iy] = TN[ix][iy];/*TO:新しい時間ステップでの初期値．TNを保存．*/
       }
     }
}
/*速度場の計算*/
void calvel()
{
     int ix,iy;
     double vv,cnvux,cnvuy,tu,buou,difu;
     double uu,cnvvx,cnvvy,tv,buov,difv;
     /*U(ix,iy)の計算*/
     for (ix=1;ix<=NX-1;ix++){
       for (iy=1;iy<=NY;iy++){
         /* vvはU(ix,iy)におけるVの補間値 */
         vv=(VO[ix][iy]+VO[ix+1][iy]+VO[ix][iy-1]+VO[ix+1][iy-1])/4.0;
         /* 対流項(cnvux,cnvuy)を１次精度風上差分にて計算 */
         if ( UO[ix][iy] >= 0.0 ){
           cnvux = UO[ix][iy]*( UO[ix][iy]-UO[ix-1][iy] ) / DX;
         }
         else if ( UO[ix][iy] < 0.0 ){
           cnvux = UO[ix][iy]*( UO[ix+1][iy] - UO[ix][iy] ) / DX;
         }
         if ( vv >= 0.0 ){
           cnvuy = vv*( UO[ix][iy] - UO[ix][iy-1] ) / DY;
         }
         else if ( vv < 0.0 ){
           cnvuy = vv*( UO[ix][iy+1] - UO[ix][iy] ) / DY;
         }
         /* x方向の浮力項(buou)はゼロ */
         tu = 0.0;
         buou = 0.0;
         /* 拡散項(difu)の計算 */
         difu = VIS*( (UO[ix-1][iy]-2.0*UO[ix][iy]+UO[ix+1][iy])/(DX*DX)
                     +(UO[ix][iy-1]-2.0*UO[ix][iy]+UO[ix][iy+1])/(DY*DY) );
         /*仮の速度(U)の計算*/
         UN[ix][iy] = UO[ix][iy]
              + DT*( -cnvux-cnvuy+difu+buou+( PO[ix][iy]-PO[ix+1][iy] )/DX );
       }
     }
     /*V(ix,iy)の計算*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY-1;iy++){
         /* uuはV(ix,iy)におけるUの補間値 */
         uu=(UO[ix-1][iy]+UO[ix][iy]+UO[ix-1][iy+1]+UO[ix][iy+1])/4.0;
         /* 対流項(cnvvx,cnvvy)を１次精度風上差分にて計算 */
         if ( uu >= 0.0 ){
           cnvvx = uu*( VO[ix][iy] - VO[ix-1][iy] ) / DX;
         }
         else if ( uu < 0.0 ){
           cnvvx = uu*( VO[ix+1][iy] - VO[ix][iy] ) / DX;
         }
         if ( VO[ix][iy] >= 0.0 ){
           cnvvy = VO[ix][iy]*( VO[ix][iy]-VO[ix][iy-1] ) / DY;
         }
         else if ( VO[ix][iy] < 0.0 ){
           cnvvy = VO[ix][iy]*( VO[ix][iy+1]-VO[ix][iy] ) / DY;
         }
         /*浮力項(buov)の計算*/
         tv = ( TO[ix][iy] + TO[ix][iy+1] )/2.0;
         buov = BUO*tv;
         /*拡散項(difv)の計算*/
         difv = VIS*( (VO[ix-1][iy]-2.0*VO[ix][iy]+VO[ix+1][iy])/(DX*DX)
                     +(VO[ix][iy-1]-2.0*VO[ix][iy]+VO[ix][iy+1])/(DY*DY) );
         /*仮の速度(V)の計算*/
         VN[ix][iy] = VO[ix][iy]
              + DT*( -cnvvx-cnvvy+difv+buov+(PO[ix][iy]-PO[ix][iy+1])/DY );
       }
     }
     velbnd();
}
/*圧力場の計算*/
void press()
{
     int ixmax,iymax,ix,iy;
     double del,div,delp,postn;
     ixmax = 0; iymax = 0; DMAX = 0.0e0;
     /*P(ix,iy)の計算*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         del = DT*( 2.0/(DX*DX) + 2.0/(DY*DY) );
         div = ( UN[ix][iy] - UN[ix-1][iy  ] )/DX
             + ( VN[ix][iy] - VN[ix  ][iy-1] )/DY;
         if ( fabs(div) >= fabs(DMAX) ){
           ixmax = ix; iymax = iy; DMAX = div;
         }
         delp = - OMG * div / del;
         PO[ix][iy] = PO[ix][iy] + delp;
         UN[ix  ][iy  ]=UN[ix  ][iy  ]+DT/DX*delp;
         UN[ix-1][iy  ]=UN[ix-1][iy  ]-DT/DX*delp;
         VN[ix  ][iy  ]=VN[ix  ][iy  ]+DT/DY*delp;
         VN[ix  ][iy-1]=VN[ix  ][iy-1]-DT/DY*delp;
       }
     }
     /* 圧力の相対性に関する処理(IRELP=1なら以下の処理を行う) */
     if (IRELP==1){
       postn = PO[1][1];
       for (ix=1;ix<=NX;ix++){
         for (iy=1;iy<=NY;iy++){
           PO[ix][iy] = PO[ix][iy]-postn;
         }
       }
     }
     /* IFLG=1なら，連続の式を満たしていないと判定し再び圧力計算を行う */
     if ( fabs(DMAX) >= EPSP ){
       IFLG = 1;
     }
     /* 圧力計算の回数を100回ごとに表示 */
     if ( (ITR%100) == 0){
       printf (" Iteration= %8d, Div(max)( %6d, %6d ) = %15.6lE\n",ITR,ixmax,iymax,DMAX);
     }
     /*新たに得られた速度を用いて境界条件を処理する*/
     velbnd();
}
/*温度場の計算*/
void caltem()
{
     int ix,iy;
     double uut,vvt,cnvtx,cnvty,dift;
     /*T(ix,iy)の計算*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         /* uut,vvtはそれぞT(IX,IY)におけるU,Vの補間値 */
         uut = ( UO[ix][iy] + UO[ix-1][iy  ] ) / 2.0;
         vvt = ( VO[ix][iy] + VO[ix  ][iy-1] ) / 2.0;
         /*対流項(cnvtx,cnvty)を１次精度風上差分にて計算*/
         if ( uut >= 0.0 ){
           cnvtx = uut*( TO[ix][iy] - TO[ix-1][iy] ) / DX;
         }
         else if ( uut < 0.0 ){
           cnvtx = uut*( TO[ix+1][iy] - TO[ix][iy] ) / DX;
         }
         if ( vvt >= 0.0 ){
           cnvty = vvt*( TO[ix][iy] - TO[ix][iy-1] ) / DY;
         }
         else if ( vvt < 0.0 ){
           cnvty = vvt*( TO[ix][iy+1] - TO[ix][iy] ) / DY;
         }
         /* 拡散項(dift)の計算 */
         dift = ALP*( (TO[ix-1][iy]-2.0*TO[ix][iy]+TO[ix+1][iy])/(DX*DX)
                     +(TO[ix][iy-1]-2.0*TO[ix][iy]+TO[ix][iy+1])/(DY*DY) );
         /* 次の時間のTの計算 */
         TN[ix][iy] = TO[ix][iy] + DT*( -cnvtx-cnvty+dift );
       }
     }
     /* 境界条件の処理 */
     tbnd();
}
/* 速度の境界条件の処理 */
void velbnd()
{
     int ix,iy;
     for (iy=1;iy<=NY;iy++){/*U（右側面）*/
       UN[NX][iy] = 0.0;
     }
     for (iy=1;iy<=NY;iy++){/*U（左側面）*/
       UN[0][iy] = 0.0;
     }
     for (ix=0;ix<=NX;ix++){/*U（上面）*/
       UN[ix][NY+1] = -UN[ix][NY];
     }
     for (ix=0;ix<=NX;ix++){/*U（下面）*/
       UN[ix][0] = -UN[ix][1];
     }
     for (iy=1;iy<=NY-1;iy++){/*V（右側面）*/
       VN[NX+1][iy] = -VN[NX][iy];
     }
     for (iy=1;iy<=NY-1;iy++){/*V（左側面）*/
       VN[0][iy] = -VN[1][iy];
     }
     for (ix=0;ix<=NX+1;ix++){/*V（上面）*/
       VN[ix][NY] = 0.0;
     }
     for (ix=0;ix<=NX+1;ix++){/*V（下面）*/
       VN[ix][0] = 0.0;
     }
}
/*温度の境界条件の処理*/
void tbnd()
{
     int ix,iy;
      for (iy=0;iy<=NY+1;iy++){/*右側面*/
        TN[NX+1][iy] = 2.0 * ( -0.5 ) - TN[NX][iy];
      }
      for (iy=0;iy<=NY+1;iy++){/*左側面*/
        TN[0][iy] = 2.0 * ( +0.5 ) - TN[1][iy];
      }
      for (ix=1;ix<=NX;ix++){/*上面*/
        TN[ix][NY+1] = TN[ix][NY];
      }
      for (ix=1;ix<=NX;ix++){/*下面*/
        TN[ix][0] = TN[ix][1];
      }
}
void prout()/*データ出力用*/
{
    fwrite(UN, sizeof(double), NX1*NY2, out_11);
    fwrite(VN, sizeof(double), NX2*NY1, out_12);
    fwrite(PO, sizeof(double), NX2*NY2, out_13);
    fwrite(TN, sizeof(double), NX2*NY2, out_14);
}
/* Tecplot用データの出力 */
void tecplt ( char *fname_tec )
{

     int ix,iy;
     double x,y,u,v,t;

     out_21=fopen(fname_tec,"wt");
     fprintf (out_21," VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"T\"\n");
     fprintf (out_21," ZONE I=%3d, J=%3d, F=POINT \n",NX1,NY1);

     for (iy=0;iy<=NY;iy++){
       for (ix=0;ix<=NX;ix++){
         x = DX * (double)ix;
         y = DY * (double)iy;
         u = ( UN[ix][iy]+UN[ix][iy+1])/2.0;
         v = ( VN[ix][iy]+VN[ix+1][iy])/2.0;
         t = ( TN[ix][iy]+TN[ix+1][iy]
              +TN[ix][iy+1]+TN[ix+1][iy+1] )/4.0;
         fprintf(out_21,"%12.3lE %12.3lE %12.3lE %12.3lE %12.3lE\n",x,y,u,v,t);
       }
     }
     fclose (out_21);
}
