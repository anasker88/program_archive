/********************************************************************
* ファイル名  ：hsmac3d.c                                           *
* タイトル    ：HSMAC法による3次元熱流動解析プログラム              *
* 製作者      ：平野　博之                                          *
* 所属        ：岡山理科大学 工学部 応用化学科                      *
* 製作日      ：2003.12.25                                          *
* 言語        ：C                                                   *
*********************************************************************
* プログラムを実行すると，"in3d.mac"(必ず小文字）を読みに行く．     *
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
#define NZ 20 /* NZ:z方向格子数 */
#define NX1 21 /* NX1=NX+1 */
#define NX2 22 /* NX2=NX+2 */
#define NY1 21 /* NY1=NY+1 */
#define NY2 22 /* NY2=NY+2 */
#define NZ1 21 /* NZ1=NZ+1 */
#define NZ2 22 /* NZ2=NZ+2 */
void cinit();
void adv();
void calvel();
void press();
void caltem();
void velbnd();
void tbnd();
void prout();
void tecplt( char *fname_tec);
double DX, DY, DZ, DT;
double VIS, ALP, BUO;
double RE, PR, GR, TIME, OMG, EPSP;
int  ICYCLE, ITR, IFLG, IRELP, METHOD;
double DMAX;
int ITYPE;
char fname[12][80];
char *fname10;
FILE *in_10, *out_11, *out_12, *out_13, *out_14, *out_15, *out_21;
FILE *in_16, *in_17, *in_18, *in_19, *in_20;
double UO[NX1][NY2][NZ2],UN[NX1][NY2][NZ2],VO[NX2][NY1][NZ2],VN[NX2][NY1][NZ2];
double WO[NX2][NY2][NZ1],WN[NX2][NY2][NZ1],PO[NX2][NY2][NZ2],TO[NX2][NY2][NZ2],TN[NX2][NY2][NZ2];
void main()
{
     int i,j,NITR,NCYCLE;
     double DLX,DLY,DLZ;
     char buff[80];
     fname10="in3d.mac";
     strncpy(fname[0],fname10,8);
     in_10=fopen(fname[0],"rt"); /*パラメータファイルのオープン*/
     for (i=1;i<=11;i++){/*出力ファイル名の読み込み*/
       fgets(buff, sizeof buff,in_10);
       for (j=0;buff[j]!=0;j++);
       strncpy(fname[i], buff, j-1);
     }
     out_11=fopen(fname[1],"wb");/*Uの計算結果出力用ファイルオープン*/
     out_12=fopen(fname[2],"wb");/*Vの計算結果出力用ファイルオープン*/
     out_13=fopen(fname[3],"wb");/*Wの計算結果出力用ファイルオープン*/
     out_14=fopen(fname[4],"wb");/*Pの計算結果出力用ファイルオープン*/
     out_15=fopen(fname[5],"wb");/*Tの計算結果出力用ファイルオープン*/
     fgets(buff, sizeof buff,in_10);/*in3d.mac中のコメント行(12行目)のスキップ*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(13行目)*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(14行目)*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(15行目)*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(16行目)*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(17行目)*/
     fscanf(in_10," %d %d %d %d ",&ITYPE,&ICYCLE,&NITR,&NCYCLE);
     fgets(buff, sizeof buff,in_10);/*in3d.mac中のコメント行(19行目)のスキップ*/
     fgets(buff, sizeof buff,in_10);/*in3d.mac中のコメント行(20行目)のスキップ*/
     fscanf(in_10," %lg %lg ",&EPSP,&OMG);
     fgets(buff, sizeof buff,in_10);/*in3d.mac中のコメント行(22行目)のスキップ*/
     fgets(buff, sizeof buff,in_10);/*in3d.mac中のコメント行(23行目)のスキップ*/
     fscanf(in_10," %lg %lg %lg %lg ",&DT,&RE,&PR,&GR);
     fgets(buff, sizeof buff,in_10);/*in3d.mac中のコメント行(25行目)のスキップ*/
     fgets(buff, sizeof buff,in_10);/*in3d.mac中のコメント行(26行目)のスキップ*/
     fscanf(in_10," %lg %lg %lg %d %d",&DLX,&DLY,&DLZ,&IRELP,&METHOD);
     printf (" ITYPE = %d ICYCLE= %d NITR= %d NCYCLE= %d \n",
             ITYPE,ICYCLE,NITR,NCYCLE);
     printf(" EPSP= %12.3lE  OMG = %12.3lE \n",EPSP,OMG);
     printf(" DT = %12.3lE  RE = %12.3lE  PR = %12.3lE  GR = %12.3lE \n",DT,RE,PR,GR);
     printf(" DLX = %12.3lE  DLY = %12.3lE  DLZ = %12.3lE  IRELP = %d  METHOD = %d \n",
              DLX,DLY,DLZ,IRELP,METHOD);
     if ( ICYCLE != 0 ){/*継続の計算の場合*/
       in_16=fopen(fname[6],"rb");/*Uデータファイルのオープン*/
       in_17=fopen(fname[7],"rb");/*Vデータファイルのオープン*/
       in_18=fopen(fname[8],"rb");/*Wデータファイルのオープン*/
       in_19=fopen(fname[9],"rb");/*Pデータファイルのオープン*/
       in_20=fopen(fname[10],"rb");/*Tデータファイルのオープン*/
     }
     DX=DLX/(double)NX; DY=DLY/(double)NY; DZ=DLZ/(double)NZ;/*x,y,z方向の格子幅DX,DY,DZ*/
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
     tecplt(fname[11]);/*Tecplot用データの出力*/
     fclose (in_10);  fclose (out_11); fclose (out_12);
     fclose (out_13); fclose (out_14); fclose(out_15);
}
/*初期設定*/
void cinit()
{
     int ix,iy,iz;
     if ( ICYCLE == 0 ){/*新規計算の場合*/
       for (ix=0;ix<=NX;ix++){/*Uの初期値設定*/
         for (iy=0;iy<=NY+1;iy++){
           for (iz=0;iz<=NZ+1;iz++){
             UN[ix][iy][iz] = 0.0;
           }
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*Vの初期値設定*/
         for (iy=0;iy<=NY;iy++){
           for (iz=0;iz<=NZ+1;iz++){
             VN[ix][iy][iz] = 0.0;
           }
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*Wの初期値設定*/
         for (iy=0;iy<=NY+1;iy++){
           for (iz=0;iz<=NZ;iz++){
             WN[ix][iy][iz] = 0.0;
           }
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*Pの初期値設定*/
         for (iy=0;iy<=NY+1;iy++){
           for (iz=0;iz<=NZ+1;iz++){
             PO[ix][iy][iz] = 0.0;
           }
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*Tの初期値設定(領域内は高温(+0.5)と低温(-0.5)の中間温度)  */
         for (iy=0;iy<=NY+1;iy++){/* （注意）浮力項の計算で温度の配列を使用しているので，  */
           for (iz=0;iz<=NZ+1;iz++){/*       等温場でもT=0として初期条件だけは設定する必要 */
             TN[ix][iy][iz] = 0.0;/*         がある．ゼロ以外の値を入れると浮力項が計算され*/
           }/*                               る可能性があるので注意．*/
         }
       }
       for (iy=0;iy<=NY+1;iy++){/*Tの境界：右側壁（冷却）T=-0.5*/
         for (iz=0;iz<=NZ+1;iz++){
           TN[NX+1][iy][iz] = 2.0 * ( -0.5 ) - TN[NX][iy][iz];
         }
       }
       for (iy=0;iy<=NY+1;iy++){/*Tの境界：左側壁（加熱）T=+0.5*/
         for (iz=0;iz<=NZ+1;iz++){
           TN[0][iy][iz] = 2.0 * ( +0.5 ) - TN[1][iy][iz];
         }
       }
       for (ix=1;ix<=NX;ix++){/*Tの境界：上面（断熱）*/
         for (iz=0;iz<=NZ+1;iz++){
           TN[ix][NY+1][iz] = TN[ix][NY][iz];
         }
       }
       for (ix=1;ix<=NX;ix++){/*Tの境界：下面（断熱）*/
         for (iz=0;iz<=NZ+1;iz++){
           TN[ix][0][iz] = TN[ix][1][iz];
         }
       }
       for (ix=1;ix<=NX;ix++){/*Tの境界：後面（断熱）*/
         for (iy=1;iy<=NY;iy++){
           TN[ix][iy][NZ+1] = TN[ix][iy][NZ];
         }
       }
       for (ix=1;ix<=NX;ix++){/*Tの境界：前面（断熱）*/
         for (iy=1;iy<=NY;iy++){
           TN[ix][iy][0] = TN[ix][iy][1];
         }
       }
     }
     else{/*継続計算（すでにある計算結果からスタート）の場合*/
       /*Uデータファイルからの読み込み[Unit No.=16]*/
       fread(UN, sizeof(double), NX1*NY2*NZ2, in_16);
       /*Vデータファイルからの読み込み[Unit No.=17]*/
       fread(VN, sizeof(double), NX2*NY1*NZ2, in_17);
       /*Wデータファイルからの読み込み[Unit No.=18]*/
       fread(WN, sizeof(double), NX2*NY2*NZ1, in_18);
       /*Pデータファイルからの読み込み[Unit No.=19]*/
       fread(PO, sizeof(double), NX2*NY2*NZ2, in_19);
       /*Tデータファイルからの読み込み[Unit No.=20]*/
       fread(TN, sizeof(double), NX2*NY2*NZ2, in_20);
       fclose (in_16); fclose (in_17); fclose (in_18); fclose (in_19); fclose (in_20);
     }
}
/*時間進行*/
void adv()
{
     int ix,iy,iz;
     TIME = DT*(double)ICYCLE; ICYCLE = ICYCLE + 1;
     if ( (ICYCLE%100) ==0 ){/*ICYCLEを100回毎に表示*/
       printf ("  CYC = %8d \n",ICYCLE);
     }
     for (ix=0;ix<=NX;ix++){/* UN -> UO : 必要なら入れ替える前にUNとUOから変動量を求める*/
       for (iy=0;iy<=NY+1;iy++){/* UN:前の時間ステップにおいて最終的に得られた値，圧力補正の度に更新される*/
         for (iz=0;iz<=NZ+1;iz++){/*UO:新しい時間ステップでの初期値．UNを保存．*/
           UO[ix][iy][iz] = UN[ix][iy][iz];
         }
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*    VN -> VO : 必要なら入れ替える前にVNとVOから変動量を求める*/
       for (iy=0;iy<=NY;iy++){/*    VN:前の時間ステップにおいて最終的に得られた値，圧力補正の度に更新される*/
         for (iz=0;iz<=NZ+1;iz++){/*VO:新しい時間ステップでの初期値．VNを保存．*/
           VO[ix][iy][iz] = VN[ix][iy][iz];
         }
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*WN -> WO : 必要なら入れ替える前にWNとWOから変動量を求める*/
       for (iy=0;iy<=NY+1;iy++){/*WN:前の時間ステップにおいて最終的に得られた値，圧力補正の度に更新される*/
         for (iz=0;iz<=NZ;iz++){/*WO:新しい時間ステップでの初期値．WNを保存．*/
           WO[ix][iy][iz] = WN[ix][iy][iz];
         }
       }
     }
     for (ix=0;ix<=NX+1;ix++){/* TN -> TO : 必要なら入れ替える前にTNとTOから変動量を求める*/
       for (iy=0;iy<=NY+1;iy++){/*  TN:前の時間ステップでの計算値*/
         for (iz=0;iz<=NZ+1;iz++){/*TO:新しい時間ステップでの初期値．TNを保存．*/
           TO[ix][iy][iz] = TN[ix][iy][iz];
         }
       }
     }
}
/*速度場の計算*/
void calvel()
{
     int ix,iy,iz;
     double vv,ww,cnvux,cnvuy,cnvuz,tu,buou,difu;
     double uu,cnvvx,cnvvy,cnvvz,tv,buov,difv;
     double cnvwx,cnvwy,cnvwz,tw,buow,difw;
     /*U(ix,iy)の計算*/
     for (ix=1;ix<=NX-1;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
         /* vvはU(ix,iy,iz)におけるVの補間値 */
         vv=(VO[ix][iy][iz]+VO[ix+1][iy][iz]+VO[ix][iy-1][iz]+VO[ix+1][iy-1][iz])/4.0;
         /* wwはU(ix,iy,iz)におけるWの補間値 */
         ww=(WO[ix][iy][iz]+WO[ix+1][iy][iz]+WO[ix][iy][iz-1]+WO[ix+1][iy][iz-1])/4.0;
         /* 対流項(cnvux,cnvuy)を１次精度風上差分にて計算 */
         if ( UO[ix][iy][iz] >= 0.0 ){
           cnvux = UO[ix][iy][iz]*( UO[ix][iy][iz]-UO[ix-1][iy][iz] ) / DX;
         }
         else if ( UO[ix][iy][iz] < 0.0 ){
           cnvux = UO[ix][iy][iz]*( UO[ix+1][iy][iz] - UO[ix][iy][iz] ) / DX;
         }
         if ( vv >= 0.0 ){
           cnvuy = vv*( UO[ix][iy][iz] - UO[ix][iy-1][iz] ) / DY;
         }
         else if ( vv < 0.0 ){
           cnvuy = vv*( UO[ix][iy+1][iz] - UO[ix][iy][iz] ) / DY;
         }
         if ( ww >= 0.0 ){
           cnvuz = ww*( UO[ix][iy][iz] - UO[ix][iy][iz-1] ) / DZ;
         }
         else if ( ww < 0.0 ){
           cnvuz = ww*( UO[ix][iy][iz+1] - UO[ix][iy][iz] ) / DZ;
         }
         /* x方向の浮力項(buou)はゼロ */
         tu = 0.0;
         buou = BUO*tu;
         /* 拡散項(difu)の計算 */
         difu = VIS*( (UO[ix-1][iy][iz]-2.0*UO[ix][iy][iz]+UO[ix+1][iy][iz])/(DX*DX)
                     +(UO[ix][iy-1][iz]-2.0*UO[ix][iy][iz]+UO[ix][iy+1][iz])/(DY*DY)
                     +(UO[ix][iy][iz-1]-2.0*UO[ix][iy][iz]+UO[ix][iy][iz+1])/(DZ*DZ) );
         /*仮の速度(U)の計算*/
         UN[ix][iy][iz] = UO[ix][iy][iz]
              + DT*( -cnvux-cnvuy-cnvuz+difu+buou+( PO[ix][iy][iz]-PO[ix+1][iy][iz] )/DX );
         }
       }
     }
     /*V(ix,iy,iz)の計算*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY-1;iy++){
         for (iz=1;iz<=NZ;iz++){
         /* uuはV(ix,iy,iz)におけるUの補間値 */
         uu=(UO[ix-1][iy][iz]+UO[ix][iy][iz]+UO[ix-1][iy+1][iz]+UO[ix][iy+1][iz])/4.0;
         /* wwはV(ix,iy,iz)におけるWの補間値 */
         ww=(WO[ix][iy][iz-1]+WO[ix][iy][iz]+WO[ix][iy+1][iz-1]+WO[ix][iy+1][iz])/4.0;
         /* 対流項(cnvvx,cnvvy,cnvvz)を１次精度風上差分にて計算 */
         if ( uu >= 0.0 ){
           cnvvx = uu*( VO[ix][iy][iz] - VO[ix-1][iy][iz] ) / DX;
         }
         else if ( uu < 0.0 ){
           cnvvx = uu*( VO[ix+1][iy][iz] - VO[ix][iy][iz] ) / DX;
         }
         if ( VO[ix][iy][iz] >= 0.0 ){
           cnvvy = VO[ix][iy][iz]*( VO[ix][iy][iz]-VO[ix][iy-1][iz] ) / DY;
         }
         else if ( VO[ix][iy][iz] < 0.0 ){
           cnvvy = VO[ix][iy][iz]*( VO[ix][iy+1][iz]-VO[ix][iy][iz] ) / DY;
         }
         if ( ww >= 0.0 ){
           cnvvz = WO[ix][iy][iz]*( VO[ix][iy][iz]-VO[ix][iy][iz-1] ) / DZ;
         }
         else if ( ww < 0.0 ){
           cnvvz = WO[ix][iy][iz]*( VO[ix][iy][iz+1]-VO[ix][iy][iz] ) / DZ;
         }
         /*浮力項(buov)の計算*/
         tv = ( TO[ix][iy][iz] + TO[ix][iy+1][iz] )/2.0;
         buov = BUO*tv;
         /*拡散項(difv)の計算*/
         difv = VIS*( (VO[ix-1][iy][iz]-2.0*VO[ix][iy][iz]+VO[ix+1][iy][iz])/(DX*DX)
                     +(VO[ix][iy-1][iz]-2.0*VO[ix][iy][iz]+VO[ix][iy+1][iz])/(DY*DY)
                     +(VO[ix][iy][iz-1]-2.0*VO[ix][iy][iz]+VO[ix][iy][iz+1])/(DZ*DZ) );
         /*仮の速度(V)の計算*/
         VN[ix][iy][iz] = VO[ix][iy][iz]
              + DT*( -cnvvx-cnvvy-cnvvz+difv+buov+(PO[ix][iy][iz]-PO[ix][iy+1][iz])/DY );
         }
       }
     }
     /*W(ix,iy,iz)の計算*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ-1;iz++){
         /* uuはW(ix,iy,iz)におけるUの補間値 */
         uu=(UO[ix-1][iy][iz]+UO[ix][iy][iz]+UO[ix-1][iy][iz+1]+UO[ix][iy][iz+1])/4.0;
         /* vvはW(ix,iy,iz)におけるVの補間値 */
         vv=(VO[ix][iy-1][iz]+VO[ix][iy][iz]+VO[ix][iy-1][iz+1]+VO[ix][iy][iz+1])/4.0;
         /* 対流項(cnvwx,cnvwy,cnvwz)を１次精度風上差分にて計算 */
         if ( uu >= 0.0 ){
           cnvwx = uu*( WO[ix][iy][iz] - WO[ix-1][iy][iz] ) / DX;
         }
         else if ( uu < 0.0 ){
           cnvwx = uu*( WO[ix+1][iy][iz] - WO[ix][iy][iz] ) / DX;
         }
         if ( vv >= 0.0 ){
           cnvwy = vv*( WO[ix][iy][iz]-WO[ix][iy-1][iz] ) / DY;
         }
         else if ( vv < 0.0 ){
           cnvwy = vv*( WO[ix][iy+1][iz]-WO[ix][iy][iz] ) / DY;
         }
         if ( WO[ix][iy][iz] >= 0.0 ){
           cnvwz = WO[ix][iy][iz]*( WO[ix][iy][iz]-WO[ix][iy][iz-1] ) / DZ;
         }
         else if ( WO[ix][iy][iz] < 0.0 ){
           cnvwz = WO[ix][iy][iz]*( WO[ix][iy][iz+1]-WO[ix][iy][iz] ) / DZ;
         }
         /*浮力項(buow)の計算*/
         tw = 0.0;
         buow = BUO*tw;
         /*拡散項(difw)の計算*/
         difw = VIS*( (WO[ix-1][iy][iz]-2.0*WO[ix][iy][iz]+WO[ix+1][iy][iz])/(DX*DX)
                     +(WO[ix][iy-1][iz]-2.0*WO[ix][iy][iz]+WO[ix][iy+1][iz])/(DY*DY)
                     +(WO[ix][iy][iz-1]-2.0*WO[ix][iy][iz]+WO[ix][iy][iz+1])/(DZ*DZ) );
         /*仮の速度(W)の計算*/
         WN[ix][iy][iz] = WO[ix][iy][iz]
              + DT*( -cnvwx-cnvwy-cnvwz+difw+buow+(PO[ix][iy][iz]-PO[ix][iy][iz+1])/DZ );
         }
       }
     }
     velbnd();
}
/*圧力場の計算*/
void press()
{
     int ixmax,iymax,izmax,ix,iy,iz;
     double del,div,delp,postn;
     ixmax = 0; iymax = 0; izmax=0; DMAX = 0.0e0;
     /*P(ix,iy,iz)の計算*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
           del = DT*( 2.0/(DX*DX) + 2.0/(DY*DY) + 2.0/(DZ*DZ) );
           div = ( UN[ix][iy][iz] - UN[ix-1][iy][iz] )/DX
               + ( VN[ix][iy][iz] - VN[ix][iy-1][iz] )/DY
               + ( WN[ix][iy][iz] - WN[ix][iy][iz-1] )/DZ;
           if ( fabs(div) >= fabs(DMAX) ){
             ixmax = ix; iymax = iy; izmax = iz; DMAX = div;
           }
           delp = - OMG * div / del;
           PO[ix][iy][iz] = PO[ix][iy][iz] + delp;
           UN[ix  ][iy  ][iz  ]=UN[ix  ][iy  ][iz  ]+DT/DX*delp;
           UN[ix-1][iy  ][iz  ]=UN[ix-1][iy  ][iz  ]-DT/DX*delp;
           VN[ix  ][iy  ][iz  ]=VN[ix  ][iy  ][iz  ]+DT/DY*delp;
           VN[ix  ][iy-1][iz  ]=VN[ix  ][iy-1][iz  ]-DT/DY*delp;
           WN[ix  ][iy  ][iz  ]=WN[ix  ][iy  ][iz  ]+DT/DZ*delp;
           WN[ix  ][iy  ][iz-1]=WN[ix  ][iy  ][iz-1]-DT/DZ*delp;
         }
       }
     }
     /* 圧力の相対性に関する処理(IRELP=1なら以下の処理を行う) */
     if (IRELP==1){
       postn = PO[1][1][1];
       for (ix=1;ix<=NX;ix++){
         for (iy=1;iy<=NY;iy++){
           for (iz=1;iz<=NZ;iz++){
              PO[ix][iy][iz] = PO[ix][iy][iz]-postn;
           }
         }
       }
     }
     /* IFLG=1なら，連続の式を満たしていないと判定し再び圧力計算を行う */
     if ( fabs(DMAX) >= EPSP ){
       IFLG = 1;
     }
     /* 圧力計算の回数を100回ごとに表示 */
     if ( (ITR%100) == 0){
       printf (" Iteration= %8d, Div(max)( %6d, %6d, %6d ) = %15.6lE\n",ITR,ixmax,iymax,izmax,DMAX);
     }
     /*新たに得られた速度を用いて境界条件を処理する*/
     velbnd();
}
/*温度場の計算*/
void caltem()
{
     int ix,iy,iz;
     double uut,vvt,wwt,cnvtx,cnvty,cnvtz,dift;
     /*T(ix,iy,iz)の計算*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
         /* uut,vvt,wwtはそれぞT(ix,iy,iz)におけるU,V,Wの補間値 */
         uut = ( UO[ix][iy][iz] + UO[ix-1][iy  ][iz  ] ) / 2.0;
         vvt = ( VO[ix][iy][iz] + VO[ix  ][iy-1][iz  ] ) / 2.0;
         wwt = ( WO[ix][iy][iz] + WO[ix  ][iy  ][iz-1] ) / 2.0;
         /*対流項(cnvtx,cnvty,cnvtz)を１次精度風上差分にて計算*/
         if ( uut >= 0.0 ){
           cnvtx = uut*( TO[ix][iy][iz] - TO[ix-1][iy][iz] ) / DX;
         }
         else if ( uut < 0.0 ){
           cnvtx = uut*( TO[ix+1][iy][iz] - TO[ix][iy][iz] ) / DX;
         }
         if ( vvt >= 0.0 ){
           cnvty = vvt*( TO[ix][iy][iz] - TO[ix][iy-1][iz] ) / DY;
         }
         else if ( vvt < 0.0 ){
           cnvty = vvt*( TO[ix][iy+1][iz] - TO[ix][iy][iz] ) / DY;
         }
         if ( wwt >= 0.0 ){
           cnvtz = wwt*( TO[ix][iy][iz] - TO[ix][iy][iz-1] ) / DZ;
         }
         else if ( wwt < 0.0 ){
           cnvtz = wwt*( TO[ix][iy][iz+1] - TO[ix][iy][iz] ) / DZ;
         }
         /* 拡散項(dift)の計算 */
         dift = ALP*( (TO[ix-1][iy][iz]-2.0*TO[ix][iy][iz]+TO[ix+1][iy][iz])/(DX*DX)
                     +(TO[ix][iy-1][iz]-2.0*TO[ix][iy][iz]+TO[ix][iy+1][iz])/(DY*DY)
                     +(TO[ix][iy][iz-1]-2.0*TO[ix][iy][iz]+TO[ix][iy][iz+1])/(DZ*DZ) );
         /* 次の時間のTの計算 */
         TN[ix][iy][iz] = TO[ix][iy][iz] + DT*( -cnvtx-cnvty-cnvtz+dift );
         }
       }
     }
     /* 境界条件の処理 */
     tbnd();
}
/* 速度の境界条件の処理 */
void velbnd()
{
     int ix,iy,iz;
     for (iy=1;iy<=NY;iy++){/*U（右側面）*/
       for (iz=1;iz<=NZ;iz++){
         UN[NX][iy][iz] = 0.0;
       }
     }
     for (iy=1;iy<=NY;iy++){/*U（左側面）*/
       for (iz=1;iz<=NZ;iz++){
         UN[0][iy][iz] = 0.0;
       }
     }
     for (ix=0;ix<=NX;ix++){/*U（後面）*/
       for (iy=1;iy<=NY;iy++){
         UN[ix][iy][NZ+1] = -UN[ix][iy][NZ];
       }
     }
     for (ix=0;ix<=NX;ix++){/*U（前面）*/
       for (iy=1;iy<=NY;iy++){
         UN[ix][iy][0] = -UN[ix][iy][1];
       }
     }
     for (ix=0;ix<=NX;ix++){/*U（上面）*/
       for (iz=0;iz<=NZ+1;iz++){
         UN[ix][NY+1][iz] = -UN[ix][NY][iz];
       }
     }
     for (ix=0;ix<=NX;ix++){/*U（下面）*/
       for (iz=0;iz<=NZ+1;iz++){
         UN[ix][0][iz] = -UN[ix][1][iz];
       }
     }
     for (iy=1;iy<=NY-1;iy++){/*V（右側面）*/
       for (iz=1;iz<=NZ;iz++){
         VN[NX+1][iy][iz] = -VN[NX][iy][iz];
       }
     }
     for (iy=1;iy<=NY-1;iy++){/*V（左側面）*/
       for (iz=1;iz<=NZ;iz++){
         VN[0][iy][iz] = -VN[1][iy][iz];
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*V（上面）*/
       for (iz=0;iz<=NZ+1;iz++){
         VN[ix][NY][iz] = 0.0;
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*V（下面）*/
       for (iz=0;iz<=NZ+1;iz++){
         VN[ix][0][iz] = 0.0;
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*V（後面）*/
       for (iy=1;iy<=NY-1;iy++){
         VN[ix][iy][NZ+1] = -VN[ix][iy][NZ];
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*V（前面）*/
       for (iy=1;iy<=NY-1;iy++){
         VN[ix][iy][0] = -VN[ix][iy][1];
       }
     }
     for (iy=1;iy<=NY;iy++){/*W（右面）*/
       for (iz=1;iz<=NZ-1;iz++){
          WN[NX+1][iy][iz] = -WN[NX][iy][iz];
       }
     }
     for (iy=1;iy<=NY;iy++){/*W（左面）*/
       for (iz=1;iz<=NZ-1;iz++){
          WN[0][iy][iz] = -WN[1][iy][iz];
       }
     }
     for (ix=0;ix<=NX;ix++){/*W（後面）*/
       for (iy=0;iy<=NY+1;iy++){
          WN[ix][iy][NZ] = 0.0;
       }
     }
     for (ix=0;ix<=NX;ix++){/*W（前面）*/
       for (iy=0;iy<=NY+1;iy++){
         WN[ix][iy][0] = 0.0;
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*W（上面）*/
       for (iz=1;iz<=NZ-1;iz++){
         WN[ix][NY+1][iz] = -WN[ix][NY][iz];
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*W（下面）*/
       for (iz=1;iz<=NZ-1;iz++){
         WN[ix][0][iz] = -WN[ix][1][iz];
       }
     }
}
/*温度の境界条件の処理*/
void tbnd()
{
     int ix,iy,iz;
      for (iy=0;iy<=NY+1;iy++){/*右側面*/
        for (iz=0;iz<=NZ+1;iz++){
          TN[NX+1][iy][iz] = 2.0 * ( -0.5 ) - TN[NX][iy][iz];
        }
      }
      for (iy=0;iy<=NY+1;iy++){/*左側面*/
        for (iz=0;iz<=NZ+1;iz++){
          TN[0][iy][iz] = 2.0 * ( +0.5 ) - TN[1][iy][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/*上面*/
        for (iz=1;iz<=NZ;iz++){
          TN[ix][NY+1][iz] = TN[ix][NY][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/*下面*/
        for (iz=1;iz<=NZ;iz++){
          TN[ix][0][iz] = TN[ix][1][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/*後面*/
        for (iy=0;iy<=NY+1;iy++){
          TN[ix][iy][NZ+1] = TN[ix][iy][NZ];
        }
      }
      for (ix=1;ix<=NX;ix++){/*前面*/
        for (iy=0;iy<=NY+1;iy++){
          TN[ix][iy][0] = TN[ix][iy][1];
        }
      }
}
void prout()/*データ出力用*/
{
    fwrite(UN, sizeof(double), NX1*NY2*NZ2, out_11);
    fwrite(VN, sizeof(double), NX2*NY1*NZ2, out_12);
    fwrite(WN, sizeof(double), NX2*NY2*NZ1, out_13);
    fwrite(PO, sizeof(double), NX2*NY2*NZ2, out_14);
    fwrite(TN, sizeof(double), NX2*NY2*NZ2, out_15);
}
/* Tecplot用データの出力 */
void tecplt ( char *fname_tec )
{

     int ix,iy,iz;
     double x,y,z,u,v,w,t;

     out_21=fopen(fname_tec,"wt");
     fprintf (out_21," VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"T\"\n");
     fprintf (out_21," ZONE I=%3d, J=%3d, K=%3d, F=POINT \n",NX1,NY1,NZ1);

     for (iz=0;iz<=NZ;iz++){
       for (iy=0;iy<=NY;iy++){
         for (ix=0;ix<=NX;ix++){
           x = DX * (double)ix;
           y = DY * (double)iy;
           z = DZ * (double)iz;
           u = (UN[ix][iy][iz]+UN[ix][iy+1][iz]+UN[ix][iy][iz+1]+UN[ix][iy+1][iz+1])/4.0;
           v = (VN[ix][iy][iz]+VN[ix+1][iy][iz]+VN[ix][iy][iz+1]+VN[ix+1][iy][iz+1])/4.0;
           w = (WN[ix][iy][iz]+WN[ix+1][iy][iz]+WN[ix][iy+1][iz]+WN[ix+1][iy+1][iz])/4.0;
           t = (TN[ix][iy][iz]+TN[ix+1][iy][iz]+TN[ix][iy+1][iz]+TN[ix+1][iy+1][iz]
               +TN[ix][iy][iz+1]+TN[ix+1][iy][iz+1]+TN[ix][iy+1][iz+1]+TN[ix+1][iy+1][iz+1] )/8.0;
         fprintf(out_21,"%12.3lE %12.3lE %12.3lE %12.3lE %12.3lE %12.3lE %12.3lE\n",x,y,z,u,v,w,t);
         }
       }
     }
     fclose (out_21);
}
