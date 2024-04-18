/********************************************************************
* ファイル名  ：smac3d.c                                            *
* タイトル    ：SMAC法による3次元熱流動解析プログラム               *
* 製作者      ：平野　博之                                          *
* 所属        ：岡山理科大学 工学部 バイオ・応用化学科              *
* 製作日      ：2011.12.25                                          *
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
#define NE 8000 /* NE:圧力補正の未知数の数 --- SMAC --- */
#define NX1 21 /* NX1=NX+1 */
#define NX2 22 /* NX2=NX+2 */
#define NY1 21 /* NY1=NY+1 */
#define NY2 22 /* NY2=NY+2 */
#define NZ1 21 /* NZ1=NZ+1 */
#define NZ2 22 /* NZ2=NZ+2 */
#define NE1 8001 /* NE1=NE+1 */
#define NXY 400 /* NXY=NX*NY --- SMAC --- */
#define NXYNXY 801 /* NXYNXY=2*NXY+1 ガウスの消去法で使用 --- SMAC --- */
void cinit();
void adv();
void calvel();
void press();
void caltem();
void velbnd();
void tbnd();
void prout();
void tecplt( char *fname_tec);
void pdbndc (double [NE1], double [NE1], double [NE1], double [NE1], double [NE1],
             double [NE1], double [NE1], double [NE1]);
void gb (double [NE1], double [NE1], double [NE1], double[NE1], double [NE1],
         double [NE1], double [NE1], double [NE1], double [NE1],
           int nx, int ny, int nz, int ne, int nxy);
void crb   (double [NE1], double [NE1], double [NE1], double[NE1], double [NE1],
            double [NE1], double [NE1], double [NE1], double [NE1],
            int nx, int ny, int nz, int ne, int nxy);
void promv (double [NE1], double [NE1], double [NE1], double [NE1],
            double [NE1], double [NE1], double [NE1],
            double [NE1], double [NE1],
            int nx, int ny, int nz, int ne, int nxy);
void bicgb (double [NE1], double [NE1], double [NE1], double[NE1], double[NE1],
            double [NE1], double [NE1],
            double [NE1], double [NE1], 
            int nx, int ny, int nz, int ne, int nxy);
double provv (double [NE1], double [NE1], int ne);
void psorb (double [NE1], double [NE1], double [NE1], double[NE1], double[NE1],
            double [NE1], double [NE1],
            double [NE1], double [NE1],
            int nx, int ny, int nz, int ne, int nxy);
void lsorb (double [NE1] ,double [NE1], double [NE1], double [NE1], double [NE1], 
            double [NE1], double [NE1],
            double [NE1], double [NE1],
            int nx, int ny, int nz, int ne, int nxy);
void pdbnd();
void pbnd();
double DX, DY, DZ, DT;
double VIS, ALP, BUO;
double RE, PR, GR, TIME, OMG, EPSP;
int  ICYCLE, ITR, IFLG, IRELP, METHOD;
double DMAX;
int ITYPE;
int NITR;
char fname[12][80];
char *fname10;
FILE *in_10, *out_11, *out_12, *out_13, *out_14, *out_15, *out_21;
FILE *in_16, *in_17, *in_18, *in_19, *in_20;
double UO[NX1][NY2][NZ2],UN[NX1][NY2][NZ2],VO[NX2][NY1][NZ2],VN[NX2][NY1][NZ2];
double WO[NX2][NY2][NZ1],WN[NX2][NY2][NZ1],PO[NX2][NY2][NZ2],PN[NX2][NY2][NZ2];
double TO[NX2][NY2][NZ2],TN[NX2][NY2][NZ2];
double DIV[NX1][NY1][NZ1], PD[NX2][NY2][NZ2];/* --- SMAC ---*/
void main()
{
     int i,j,NCYCLE;
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
     printf(" DLX=%12.3lE DLY=%12.3lE DLZ=%12.3lE IRELP=%d METHOD=%d \n",DLX,DLY,DLZ,IRELP,METHOD);
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
     IFLG=0;
     /*  --- SMAC ---
      本計算においては，ICYCLE=0のときの初期条件は速度場と圧力場をゼロと
      しているので，圧力補正の線形システムを解くのが困難となるため，
      ICYCLE=1のときだけは温度の計算へジャンプするようにする．
      計算条件や問題に応じて適宜削除あるいは変更 */
      if ( ICYCLE==1 ){
        goto label_720;
      }
     calvel();/*速度場の計算*/
     press();/*仮の圧力場の計算*/
     label_720:{};/* 速度場と圧力場がゼロとなるときのスキップ先 --- SMAC --- */
     if ( IFLG==0 ){/*圧力場の計算が収束したとき*/
       if ( ITYPE==2 ){/*非等温場計算の場合*/
         caltem();/*温度場を計算*/
       }
     }
     if ( IFLG==1 ){/*圧力場の計算が収束していないとき*/
       /*データを出力して強制終了*/
       printf(" calculation has diverged \n");
       prout();
       goto label_900;
     }
     if ( ICYCLE < NCYCLE ){/*時間進行カウンタ(ICYCLE)がNCYCLEより小さい時*/
       goto label_700;
     }
     else{/*時間進行カウンタがNCYCLEになったら->計算終了*/
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
             PD[ix][iy][iz] = 0.0; PN[ix][iy][iz] = 0.0;/* --- SMAC ---*/
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
       fread(PN, sizeof(double), NX2*NY2*NZ2, in_19);
       fread(PD, sizeof(double), NX2*NY2*NZ2, in_19);
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
         for (iz=0;iz<=NZ+1;iz++){/*UO:新しい時間ステップでの初期値，UNを保存．*/
           UO[ix][iy][iz] = UN[ix][iy][iz];
         }
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*    VN -> VO : 必要なら入れ替える前にVNとVOから変動量を求める*/
       for (iy=0;iy<=NY;iy++){/*    VN:前の時間ステップにおいて最終的に得られた値，圧力補正の度に更新される*/
         for (iz=0;iz<=NZ+1;iz++){/*VO:新しい時間ステップでの初期値，VNを保存．*/
           VO[ix][iy][iz] = VN[ix][iy][iz];
         }
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*WN -> WO : 必要なら入れ替える前にWNとWOから変動量を求める*/
       for (iy=0;iy<=NY+1;iy++){/*WN:前の時間ステップにおいて最終的に得られた値，圧力補正の度に更新される*/
         for (iz=0;iz<=NZ;iz++){/*WO:新しい時間ステップでの初期値，WNを保存．*/
           WO[ix][iy][iz] = WN[ix][iy][iz];
         }
       }
     }
     for (ix=0;ix<=NX+1;ix++){/* TN -> TO : 必要なら入れ替える前にTNとTOから変動量を求める*/
       for (iy=0;iy<=NY+1;iy++){/*  TN:前の時間ステップでの計算値*/
         for (iz=0;iz<=NZ+1;iz++){/*TO:新しい時間ステップでの初期値．TNを保存．*/
           TO[ix][iy][iz] = TN[ix][iy][iz];
           PO[ix][iy][iz] = PN[ix][iy][iz];/* --- SMAC --- PO:新しい時間ステップでの初期値．PNを保存．*/
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
/*圧力場の計算 Fortran プログラムのコメントを参照してください．*/
void press()
{
     /* --- SMAC --- 圧力補正の線形システム用 */
     double a1[NE1],a2[NE1],a3[NE1],a4[NE1],a5[NE1],a6[NE1],a7[NE1],b[NE1],x[NE1];
     int ix,iy,iz,k,nnx,nny,nnz,nne,nnxy;
     /* 線形システム関係の配列の初期化 --- SMAC --- */
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
         k = iy + (ix-1)*NY + (iz-1)*NXY;
         a1[k]=0.0; a2[k]=0.0; a3[k]=0.0; a4[k]=0.0; a5[k]=0.0; a6[k]=0.0; a7[k]=0.0;
         b[k]=0.0; x[k]=0.0;
         }
       }
     }
     /*P(ix,iy,iz)の計算*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
           /* --- SMAC ---*/
           DIV[ix][iy][iz] = ( UN[ix][iy][iz] - UN[ix-1][iy][iz] )/DX
                           + ( VN[ix][iy][iz] - VN[ix][iy-1][iz] )/DY
                           + ( WN[ix][iy][iz] - WN[ix][iy][iz-1] )/DZ;
         }
       }
     }
     /* 係数行列の作成 --- SMAC ---
        計算領域の1点さらに内側の点（仮想セルを含めて考えると2点を除く内側の点）
        に関して係数行列を作成．残りは境界条件を反映させて pdbndc で設定する */
     for (ix=2;ix<=NX-1;ix++){
       for (iy=2;iy<=NY-1;iy++){
         for (iz=2;iz<=NZ-1;iz++){
         k = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[k] = DT*( 2.0/(DX*DX) + 2.0/(DY*DY) + 2.0/(DZ*DZ));
         a4[k] = - 1.0 / (DY*DY) * DT;
         a2[k] = - 1.0 / (DY*DY) * DT;
         a1[k] = - 1.0 / (DX*DX) * DT;
         a5[k] = - 1.0 / (DX*DX) * DT;
         a6[k] = - 1.0 / (DZ*DZ) * DT;
         a7[k] = - 1.0 / (DZ*DZ) * DT;
         b[k] = -DIV[ix][iy][iz];
         x[k] = PD[ix][iy][iz];
         }
       }
     }
     /* --- SMAC --- 境界条件を係数行列に反映させる */
     pdbndc(a1,a2,a3,a4,a5,a6,a7,b);
     /* --- SMAC --- 圧力補正 P'(PD) に関するポアソン方程式の解法 */
     nnx = NX; nny = NY; nnz= NZ; nne = NE; nnxy = NXY;
     /* 1. 直接法 : バンドマトリックスによるガウスの消去法 */
     if (METHOD==1) gb    (a1,a2,a3,a4,a5,a6,a7,b,x,nnx,nny,nnz,nne,nnxy);
     /* 2. 反復法1 : point-SOR 法 */
     if (METHOD==2) psorb (a1,a2,a3,a4,a5,a6,a7,b,x,nnx,nny,nnz,nne,nnxy);
     /* 3. 反復法2 : line-SOR 法 : OMG=1で十分．大きくしすぎると発散する */
     if (METHOD==3) lsorb (a1,a2,a3,a4,a5,a6,a7,b,x,nnx,nny,nnz,nne,nnxy);
     /* 4. クリロフ部分空間法1 : 共役残差法 */
     if (METHOD==4) crb   (a1,a2,a3,a4,a5,a6,a7,b,x,nnx,nny,nnz,nne,nnxy);
     /* 5. クリロフ部分空間法2 : BiCGSTAB */
     if (METHOD==5) bicgb (a1,a2,a3,a4,a5,a6,a7,b,x,nnx,nny,nnz,nne,nnxy);
     /* METHOD=4,5における探索ベクトル計算時のゼロ除算対策 */
     if (IFLG==2)   psorb (a1,a2,a3,a4,a5,a6,a7,b,x,nnx,nny,nnz,nne,nnxy);
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
         k = iy + (ix-1)*NY + (iz-1)*NXY;
         /* 圧力の相対性の処理 : 基準値の設定
            1次独立な解を求める場合は圧力の基準点が強制的に設定される */
         if (IRELP != 2) PD[ix][iy][iz] = x[k]; /* 圧力の基準点を設けない場合 */
         if (IRELP == 2) PD[ix][iy][iz] = x[k]-x[1]; /* P'(1,1,1)=0 ---> P(1,1,1)=0 */
         }
       }
     }
     /* 圧力補正の境界条件の処理 */
     pdbnd();
     /* 速度の修正 */
     for (ix=1;ix<=NX-1;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
          UN[ix][iy][iz]=UN[ix][iy][iz]+(PD[ix][iy][iz]-PD[ix+1][iy][iz])/DX*DT;
         }
       }
     }
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY-1;iy++){
         for (iz=1;iz<=NZ;iz++){
          VN[ix][iy][iz]=VN[ix][iy][iz]+(PD[ix][iy][iz]-PD[ix][iy+1][iz])/DY*DT;
         }
       }
     }
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ-1;iz++){
          WN[ix][iy][iz]=WN[ix][iy][iz]+(PD[ix][iy][iz]-PD[ix][iy][iz+1])/DZ*DT;
         }
       }
     }
     /*新たに得られた速度を用いて境界条件を処理する*/
     velbnd();
     /* 圧力の修正 */
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
          PN[ix][iy][iz] = PO[ix][iy][iz] + PD[ix][iy][iz];
         }
       }
     }
     /* 新たに得られた圧力を用いて境界条件を処理する */
     pbnd();
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
/* 圧力補正の係数行列と境界条件
   本計算プログラムでは，圧力補正の境界条件を，係数行列作成時に考慮する．
   また，解の1次独立性の確保もここで処理する．*/
void pdbndc(double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1], 
            double a5[NE1], double a6[NE1], double a7[NE1], double b[NE1] )
{
     int ix,iy,iz,m;
     /* --- SMAC --- 係数行列に境界条件を反映させる
     計算領域の左側の係数行列（境界条件を反映）*/
     ix = 1;
     for (iy=2;iy<=NY-1;iy++){
       iz=1;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       for (iz=2;iz<=NZ-1;iz++){
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) + 2.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       }
       iz=NZ;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
     }
     /* 計算領域の右側の係数行列（境界条件を反映）*/
     ix = NX;
     for (iy=2;iy<=NY-1;iy++){
       iz=1;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) +1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       for (iz=2;iz<=NZ-1;iz++){
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) + 2.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       }
       iz=NZ;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
     }
     /* 計算領域の下側の係数行列（境界条件を反映）*/
     iy = 1;
     for (ix=2;ix<=NX-1;ix++){
       iz=1;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       for (iz=2;iz<=NZ-1;iz++){
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 2.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       }
       iz=NZ;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
     }
     /* 計算領域の上側の係数行列（境界条件を反映）*/
     iy = NY;
     for (ix=2;ix<=NX-1;ix++){
       iz=1;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) +1.0/(DZ*DZ) );
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       for (iz=2;iz<=NZ-1;iz++){
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 2.0/(DZ*DZ) );
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       }
       iz=NZ;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
     }
     /* 計算領域の後面の係数行列（境界条件を反映）*/
     iz = 1;
     for (ix=2;ix<=NX-1;ix++){
       iy=1;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       for (iy=2;iy<=NY-1;iy++){
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 2.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       }
       iy=NY;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
     }
     /* 計算領域の前面の係数行列（境界条件を反映）*/
     iz = NZ;
     for (ix=2;ix<=NX-1;ix++){
       iy=1;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       for (iy=2;iy<=NY-1;iy++){
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 2.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       }
       iy=NY;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
     }
     /* 計算領域の左上点列の係数行列（境界条件を反映）*/
     ix=1; iy=NY;
     iz=1;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a2[m] = - 1.0 / (DY*DY) * DT;
       a5[m] = - 1.0 / (DX*DX) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     for (iz=2;iz<=NZ-1;iz++){
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 2.0/(DZ*DZ) );
       a2[m] = - 1.0 / (DY*DY) * DT;
       a5[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     }
     iz=NZ;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a2[m] = - 1.0 / (DY*DY) * DT;
       a5[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     /* 計算領域の左下点列の係数行列（境界条件を反映）*/
     ix=1; iy=1;
     iz=1;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a4[m] = - 1.0 / (DY*DY) * DT;
       a5[m] = - 1.0 / (DX*DX) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     for (iz=2;iz<=NZ-1;iz++){
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 2.0/(DZ*DZ) );
       a4[m] = - 1.0 / (DY*DY) * DT;
       a5[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     }
     iz=NZ;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a4[m] = - 1.0 / (DY*DY) * DT;
       a5[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     /* 計算領域の右上点列の係数行列（境界条件を反映）*/
     ix=NX; iy=NY;
     iz=1;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a2[m] = - 1.0 / (DY*DY) * DT;
       a1[m] = - 1.0 / (DX*DX) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     for (iz=2;iz<=NZ-1;iz++){
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 2.0/(DZ*DZ) );
       a2[m] = - 1.0 / (DY*DY) * DT;
       a1[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     }
     iz=NZ;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a2[m] = - 1.0 / (DY*DY) * DT;
       a1[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     /* 計算領域の右下点列の係数行列（境界条件を反映）*/
     ix=NX; iy=1;
     iz=1;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a4[m] = - 1.0 / (DY*DY) * DT;
       a1[m] = - 1.0 / (DX*DX) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     for (iz=2;iz<=NZ-1;iz++){
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 2.0/(DZ*DZ) );
       a4[m] = - 1.0 / (DY*DY) * DT;
       a1[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     }
     iz=NZ;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a4[m] = - 1.0 / (DY*DY) * DT;
       a1[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     /* １次独立な解を得るための処理 : 直接法においては必須
        (ix=1,iy=1,iz ---> m=1を基準点とし，常にPN[1][1][1]=PD[1][1][1]=0とする) */
     if (IRELP ==1) {
       a3[1] = 1.0; a4[1] = 0.0; a5[1] = 0.0; a7[1]=0.0; b[1] = 0.0;
       /* m=2 の点の処理(m=1とのリンクを断つ) */
       a3[2] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) + 1.0/(DZ*DZ) );
       a4[2] = - 1.0 / (DY*DY) * DT;
       a2[2] = 0.0;
       a5[2] = - 1.0 / (DX*DX) * DT;
       a7[2] = - 1.0 / (DZ*DZ) * DT;
       /* m=1+NY の点の処理(m=1とのリンクを断つ) */
       a3[1+NY] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a4[1+NY] = - 1.0 / (DY*DY) * DT;
       a1[1+NY] = 0.0;
       a5[1+NY] = - 1.0 / (DX*DX) * DT;
       a7[1+NY] = - 1.0 / (DZ*DZ) * DT;
       /* m=1+NX*NY の点の処理(m=1とのリンクを断つ) */
       a3[1+NXY] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a4[1+NXY] = - 1.0 / (DY*DY) * DT;
       a5[1+NXY] = - 1.0 / (DX*DX) * DT;
       a6[1+NXY] = 0.0;
       a7[1+NXY] = - 1.0 / (DZ*DZ) * DT;
     }
}
/********************************************************************
*                                                                   *
*                  圧力補正の境界条件の処理                         *
*                                                                   *
* 本計算プログラムでは，係数行列作成時に境界条件を考慮しているので，*
* 実質的には，最終的に仮想セルの値を配列に格納しているにすぎない．  *
* このルーチンはなくともよい．                                      *
* 反復過程で境界条件を考慮する場合は必須で，反復のたびにこの処理を  *
* 行う必要が生じる．                                                *
*                                                                   *
********************************************************************/
void pdbnd()
{
     int ix,iy,iz;
      for (iy=0;iy<=NY+1;iy++){/*右側面*/
        for (iz=0;iz<=NZ+1;iz++){
          PD[NX+1][iy][iz] = PD[NX][iy][iz];
        }
      }
      for (iy=0;iy<=NY+1;iy++){/*左側面*/
        for (iz=1;iz<=NZ+1;iz++){
          PD[0][iy][iz] = PD[1][iy][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/*上面*/
        for (iz=0;iz<=NZ;iz++){
          PD[ix][NY+1][iz] = PD[ix][NY][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/*下面*/
        for (iz=1;iz<=NZ;iz++){
          PD[ix][0][iz] = PD[ix][1][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/* 後面 */
        for (iy=1;iy<=NY;iy++){
          PD[ix][iy][NZ+1] = PD[ix][iy][NZ];
        }
      }
      for (ix=1;ix<=NX;ix++){/* 前面 */
        for (iy=1;iy<=NY;iy++){
          PD[ix][iy][0] = PD[ix][iy][1];
        }
      }
}
/*********************************************************************
*                     圧力の境界条件の処理
*********************************************************************/
void pbnd()
{
     int ix,iy,iz;
      for (iy=0;iy<=NY+1;iy++){/*右側面*/
        for (iz=0;iz<=NZ+1;iz++){
          PN[NX+1][iy][iz] = PN[NX][iy][iz];
        }
      }
      for (iy=0;iy<=NY+1;iy++){/*左側面*/
        for (iz=1;iz<=NZ+1;iz++){
          PN[0][iy][iz] = PN[1][iy][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/*上面*/
        for (iz=0;iz<=NZ;iz++){
          PN[ix][NY+1][iz] = PN[ix][NY][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/*下面*/
        for (iz=1;iz<=NZ;iz++){
          PN[ix][0][iz] = PN[ix][1][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/* 後面 */
        for (iy=1;iy<=NY;iy++){
          PN[ix][iy][NZ+1] = PN[ix][iy][NZ];
        }
      }
      for (ix=1;ix<=NX;ix++){/* 前面 */
        for (iy=1;iy<=NY;iy++){
          PN[ix][iy][0] = PN[ix][iy][1];
        }
      }
}
void prout()/*データ出力用*/
{
    fwrite(UN, sizeof(double), NX1*NY2*NZ2, out_11);
    fwrite(VN, sizeof(double), NX2*NY1*NZ2, out_12);
    fwrite(WN, sizeof(double), NX2*NY2*NZ1, out_13);
    fwrite(PN, sizeof(double), NX2*NY2*NZ2, out_14);
    fwrite(PD, sizeof(double), NX2*NY2*NZ2, out_14);
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
/********************************************************************
*                    各種の線形システム解法                         *
*                                                                   *
*  1. ガウスの消去法                                                *
*  2. point-SOR 法                                                  *
*  3. line-SOR 法                                                   *
*  4. 共役残差法                                                    *
*  5. Bi-CGSTAB法                                                   *
*                                                                   *
*  いずれも，3次元のポアソン方程式を7点差分近似にて離散化した       *
*  線形システムを解くためのもので，最適化してある                   *
*  いずれのサブルーチンも同一引数としてある                         *
*                                                                   *
********************************************************************/
/********************************************************************
*                                                                   *
*   3次元ラプラシアン離散化による7点差分近似                        *
*   にて得られた規則的非対称行列Aを含んだ線形システム               *
*                        AX=B                                       *
*   をGaussの消去法を用いて解くサブルーチン．(軸選択無)             *
*   Aはバンドマトリックス                                           *
*                                                                   *
*    線形システム --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    係数行列の計算容量節約：詳細はサブルーチン PRESS を参照        *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE)  *
*            ---> A(-NXY:NXY,NE)に格納しなおす                      *
*                                                                   *
*    B : 既知ベクトル                                               *
*    X : 未知ベクトル ---> これを求める                             *
*                                                                   *
*［変数の説明］                                                     *
*    NX : x方向格子分割数                                           *
*    NY : y方向格子分割数                                           *
*    NZ : z方向格子分割数                                           *
*    NXY: NX * NY                                                   *
*    NE : 総格子点数 = NX * NY * NZ                                 *
*                                                                   *
********************************************************************/
void gb (double a1[NE1], double a2[NE1], double a3[NE1],
         double a4[NE1], double a5[NE1], double a6[NE1], double a7[NE1],
         double b[NE1], double x[NE1],
         int nx, int ny, int nz, int ne, int nxy)
{
     double a[NXYNXY][NE1];
     int ine,i,j,n,k;
     double aa,s;
     /* 直接法のときは(係数行列が特異でなければ)反復なしで，必ず解を得る */
     IFLG = 0;
     /* マトリックスAのゼロクリア */
     for (ine=1;ine<=ne;ine++){
       for (i=-nxy;i<=nxy;i++){
         a[i+nxy][ine] = 0.0;
       }
     }
     /* 必要なところにA1からA7までを格納する */
     for (ine=1;ine<=ne;ine++){
       a[ -ny+nxy][ine] = a1[ine];
       a[  -1+nxy][ine] = a2[ine];
       a[   0+nxy][ine] = a3[ine];
       a[   1+nxy][ine] = a4[ine];
       a[  ny+nxy][ine] = a5[ine];
       a[-nxy+nxy][ine] = a6[ine];
       a[ nxy+nxy][ine] = a7[ine];
     }
/* 前進消去 */
     for (i=1;i<=ne-1;i++){
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
/* 後退代入 */
     /* 係数行列の特異性を判定 */
     if ( fabs( a[0+nxy][ne] ) <= 1.0e-50 ){
        printf("  Matrix is singular : |A(0,NE)| < 1E-50 \n");
        IFLG = 1; /* 圧力補正の計算を発散に設定 */
     }
     x[ne] = b[ne] / a[0+nxy][ne];
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
}
/********************************************************************
*                                                                   *
*  point-SOR 法による非対称行列 A を含む線形システム解法サブルーチン*
*  3次元ラプラシアン離散化による7点差分近似用                       *
*                                                                   *
*    線形システム --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    係数行列の計算容量節約：詳細はサブルーチン PRESS を参照        *
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],                  *
*                 a5[NE1],a6[NE1],a7[NE1]                           *
*                                                                   *
*    B : 既知ベクトル                                               *
*    X : 未知ベクトル ---> これを求める                             *
*                                                                   *
*［変数の説明］                                                     *
*    NX : x方向格子分割数                                           *
*    NY : y方向格子分割数                                           *
*    NZ : z方向格子分割数                                           *
*    NXY: NX * NY                                                   *
*    NE : 総格子点数 = NX * NY * NZ                                 *
*    NITR : 許容反復回数(in3d.mac)にて設定                          *
*    EPSP : 収束判定条件で用いる値(in3d.mac)にて設定                *
*                                                                   *
* [収束判定条件]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : 前の反復による値                                 *
*  \vec{x}^{new} : 新しい反復による値                               *
*                                                                   *
********************************************************************/
void psorb (double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1], 
            double a5[NE1], double a6[NE1], double a7[NE1],
            double b[NE1], double x[NE1],
            int nx, int ny, int nz, int ne, int nxy)
{
     int i,j;
     double rnorm1, xold, sum, xnew;

     IFLG=1;/* FLGの初期値は"収束せず" */

     for (j=1;j<=NITR;j++){
       rnorm1 = 0.0;

       i=1;
         xold = x[i];
         sum =                           +a4[i]*x[i+1]+a5[i]*x[i+ny]               +a7[i]*x[i+nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );

       for (i=2;i<=ny;i++){
         xold = x[i];
         sum =               a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny]               +a7[i]*x[i+nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }

       for (i=ny+1;i<=nxy;i++){
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny]               +a7[i]*x[i+nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }

       for (i=nxy+1;i<=ne-nxy;i++){
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny]+a6[i]*x[i-nxy]+a7[i]*x[i+nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }

       for (i=ne-nxy+1;i<=ne-ny;i++){
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny]+a6[i]*x[i-nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }

       for (i=ne-ny+1;i<=ne-1;i++){
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1]              +a6[i]*x[i-nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }

       i=ne;
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]                           +a6[i]*x[i-nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
 
       if ( rnorm1 <= EPSP ){/* 収束判定; 収束なら IFLG=0 に設定 */
         IFLG=0; ITR = j; goto label_700;
       }
     }
 
     label_700:{};/* 収束と判定されたときの分岐点 */
}
/********************************************************************
*                                                                   *
*  line-SOR 法による非対称行列 A を含む線形システム解法サブルーチン *
*  3次元ラプラシアン離散化による7点差分近似用                       *
*                                                                   *
*    線形システム --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    係数行列の計算容量節約：詳細はサブルーチン PRESS を参照        *
*    A_{i,j} ---> at1[NE1],at2[NE1],at3[NE1],at4[NE1],at5[NE1],     *
*                 at6[NE1],at7[NE1]                                 *
*                                                                   *
*    B -> bx[NE1]: 既知ベクトル                                     *
*    X -> xn[NE1]: 未知ベクトル ---> これを求める                   *
*                                                                   *
*［変数の説明］                                                     *
*    NX : x方向格子分割数                                           *
*    NY : y方向格子分割数                                           *
*    NZ : z方向格子分割数                                           *
*    NXY: NX * NY                                                   *
*    NE : 総格子点数 = NX * NY * NZ                                 *
*    NITR : 許容反復回数(in3d.mac)にて設定                          *
*    EPSP : 収束判定条件で用いる値(in3d.mac)にて設定                *
*           1.0で十分．                                             *
*           注意：Point-SORと異なり，あまり大きくしすぎると発散する *
*                                                                   *
* [収束判定条件]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : 前の反復による値                                 *
*  \vec{x}^{new} : 新しい反復による値                               *
*                                                                   *
********************************************************************/
void lsorb (double at1[NE1] ,double at2[NE1], double at3[NE1],
            double at4[NE1], double at5[NE1], double at6[NE1], double at7[NE1],
            double bx[NE1], double xn[NE1],
            int nx, int ny, int nz, int ne, int nxy)
{
     double x1[NE1],xo[NE1],xold[NE1],a[NE1],b[NE1],c[NE1],d[NE1],u[NE1],y[NE1];
     double rnorm;
     int k, inx, iy, ix, iz, iny, j, i, inz, inzz;

     IFLG=1; /* IFLG の初期値は"収束せず" */

     for (k=1;k<=NITR;k++){
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
           xn[inz]=(1.0-OMG)*xo[inz]+OMG*x1[inx];
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
           xn[inz] = (1.0-OMG)*xo[inz]+OMG*x1[iny];
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
           xn[inz] = (1.0-OMG)*xo[inz]+OMG*x1[inzz];
           inzz = inzz + 1;
           }
         }
       }
       rnorm= 0.0;
       for (i=1;i<=ne;i++){
         rnorm= rnorm + (xn[i]-xold[i])*(xn[i]-xold[i]);
       }
       if (rnorm <= EPSP){/* 収束判定 */
         IFLG=0; ITR=k; goto label_900;
       }
     }
     label_900:{};/* 収束と判定されたときの分岐点 */
}
/********************************************************************
*                                                                   *
*    共役残差(Conjugate Residual)法による非対称行列 A を含む        *
*  線形システム解法サブルーチン                                     *
*                        AX=B                                       *
*                                                                   *
*    係数行列の計算容量節約：詳細はサブルーチン PRESS を参照        *
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],                  *
*                 a5[NE1],a6[NE1],a7[NE1]                           *
*                                                                   *
* B : 既知ベクトル                                                  *
* X : 未知ベクトル ---> これを求める ---> ここでは便宜上配列 XP[NE1]*
*                                                                   *
*［変数の説明］                                                     *
*    NX : x方向格子分割数                                           *
*    NY : y方向格子分割数                                           *
*    NZ : z方向格子分割数                                           *
*    NXY: NX * NY                                                   *
*    NE : 総格子点数 = NX * NY * NZ                                 *
*    NITR : 許容反復回数(in3d.mac)にて設定                          *
*    EPSP : 収束判定条件で用いる値(in3d.mac)にて設定                *
*                                                                   *
*［配列の説明］                                                     *
*    R[NE1]: r_{k} = B - A x_{k}                                    *
*    P[NE1]:  p_{k+1} = r_{k+1} + β_{k} p_{k}, p_{0} = r_{0}       *
*    AP[NE1]: A * P                                                 *
*    AR[NE1]: A * R                                                 *
*                                                                   *
* [収束判定条件]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : 前の反復による値                                 *
*  \vec{x}^{new} : 新しい反復による値                               *
*                                                                   *
********************************************************************/
void crb (double a1[NE1], double a2[NE1], double a3[NE1],
          double a4[NE1], double a5[NE1], double a6[NE1], double a7[NE1],
          double b[NE1], double xp[NE1],
          int nx, int ny, int nz, int ne, int nxy)
{
     double r[NE1], p[NE1], ap[NE1], ar[NE1], x[NE1], xold[NE1];/* 作業用配列 */
     double rap, apap, alp, rnorm, beta, arap;
     int i, k;

     promv (a1,a2,a3,a4,a5,a6,a7,x,r,nx,ny,nz,ne,nxy);/* R に AX を代入 */
     for (i=1;i<=ne;i++){
       r[i]=b[i]-r[i]; p[i]=r[i];/* r_{0}とp_{0}(初期値)の設定 */
       xold[i]=xp[i];/* 前の時刻のXをXOLDに代入 */
     }
     promv (a1,a2,a3,a4,a5,a6,a7,p,ap,nx,ny,nz,ne,nxy);/* APに A p_{0} を代入 */
/* 反復計算 */
     for (k=1;k<=NITR;k++){
       rap = provv (r,ap,ne);/* ( r_{k}, A p_{k} )の計算 => RAP */
       apap = provv (ap,ap,ne);/* ( A p_{k}, A p_{k} )の計算 => APAP */
       if (fabs(apap) < 1.0e-50){
         printf(" 0 division : ALPHA_{K} in Conjugate Residual \n");
         IFLG = 2; return;/* 線形システムの計算を発散に設定 */
       }
       else {
         alp = rap / apap;
       }
       rnorm = 0.0;
       for (i=1;i<=ne;i++){
         x[i]=x[i]+alp*p[i];/* x_{k+1}=x_{k}+α_{k}p_{k} */
         r[i]=r[i]-alp*ap[i];/* r_{k+1}=r_{k}-α_{k}Ap_{k} */
         /* 前の反復との差のノルムの計算 */
         rnorm=rnorm+(x[i]-xold[i])*(x[i]-xold[i]);
         xold[i]=x[i];/* 得られたXをXOLDに代入 */
       }
       /* rnorm が EPSP 以下なら収束とみなして 700 へ */
       if (rnorm <= EPSP){
         IFLG=0; ITR=k; goto label_700;
       }
       /* 収束せずの場合 */
       promv (a1,a2,a3,a4,a5,a6,a7,r,ar,nx,ny,nz,ne,nxy);/* A r_{k+1} の計算 => AR(NE) */
       arap =  provv( ar,ap,ne);/* ( A r_{k+1}, A p_{k} )の計算 => ARAP */
       if (fabs(apap) < 1.0e-50){
         printf(" 0 division : BETA_{K} in Conjugate Residual \n");
         IFLG = 2; return;/* 線形システムの計算を発散に設定 */
       }
       else {
         beta = -arap / apap;/* β_{k} = - ( A r_{k+1}, A p_{k} ) / ( A p_{k}, A p_{k} ) */
       }
       for (i=1;i<=ne;i++){
         p[i]=r[i]+beta*p[i];/* p_{k+1} = r_{k+1} + β_{k} p_{k} */
         ap[i]=ar[i]+beta*ap[i];/* A p_{k+1} = A r_{k+1} + β_{k}A p_{k} */
       }
     }
     IFLG=1;/* NITR まで計算しても収束せず */
     label_700:{}; /* 収束と判定されたときの分岐点 */
     for (i=1;i<=ne;i++){
       xp[i]=x[i];
     }
}
/********************************************************************
*                                                                   *
*    ベクトル A とベクトル B の積の計算サブルーチン                 *
*                        AB=C                                       *
*                                                                   *
*［変数の説明］                                                     *
*    NE : 総格子点数(ベクトル A,B のサイズ)                         *
*    C  : A と B の積(スカラー)                                     *
*                                                                   *
********************************************************************/
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
/********************************************************************
*                                                                   *
*    マトリックス A とベクトル B の積の計算サブルーチン             *
*                        AB=C                                       *
*                                                                   *
*［変数の説明］                                                     *
*    NE : 総格子点数(正方マトリックス A,B,C のサイズ)               *
*    C  : A と B の積(ベクトル)                                     *
*                                                                   *
********************************************************************/
void promv (double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1],
            double a5[NE1], double a6[NE1], double a7[NE1],
            double b[NE1], double c[NE1],
            int nx, int ny, int nz, int ne, int nxy)
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
/********************************************************************
*                                                                   *
*  Bi-CGSTAB 法による非対称行列 A を含む                            *
*  線形システム解法サブルーチン                                     *
*                        AX=B                                       *
*                                                                   *
*    係数行列の計算容量節約：詳細はサブルーチン PRESS を参照        *
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],                  *
*                 a5[NE1],a6[NE1],a7[NE1]                           *
*                                                                   *
*    B : 既知ベクトル                                               *
*    X : 未知ベクトル ---> これを求める                             *
*                                                                   *
*［変数の説明］                                                     *
*    NX : x方向格子分割数                                           *
*    NY : y方向格子分割数                                           *
*    NZ : z方向格子分割数                                           *
*    NXY: NX * NY                                                   *
*    NE : 総格子点数 = NX * NY * NZ                                 *
*    NITR : 許容反復回数(in2d.mac)にて設定                          *
*    EPSP : 収束判定条件で用いる値(in3d.mac)にて設定                *
*                                                                   *
* [配列の説明]                                                      *
*    T[NE] : t_{k} = r_{k} - α_{k} A p_{k}                         *
*    X[NE] : x_{k+1} = x_{k} + α_{k} p_{k} + ξ_{K} t_{k}          *
*    R[NE] : r_{k+1} = t_{k} - ξ_{k} A t_{k}                       *
*            r_{0} = B - A x_{0}                                    *
*    P[NE] : p_{k+1} = r_{k+1} + β_{k} ( p_{k}-ξ_{k} A p_{k} )    *
*            p_{0} = r_{0}                                          *
*    AP[NE] : A * P                                                 *
*    AR[NE] : A * T                                                 *
*                                                                   *
* [収束判定条件]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : 前の反復による値                                 *
*  \vec{x}^{new} : 新しい反復による値                               *
*                                                                   *
********************************************************************/
void bicgb ( double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1], double a5[NE1],
             double a6[NE1], double a7[NE1], double b[NE1], double x[NE1],
             int nx, int ny, int nz, int ne, int nxy)
{
     double r[NE1], ap[NE1], at[NE1], p[NE1], s[NE1], t[NE1], xold[NE1];
     int i,j;
     double sr1, sap, alpha, att, atat, xi, sr2, beta, rnorm;

     for (i=1;i<=ne;i++){
       xold[i] = x[i];
     }

      promv (a1,a2,a3,a4,a5,a6,a7,x,r,nx,ny,nz,ne,nxy);/* R に AX を代入 */
      for (i=1;i<=ne;i++){/* r_{0}とp_{0}(初期値)，そして s=r_{0} の設定 */
        r[i] = b[i] - r[i]; p[i] = r[i]; s[i] = r[i];
      }

      /* 繰り返し計算 */
      for (j=1;j<=NITR;j++){
        sr1=provv (s,r,ne);/* ( s, r_{k} ) の計算 => SR1 */
        promv (a1,a2,a3,a4,a5,a6,a7,p,ap,nx,ny,nz,ne,nxy);/* A p_{k} の計算 => AP(NE) */
        sap=provv (s,ap,ne);/* ( s, A p_{k} ) の計算 => SAP */
        if (fabs(sap) < 1.0e-50){
          printf(" 0 division : ALPHA_{K} in Bi-CGSTAB \n");
          IFLG = 2; return;/* 線形システムの計算を発散に設定 */
        }
        else {
          alpha = sr1/sap;/* α_{k} = ( s, r_{k} ) / ( s, A p_{k} ) */
        }
        for (i=1;i<=ne;i++){
          t[i] = r[i] - alpha*ap[i];/* t_{k} = r_{k} - α_{k} A p_{k} */
        }
        promv (a1,a2,a3,a4,a5,a6,a7,t,at,nx,ny,nz,ne,nxy);/* A t_{k} の計算 => AT(NE) */
        att=provv (at,t,ne);/* ( A t_{k}, t_{k} ) の計算 => ATT */
        atat=provv (at,at,ne);/* ( A t_{k}, A t_{k} ) の計算 => ATAT */
        if (fabs(atat) < 1.0e-50){
          printf(" 0 division : XI_{K} in Bi-CGSTAB \n");
          IFLG = 2; return;/* 線形システムの計算を発散に設定 */
        }
        else {
          xi = att / atat;/* ξ_{k} = ( A t_{k}, t_{k} ) / ( A t_{k}, A t_{k} ) */
        }
        rnorm = 0.0;
        for (i=1;i<=ne;i++){
          x[i] = x[i] + alpha*p[i] + xi*t[i];/* x_{k+1}=x_{k}+α_{k}p_{k}+ξ_{k}t_{k} */
          r[i] = t[i] - xi*at[i];/* r_{k+1}=t_{k}-ξ_{k} A t_{k} */
          /* 前の反復との差のノルムの計算 */
          rnorm = rnorm + (x[i]-xold[i])*(x[i]-xold[i]);
          xold[i] = x[i];/* 得られたXをXOLDに代入 */
        }
        if (rnorm <= EPSP){/* RNORM が EPSP 以下なら収束とみなして 900 へ */
          IFLG=0; ITR=j; goto label_900;
        }
        /* 収束せずの場合 */
        sr2=provv (s,r,ne);/* ( s, r_{k+1} ) の計算 => SR2 */
        /* β_{k} = ( α_{k} / ξ_{k} ) * ( s, r_{k+1} )/( s, r_{k} ) */
        beta = (alpha / xi) * (sr2 / sr1);
        for (i=1;i<=ne;i++){
          /* p_{k+1} = r_{k+1} + β_{k}( p_{k}-ξ_{k} A p_{k} ) */
          p[i] = r[i] + beta * ( p[i] - xi*ap[i] );
        }
      }
      /* NITR まで計算しても収束せず */
      IFLG=1;

      label_900:{};/* 収束と判定されたときの分岐点 */
}
