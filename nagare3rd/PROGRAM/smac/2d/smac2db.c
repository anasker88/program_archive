/********************************************************************
* ファイル名  ：smac2d.c                                            *
* タイトル    ：SMAC法による2次元熱流動解析プログラム               *
* 製作者      ：平野　博之                                          *
* 所属        ：岡山理科大学 工学部 バイオ・応用化学科              *
* 製作日      ：2011.11.01                                          *
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
#define NE 400 /* NE:圧力補正の未知数の数 */
#define NX1 21 /* NX1=NX+1 */
#define NX2 22 /* NX2=NX+2 */
#define NY1 21 /* NY1=NY+1 */
#define NY2 22 /* NY2=NY+2 */
#define NYNY 41 /* NYNY=2*NY+1 */
#define NE1 401 /* NE1=NE+1 */
void cinit();
void adv();
void calvel();
void press();
void caltem();
void velbnd();
void tbnd();
void prout();
void tecplt( char *fname_tec);
void pdbndc (double [NE1], double [NE1], double [NE1], double [NE1], double [NE1], double [NE1]);
void gb (double [NE1], double [NE1], double [NE1], double[NE1], double [NE1],
         double [NE1], double [NE1], int nx, int ny, int ne);
void crb   (double [NE1], double [NE1], double [NE1], double[NE1], double [NE1],
            double [NE1], double [NE1], int nx, int ny, int ne);
void promv (double [NE1], double [NE1], double [NE1], double [NE1],
            double [NE1], double [NE1], double [NE1],
            int nx, int ny, int ne);
void bicgb (double [NE1], double [NE1], double [NE1], double[NE1], double[NE1],
            double [NE1], double [NE1], int nx, int ny, int ne);
double provv (double [NE1], double [NE1], int ne);
void psorb (double [NE1], double [NE1], double [NE1], double[NE1], double[NE1],
            double [NE1], double [NE1], int nx, int ny, int ne);
void lsorb (double [NE1] ,double [NE1], double [NE1],
            double [NE1], double [NE1], double [NE1], double [NE1],
           int nx, int ny, int ne);
void pdbnd();
void pbnd();
double DX, DY, DT;
double VIS, ALP, BUO;
double RE, PR, GR, TIME, OMG, EPSP;
int  ICYCLE, ITR, IFLG, IRELP, METHOD;
double DMAX;
int ITYPE;
int NITR;
char fname[10][80];
char *fname10;
FILE *in_10, *out_11, *out_12, *out_13, *out_14, *out_21;
FILE *in_15, *in_16, *in_17, *in_18;
double UO[NX1][NY2], UN[NX1][NY2], VO[NX2][NY1], VN[NX2][NY1];
double PO[NX2][NY2], PN[NX2][NY2], TO[NX2][NY2], TN[NX2][NY2];
double DIV[NX1][NY1], PD[NX2][NY2];/* --- SMAC ---*/
void main()
{
     int i,j,NCYCLE;
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
     fscanf(in_10," %lg %lg %d %d ",&DLX,&DLY,&IRELP,&METHOD);
     printf (" ITYPE = %d ICYCLE= %d NITR= %d NCYCLE= %d \n",
             ITYPE,ICYCLE,NITR,NCYCLE);
     printf(" EPSP= %12.3lE  OMG = %12.3lE \n",EPSP,OMG);
     printf(" DT = %12.3lE  RE = %12.3lE  PR = %12.3lE  GR = %12.3lE \n",DT,RE,PR,GR);
     printf(" DLX = %12.3lE  DLY = %12.3lE IRELP = %d  METHOD = %d \n",DLX,DLY,IRELP,METHOD);
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
     else{/*時間進行カウンタがNCYCLEになったら計算終了*/
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
           PD[ix][iy] = 0.0; PN[ix][iy] = 0.0;/* --- SMAC ---*/
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
       fread(PN, sizeof(double), NX2*NY2, in_17);
       fread(PD, sizeof(double), NX2*NY2, in_17);
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
         PO[ix][iy] = PN[ix][iy];/* --- SMAC --- PO:新しい時間ステップでの初期値．PNを保存．*/
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
     velbnd();/* 速度の境界条件の処理 */
}
/*圧力場の計算 Fortran プログラムのコメントを参照してください．*/
void press()
{
     double a1[NE1],a2[NE1],a3[NE1],a4[NE1],a5[NE1],b[NE1],x[NE1];/* --- SMAC --- 圧力補正の線形システム用 */
     int ix,iy,k,nnx,nny,nne;
     /* 線形システム関係の配列の初期化 :*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         k = iy + (ix-1)*NY;
         a1[k]=0.0; a2[k]=0.0; a3[k]=0.0; a4[k]=0.0; a5[k]=0.0;
         b[k]=0.0; x[k]=0.0;
       }
     }
     /*P(ix,iy)の計算*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         /* --- SMAC ---*/
         DIV[ix][iy] = ( UN[ix][iy] - UN[ix-1][iy  ] )/DX
                     + ( VN[ix][iy] - VN[ix  ][iy-1] )/DY;
       }
     }
     /* 係数行列の作成
        計算領域の1点さらに内側の点（仮想セルを含めて考えると2点を除く内側の点）
        に関して係数行列を作成．残りは境界条件を反映させて pdbndc で設定する */
     for (ix=2;ix<=NX-1;ix++){
       for (iy=2;iy<=NY-1;iy++){
         k = iy + (ix-1)*NY;
         a3[k] = DT*( 2.0/(DX*DX) + 2.0/(DY*DY) );
         a4[k] = - 1.0 / (DY*DY) * DT;
         a2[k] = - 1.0 / (DY*DY) * DT;
         a1[k] = - 1.0 / (DX*DX) * DT;
         a5[k] = - 1.0 / (DX*DX) * DT;
         b[k] = -DIV[ix][iy];
         x[k] = PD[ix][iy];
       }
     }
     /* --- SMAC --- 境界条件を係数行列に反映させる */
     pdbndc(a1,a2,a3,a4,a5,b);
     /* --- SMAC --- 圧力補正 P'(PD) に関するポアソン方程式の解法 */
     nnx = NX; nny = NY; nne = NE;
     /* 1. 直接法 : バンドマトリックスによるガウスの消去法 */
     if (METHOD==1) gb    (a1,a2,a3,a4,a5,b,x,nnx,nny,nne);
     /* 2. 反復法1 : point-SOR 法 */
     if (METHOD==2) psorb (a1,a2,a3,a4,a5,b,x,nnx,nny,nne);
     /* 3. 反復法2 : line-SOR 法 : OMG=1で十分．大きくしすぎると発散する */
     if (METHOD==3) lsorb (a1,a2,a3,a4,a5,b,x,nnx,nny,nne);
     /* 4. クリロフ部分空間法1 : 共役残差法 */
     if (METHOD==4) crb   (a1,a2,a3,a4,a5,b,x,nnx,nny,nne);
     /* 5. クリロフ部分空間法2 : BiCGSTAB */
     if (METHOD==5) bicgb (a1,a2,a3,a4,a5,b,x,nnx,nny,nne);
     /* METHOD=4,5における探索ベクトル計算時のゼロ除算対策 */
     if (IFLG==2)   psorb (a1,a2,a3,a4,a5,b,x,nnx,nny,nne);
     /* 線形システム解法で得られた1次元配列の解を2次元配列に置き換える */
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         k=iy+(ix-1)*NY;
         /* 圧力の相対性の処理 : 基準値の設定
            1次独立な解を求める場合は圧力の基準点が強制的に設定される */
         if (IRELP != 2)  PD[ix][iy] = x[k]; /* 圧力の基準点を設けない場合 */
         /* 1次従属な解のうちの1つを求めた後，圧力基準を設ける場合 (IRELP=2) */
         if (IRELP == 2) PD[ix][iy] = x[k]-x[1]; /* P'(1,1)=0 ---> P(1,1)=0 */
       }
     }
     /* 圧力補正の境界条件の処理 */
     pdbnd();
     /* 速度の修正 */
     for (ix=1;ix<=NX-1;ix++){
       for (iy=1;iy<=NY;iy++){
          UN[ix][iy] = UN[ix][iy] + ( PD[ix][iy]-PD[ix+1][iy] )/DX*DT;
       }
     }
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY-1;iy++){
          VN[ix][iy] = VN[ix][iy] + ( PD[ix][iy]-PD[ix][iy+1] )/DY*DT;
       }
     }
     /*新たに得られた速度を用いて境界条件を処理する*/
     velbnd();
     /* 圧力の修正 */
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
          PN[ix][iy] = PO[ix][iy] + PD[ix][iy];
       }
     }
     /* 新たに得られた圧力を用いて境界条件を処理する */
     pbnd();
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
/* 圧力補正の係数行列と境界条件
   本計算プログラムでは，圧力補正の境界条件を，係数行列作成時に考慮する．
   また，解の1次独立性の確保もここで処理する．*/
void pdbndc(double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1], double a5[NE1],
            double b[NE1] )
{
     int ix,iy,k;
     /* --- SMAC --- 係数行列に境界条件を反映させる
     計算領域の左側の係数行列（境界条件を反映）*/
     ix = 1;
     for (iy=2;iy<=NY-1;iy++){
       k = iy + (ix-1)*NY;
       a3[k] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) );
       a4[k] = - 1.0 / (DY*DY) * DT;
       a2[k] = - 1.0 / (DY*DY) * DT;
       a5[k] = - 1.0 / (DX*DX) * DT;
       b[k] = -DIV[ix][iy];
     }
     /* 計算領域の右側の係数行列（境界条件を反映）*/
     ix = NX;
     for (iy=2;iy<=NY-1;iy++){
       k = iy + (ix-1)*NY;
       a3[k] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) );
       a4[k] = - 1.0 / (DY*DY) * DT;
       a2[k] = - 1.0 / (DY*DY) * DT;
       a1[k] = - 1.0 / (DX*DX) * DT;
       b[k] = -DIV[ix][iy];
     }
     /* 計算領域の下側の係数行列（境界条件を反映）*/
     iy = 1;
     for (ix=2;ix<=NX-1;ix++){
       k = iy + (ix-1)*NY;
       a3[k] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) );
       a4[k] = - 1.0 / (DY*DY) * DT;
       a1[k] = - 1.0 / (DX*DX) * DT;
       a5[k] = - 1.0 / (DX*DX) * DT;
       b[k] = -DIV[ix][iy];
     }
     /* 計算領域の上側の係数行列（境界条件を反映）*/
     iy = NY;
     for (ix=2;ix<=NX-1;ix++){
       k = iy + (ix-1)*NY;
       a3[k] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) );
       a2[k] = - 1.0 / (DY*DY) * DT;
       a1[k] = - 1.0 / (DX*DX) * DT;
       a5[k] = - 1.0 / (DX*DX) * DT;
       b[k] = -DIV[ix][iy];
     }
     /* 左下点の係数行列（境界条件を反映）*/
     ix = 1; iy = 1; k = iy + (ix-1)*NY;
     a3[k] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) );
     a4[k] = - 1.0 / (DY*DY) * DT;
     a5[k] = - 1.0 / (DX*DX) * DT;
     b[k] = -DIV[ix][iy];
     /* 左上点の係数行列（境界条件を反映）*/
     ix = 1; iy = NY; k = iy + (ix-1)*NY;
     a3[k] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) );
     a2[k] = - 1.0 / (DY*DY) * DT;
     a5[k] = - 1.0 / (DX*DX) * DT;
     b[k] = -DIV[ix][iy];
     /* 右下点の係数行列（境界条件を反映）*/
     ix = NX; iy = 1; k = iy + (ix-1)*NY;
     a3[k] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) );
     a4[k] = - 1.0 / (DY*DY) * DT;
     a1[k] = - 1.0 / (DX*DX) * DT;
     b[k] = -DIV[ix][iy];
     /* 右上点の係数行列（境界条件を反映）*/
     ix = NX; iy = NY; k = iy + (ix-1)*NY;
     a3[k] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) );
     a2[k] = - 1.0 / (DY*DY) * DT;
     a1[k] = - 1.0 / (DX*DX) * DT;
     b[k] = -DIV[ix][iy];
     /* １次独立な解を得るための処理 : 直接法においては必須
        (ix=1,iy=1 ---> k=1を基準点とし，常にPN[1][1]=PD[1][1]=0とする) */
     if (IRELP==1) {
       a3[1] = 1.0; a4[1] = 0.0; a5[1] = 0.0; b[1] = 0.0;
       /* k=2 の点の処理(k=1とのリンクを断つ) */
       a3[2] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) );
       a4[2] = - 1.0 / (DY*DY) * DT;
       a2[2] = 0.0;
       a5[2] = - 1.0 / (DX*DX) * DT;
       /* k=1+NY の点の処理(k=1とのリンクを断つ) */
       a3[1+NY] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) );
       a4[1+NY] = - 1.0 / (DY*DY) * DT;
       a1[1+NY] = 0.0;
       a5[1+NY] = - 1.0 / (DX*DX) * DT;
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
     int ix,iy;
      for (iy=1;iy<=NY;iy++){/*右側面*/
        PD[NX+1][iy] = PD[NX][iy];
      }
      for (iy=1;iy<=NY;iy++){/*左側面*/
        PD[0][iy] = PD[1][iy];
      }
      for (ix=0;ix<=NX+1;ix++){/*上面*/
        PD[ix][NY+1] = PD[ix][NY];
      }
      for (ix=0;ix<=NX+1;ix++){/*下面*/
        PD[ix][0] = PD[ix][1];
      }
}
/*********************************************************************
*                     圧力の境界条件の処理
*********************************************************************/
void pbnd()
{
     int ix,iy;
      for (iy=1;iy<=NY;iy++){/*右側面*/
        PN[NX+1][iy] = PN[NX][iy];
      }
      for (iy=1;iy<=NY;iy++){/*左側面*/
        PN[0][iy] = PN[1][iy];
      }
      for (ix=0;ix<=NX+1;ix++){/*上面*/
        PN[ix][NY+1] = PN[ix][NY];
      }
      for (ix=0;ix<=NX+1;ix++){/*下面*/
        PN[ix][0] = PN[ix][1];
      }
}
void prout()/*データ出力用*/
{
    fwrite(UN, sizeof(double), NX1*NY2, out_11);
    fwrite(VN, sizeof(double), NX2*NY1, out_12);
    fwrite(PN, sizeof(double), NX2*NY2, out_13);
    fwrite(PD, sizeof(double), NX2*NY2, out_13);
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
/********************************************************************
*                    各種の線形システム解法                         *
*                                                                   *
*  1. ガウスの消去法                                                *
*  2. point-SOR 法                                                  *
*  3. line-SOR 法                                                   *
*  4. 共役残差法                                                    *
*  5. Bi-CGSTAB法                                                   *
*                                                                   *
*  いずれも，2次元のポアソン方程式を5点差分近似にて離散化した       *
*  線形システムを解くためのもので，最適化してある                   *
*  いずれのサブルーチンも同一引数としてある                         *
*                                                                   *
********************************************************************/
/********************************************************************
*                                                                   *
*   2次元ラプラシアン離散化による5点差分近似                        *
*   にて得られた規則的非対称行列Aを含んだ線形システム               *
*                        AX=B                                       *
*   をGaussの消去法を用いて解くサブルーチン．(軸選択無)             *
*   Aはバンドマトリックス                                           *
*                                                                   *
*    線形システム --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    係数行列の計算容量節約：詳細はサブルーチン PRESS を参照        *
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],a5[NE1]           *
*            ---> a[-NY:NY][NE1] -> a[NYNY][NE1]に格納しなおす      *
*                                                                   *
*    B : 既知ベクトル                                               *
*    X : 未知ベクトル ---> これを求める                             *
*                                                                   *
*［変数の説明］                                                     *
*    NX : x方向格子分割数                                           *
*    NY : y方向格子分割数                                           *
*    NE : 総格子点数 = NX * NY                                      *
*                                                                   *
********************************************************************/
void gb (double a1[NE1], double a2[NE1], double a3[NE1],
         double a4[NE1], double a5[NE1], double b[NE1],
         double x[NE1], int nx, int ny, int ne)
{
     double a[NYNY][NE1];
     int ine,i,j,n,k;
     double aa,s;

     /* 直接法のときは(係数行列が特異でなければ)反復なしで，必ず解を得る */
     IFLG = 0;

     /* マトリックスAのゼロクリア */
     for (ine=1;ine<=ne;ine++){
       for (i=-ny;i<=ny;i++){
         a[i+ny][ine] = 0.0;
       }
     }

     /* 必要なところにAT1からAT5までを格納する */
     for (ine=1;ine<=ne;ine++){
       a[-ny+ny][ine] = a1[ine];
       a[ -1+ny][ine] = a2[ine];
       a[  0+ny][ine] = a3[ine];
       a[  1+ny][ine] = a4[ine];
       a[ ny+ny][ine] = a5[ine];
     }

/* 前進消去 */
     for (i=1;i<=ne-1;i++){
       if ( i <= ne-ny ){
         for (j=1;j<=ny;j++){
           aa = a[-j+ny][i+j]/a[0+ny][i];
           b[i+j] = b[i+j] - b[i]*aa;
           n=1;
           for (k=-j+1;k<=ny-j;k++){
             a[k+ny][i+j] = a[k+ny][i+j]-a[n+ny][i]*aa;
             n=n+1;
           }
         }
       }
       else{
         for (j=1;j<=ne-i;j++){
           aa = a[-j+ny][i+j]/a[0+ny][i];
           b[i+j] = b[i+j] - b[i]*aa;
           n=1;
           for (k=-j+1;k<=ne-i-j;k++){
             a[k+ny][i+j] = a[k+ny][i+j]-a[n+ny][i]*aa;
             n=n+1;
           }
         }
       }
     }

/* 後退代入 */
     /* 係数行列の特異性を判定 */
     if ( fabs( a[0+ny][ne] ) <= 1.0e-50 ){
        printf("  Matrix Singular : |A(0,NE)| < 1E-50 \n");
        IFLG = 1; return;/* 圧力補正の計算を発散に設定 */
     }
     x[ne] = b[ne] / a[0+ny][ne];
     for (i=ne-1;i>=1;i--){
       s = 0.0;
       if ( i >  ne-ny ){
         for (n=1;n<=ne-i;n++){
           s = s + a[n+ny][i] * x[i+n];
         }
         x[i] = ( b[i]-s ) / a[0+ny][i];
       }
       else{
         for (n=1;n<=ny;n++){
           s = s + a[n+ny][i] * x[i+n];
         }
         x[i] = ( b[i]-s ) / a[0+ny][i];
       }
     }

}
/********************************************************************
*                                                                   *
*  point-SOR 法による非対称行列 A を含む線形システム解法サブルーチン*
*  2次元ラプラシアン離散化による5点差分近似用                       *
*                                                                   *
*    線形システム --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    係数行列の計算容量節約：詳細はサブルーチン PRESS を参照        *
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],a5[NE1]           *
*                                                                   *
*    B : 既知ベクトル                                               *
*    X : 未知ベクトル ---> これを求める                             *
*                                                                   *
*［変数の説明］                                                     *
*    NX : x方向格子分割数                                           *
*    NY : y方向格子分割数                                           *
*    NE : 総格子点数 = NX * NY                                      *
*    NITR : 許容反復回数(in2d.mac)にて設定                          *
*    EPSP : 収束判定条件で用いる値(in2d.mac)にて設定                *
*                                                                   *
* [収束判定条件]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : 前の反復による値                                 *
*  \vec{x}^{new} : 新しい反復による値                               *
*                                                                   *
********************************************************************/
void psorb (double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1], 
            double a5[NE1], double b[NE1], double x[NE1],
            int nx, int ny, int ne)
{
     int i,j;
     double rnorm1, xold, sum, xnew;

     IFLG=1;/* FLGの初期値は"収束せず" */
 
     for (j=1;j<=NITR;j++){
       rnorm1 = 0.0;
 
       i=1;
         xold = x[i];
         sum =                           +a4[i]*x[i+1]+a5[i]*x[i+ny];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
 
       for (i=2;i<=ny;i++){
         xold = x[i];
         sum =               a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }
 
       for (i=ny+1;i<=ne-ny;i++){
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }
 
       for (i=ne-ny+1;i<=ne-1;i++){
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }
 
       i=ne;
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1];
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
*  2次元ラプラシアン離散化による5点差分近似用                       *
*                                                                   *
*    線形システム --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    係数行列の計算容量節約：詳細はサブルーチン PRESS を参照        *
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],a5[NE1]           *
*                                                                   *
*    B -> bx[NE1] : 既知ベクトル                                    *
*    X -> xn[NE1] : 未知ベクトル ---> これを求める                  *
*                                                                   *
*［変数の説明］                                                     *
*    NX : x方向格子分割数                                           *
*    NY : y方向格子分割数                                           *
*    NE : 総格子点数 = NX * NY                                      *
*    NITR : 許容反復回数(in2d.mac)にて設定                          *
*    EPSP : 収束判定条件で用いる値(in2d.mac)にて設定                *
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
            double at4[NE1], double at5[NE1], double bx[NE1], double xn[NE1],
            int nx, int ny, int ne)
{
     double x1[NE1],xo[NE1],xold[NE1],a[NE1],b[NE1],c[NE1],d[NE1],p[NE1],q[NE1];
     double rnorm;
     int k, inx, iy, ix, iny, j, i ;

     IFLG=1; /* IFLG の初期値は"収束せず" */

     for (k=1;k<=NITR;k++){/* x 軸方向への掃引 : トーマス法による*/
       inx = 1;
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
           xold[iny]=xn[iny];/* x方向へのトーマス法で答えを求める前のXNをXOLDに保存 */
           xo[iny]=xn[iny];/* トーマス法で答えを求める前のXNをXOに保存 */
           inx=inx+1;
         }
       }
       p[1]=c[1]/b[1];/* Ly=b を解く */
       for (j=2;j<=ne-1;j++){
         p[j]=c[j]/( b[j]-a[j]*p[j-1] );
       }
       q[1]=d[1]/b[1];
       for (j=2;j<=ne;j++){
         q[j]=( d[j]-a[j]*q[j-1] ) / ( b[j]-a[j]*p[j-1] );
       }
       x1[ne] = q[ne];/* Ux=y を解く */
       for (j=ne-1;j>=1;j--){
         x1[j] = q[j] - p[j]*x1[j+1];
       }
       inx = 1;
       for (iy=1;iy<=ny;iy++){
         for (ix=1;ix<=nx;ix++){
           iny=iy+(ix-1)*ny;
           /* 得られたX1と反復前のXOにより最新のXNを緩和 */
           xn[iny]=(1.0-OMG)*xo[iny]+OMG*x1[inx];
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
       p[1] = c[1] / b[1];/* Ly=b を解く */
       for (j=2;j<=ne-1;j++){
         p[j] = c[j]/( b[j]-a[j]*p[j-1] );
       }
       q[1] = d[1] / b[1];
       for (j=2;j<=ne;j++){
         q[j] = ( d[j]-a[j]*q[j-1] ) / ( b[j]-a[j]*p[j-1] );
       }
       x1[ne] = q[ne];/* Ux=y を解く */
       for (j=ne-1;j>=1;j--){
         x1[j] = q[j] - p[j]*x1[j+1];
       }
       iny = 1;
       for (ix=1;ix<=nx;ix++){
         for (iy=1;iy<=ny;iy++){
           /* 得られたX1と反復前のXOにより最新のXNを緩和 */
           xn[iny] = (1.0-OMG)*xo[iny]+OMG*x1[iny];
           iny = iny + 1;
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
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],a5[NE1]           *
*                                                                   *
*  B : 既知ベクトル                                                 *
*  X : 未知ベクトル ---> これを求める ---> ここでは便宜上配列xp[NE1]*
*                                                                   *
*［変数の説明］                                                     *
*    NX : x方向格子分割数                                           *
*    NY : y方向格子分割数                                           *
*    NE : 総格子点数 = NX * NY                                      *
*    NITR : 許容反復回数(in2d.mac)にて設定                          *
*    EPSP : 収束判定条件で用いる値(in2d.mac)にて設定                *
*                                                                   *
*［配列の説明］                                                     *
*    r[NE] : r_{k} = B - A x_{k}                                    *
*    p[NE] :  p_{k+1} = r_{k+1} + β_{k} p_{k}, p_{0} = r_{0}       *
*    ap[NE] : A * P                                                 *
*    ar[NE] : A * R                                                 *
*                                                                   *
* [収束判定条件]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : 前の反復による値                                 *
*  \vec{x}^{new} : 新しい反復による値                               *
*                                                                   *
********************************************************************/
void crb (double a1[NE1], double a2[NE1], double a3[NE1],
          double a4[NE1], double a5[NE1], double b[NE1],
          double xp[NE1], int nx, int ny, int ne)
{
     double r[NE1], p[NE1], ap[NE1], ar[NE1], x[NE1], xold[NE1];/* 作業用配列 */
     double rap, apap, alp, rnorm, beta, arap;
     int i, k;

     promv (a1,a2,a3,a4,a5,x,r,nx,ny,ne);/* R に AX を代入 */
     for (i=1;i<=ne;i++){
       r[i]=b[i]-r[i]; p[i]=r[i];/* r_{0}とp_{0}(初期値)の設定 */
       xold[i]=xp[i];/* 前の時刻のXをXOLDに代入 */
     }
     promv (a1,a2,a3,a4,a5,p,ap,nx,ny,ne);/* APに A p_{0} を代入 */
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
       promv (a1,a2,a3,a4,a5,r,ar,nx,ny,ne);/* A r_{k+1} の計算 => AR(NE) */
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
*    ne : 総格子点数(ベクトル A,B のサイズ)                         *
*    c  : a と b の積(スカラー)                                     *
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
*    ne : 総格子点数(正方マトリックス A, ベクトルB,C のサイズ)      *
*    c  : a と b の積(ベクトル)                                     *
*                                                                   *
********************************************************************/
void promv (double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1],
            double a5[NE1], double b[NE1], double c[NE1],
            int nx, int ny, int ne)
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
/********************************************************************
*                                                                   *
*  Bi-CGSTAB 法による非対称行列 A を含む                            *
*  線形システム解法サブルーチン                                     *
*                        AX=B                                       *
*                                                                   *
*    係数行列の計算容量節約：詳細はサブルーチン PRESS を参照        *
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],a5[NE1]           *
*                                                                   *
*    B : 既知ベクトル                                               *
*    X : 未知ベクトル ---> これを求める                             *
*                                                                   *
*［変数の説明］                                                     *
*    NX : x方向格子分割数                                           *
*    NY : y方向格子分割数                                           *
*    NE : 総格子点数 = NX * NY                                      *
*    NITR : 許容反復回数(in2d.mac)にて設定                          *
*    EPSP : 収束判定条件で用いる値(in2d.mac)にて設定                *
*                                                                   *
*［配列の説明］                                                     *
*    T[NE] : t_{k} = r_{k} - α_{k} A p_{k}                         *
*    X[NE] : x_{k+1} = x_{k} + α_{k} p_{k} + ξ_{K} t_{k}          *
*    R[NE] : r_{k+1} = t_{k} - ξ_{k} A t_{k}                       *
*            r_{0} = B - A x_{0}                                    *
*    P[NE] : p_{k+1} = r_{k+1} + β_{k} ( p_{k}-ξ_{k} A p_{k} )    *
*            p_{0} = r_{0}                                          *
*    AP[NE] : A * P                                                 *
*    AT[NE] : A * T                                                 *
*                                                                   *
* [収束判定条件]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : 前の反復による値                                 *
*  \vec{x}^{new} : 新しい反復による値                               *
*                                                                   *
********************************************************************/
void bicgb ( double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1], double a5[NE1],
             double b[NE1], double x[NE1], int nx, int ny, int ne)
{
     double r[NE1], ap[NE1], at[NE1], p[NE1], s[NE1], t[NE1], xold[NE1];
     int i,j;
     double sr1, sap, alpha, att, atat, xi, sr2, beta, rnorm;

     for (i=1;i<=ne;i++){
       xold[i] = x[i];
     }

      promv (a1,a2,a3,a4,a5,x,r,nx,ny,ne);/* R に AX を代入 */
      for (i=1;i<=ne;i++){/* r_{0}とp_{0}(初期値)，そして s=r_{0} の設定 */
        r[i] = b[i] - r[i]; p[i] = r[i]; s[i] = r[i];
      }

      /* 繰り返し計算 */
      for (j=1;j<=NITR;j++){
        sr1=provv (s,r,ne);/* ( s, r_{k} ) の計算 => SR1 */
        promv (a1,a2,a3,a4,a5,p,ap,nx,ny,ne);/* A p_{k} の計算 => AP(NE) */
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
        promv (a1,a2,a3,a4,a5,t,at,nx,ny,ne);/* A t_{k} の計算 => AT(NE) */
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
