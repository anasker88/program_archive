/********************************************************************
* ファイル名  ：simple3d.C                                          *
* タイトル    ：SIMPLE法による3次元熱流動解析プログラム             *
* 製作者      ：平野　博之                                          *
* 所属        ：岡山理科大学 工学部 バイオ・応用化学科              *
* 製作日      ：2011.11.01                                          *
* 言語        ：C                                                   *
*********************************************************************
* プログラムを実行すると，"in3d.sim"(必ず小文字）を読みに行く．     *
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
#define NEU 7600 /* NEU:Uの未知数の数 */
#define NEV 7600 /* NEV:Vの未知数の数 */
#define NEW 7600 /* NEW:Wの未知数の数 */
#define NEP 8000 /* NEP:P, PD, Tの未知数の数 */
#define NX1 21 /* NX1=NX+1 */
#define NX2 22 /* NX2=NX+2 */
#define NY1 21 /* NY1=NY+1 */
#define NY2 22 /* NY2=NY+2 */
#define NZ1 21 /* NZ1=NZ+1 */
#define NZ2 22 /* NZ2=NZ+2 */
#define NEU1 7601 /* NEU1=NEU+1 */
#define NEV1 7601 /* NEV1=NEV+1 */
#define NEW1 7601 /* NEW1=NEW+1 */
#define NEP1 8001 /* NEP1=NEP+1 */
#define NXY 400 /* NXY=NX*NY */
#define NXYNXY 801 /* NXYNXY=2*NXY+1 ガウスの消去法で使用 */
void cinit();
void adv();
void calvel();
void press();
void caltem();
void ubnd(double [NX1][NY2][NZ2]);
void vbnd(double [NX2][NY1][NZ2]);
void wbnd(double [NX2][NY2][NZ1]);
void tbnd(double [NX2][NY2][NZ2]);
void pdbnd(double [NX2][NY2][NZ2]);
void psbnd(double [NX2][NY2][NZ2]);
void prout();
void tecplt( char *fname_tec);
void gb (double [NEP1], double [NEP1], double [NEP1], double [NEP1], double [NEP1],
         double [NEP1], double [NEP1], double [NEP1], double [NEP1],    int nx, int ny, int nz, int ne, int nxy);
void psorb (double [NEP1], double [NEP1], double [NEP1], double [NEP1], double [NEP1],
            double [NEP1], double [NEP1], double [NEP1], double [NEP1], int nx, int ny, int nz, int ne, int nxy);
void lsorb (double [NEP1] ,double [NEP1], double [NEP1], double [NEP1], double [NEP1],
            double [NEP1], double [NEP1], double [NEP1], double [NEP1], int nx, int ny, int nz, int ne, int nxy);
void crb   (double [NEP1], double [NEP1], double [NEP1], double [NEP1], double [NEP1],
            double [NEP1], double [NEP1], double [NEP1], double [NEP1], int nx, int ny, int nz, int ne, int nxy);
void promv (double [NEP1], double [NEP1], double [NEP1], double [NEP1], double [NEP1],
            double [NEP1], double [NEP1], double [NEP1], double [NEP1], int nx, int ny, int nz, int ne, int nxy);
double provv (double [NEP1], double [NEP1], int ne);
void bicgb (double [NEP1], double [NEP1], double [NEP1], double [NEP1], double [NEP1],
            double [NEP1], double [NEP1], double [NEP1], double [NEP1], int nx, int ny, int nz, int ne, int nxy);
double Max(double x,  double y);
double DX, DY, DZ, DT;
double VIS, ALP, BUO;
double RE, PR, GR, TIME, OMG, EPSP, EPSC, ALPHAP, ALPHAU, ALPHAV, ALPHAW, ALPHAT;
double DMAX, RO;
int ICYCLE, ITR, IFLG, IFLGC, IRELP, METHOD, NMT, ITYPE;
int NITR, NCNT;
char fname[12][80];
char *fname10;
FILE *in_10, *out_11, *out_12, *out_13, *out_14, *out_15, *out_21;
FILE *in_16, *in_17, *in_18, *in_19, *in_20;
double AU3[NX1][NY2][NZ2], AV3[NX2][NY1][NZ2], AW3[NX2][NY2][NZ1];
double  U[NX1][NY2][NZ2],  V[NX2][NY1][NZ2],  W[NX2][NY2][NZ1];
double US[NX1][NY2][NZ2], VS[NX2][NY1][NZ2], WS[NX2][NY2][NZ1];
double PS[NX2][NY2][NZ2], PD[NX2][NY2][NZ2], T[NX2][NY2][NZ2];
double UO[NX1][NY2][NZ2], VO[NX2][NY1][NZ2], WO[NX2][NY2][NZ1], TO[NX2][NY2][NZ2], XX[NEP1];
void main()
{
     int i,j,NCYCLE,ITRCNT,ITRALL,ICNTAL;
     double DLX,DLY,DLZ;
     char buff[80];
     ITRALL=0; ICNTAL=0;
     fname10="in3d.sim";
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
     fgets(buff, sizeof buff,in_10);/*in3d.sim中のコメント行(12行目)のスキップ*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(13行目)*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(14行目)*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(15行目)*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(16行目)*/
     fgets(buff, sizeof buff,in_10);/*                   ↓(17行目)*/
     fscanf(in_10," %d %d %d %d %d ",&ITYPE,&ICYCLE,&NITR,&NCYCLE,&NCNT);
     fgets(buff, sizeof buff,in_10);/*in3d.mac中のコメント行(19行目)のスキップ*/
     fgets(buff, sizeof buff,in_10);/*in3d.mac中のコメント行(20行目)のスキップ*/
     fscanf(in_10," %lg %lg %lg %lg %lg %lg %lg %lg ",&EPSP,&EPSC,&OMG,&ALPHAP,&ALPHAU,&ALPHAV,&ALPHAW,&ALPHAT);
     fgets(buff, sizeof buff,in_10);/*in3d.mac中のコメント行(22行目)のスキップ*/
     fgets(buff, sizeof buff,in_10);/*in3d.mac中のコメント行(23行目)のスキップ*/
     fscanf(in_10," %lg %lg %lg %lg ",&DT,&RE,&PR,&GR);
     fgets(buff, sizeof buff,in_10);/*in2d.mac中のコメント行(25行目)のスキップ*/
     fgets(buff, sizeof buff,in_10);/*in2d.mac中のコメント行(26行目)のスキップ*/
     fscanf(in_10," %lg %lg %lg %d %d %d ",&DLX,&DLY,&DLZ,&IRELP,&METHOD,&NMT);
     printf (" ITYPE= %d ICYCLE= %d NITR= %d NCYCLE= %d NCNT= %d \n",ITYPE,ICYCLE,NITR,NCYCLE,NCNT);
     printf(" EPSP= %12.3lE EPSC= %12.3lE OMG= %12.3lE \n",EPSP,EPSC,OMG);
     printf(" ALP_P= %12.3lE ALP_U= %12.3lE ALP_V= %12.3lE ALP_W= %12.3lE ALP_T= %12.3lE \n",ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT);
     printf(" DT= %12.3lE  RE= %12.3lE  PR= %12.3lE  GR= %12.3lE \n",DT,RE,PR,GR);
     printf(" DLX= %12.3lE  DLY= %12.3lE DLZ= %12.3lE IRELP= %d  METHOD= %d \n",DLX,DLY,DLZ,IRELP,METHOD);
     if ( ICYCLE != 0 ){/*継続の計算の場合*/
       in_16=fopen(fname[6],"rb");/*Uデータファイルのオープン*/
       in_17=fopen(fname[7],"rb");/*Vデータファイルのオープン*/
       in_18=fopen(fname[8],"rb");/*Wデータファイルのオープン*/
       in_19=fopen(fname[9],"rb");/*P,PDデータファイルのオープン*/
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
     /*文献[1]にしたがった離散化を考慮してROを定義
       本書では無次元化をおこなってから計算しており，RO=1とする*/
     RO = 1.0;
     cinit();/*初期値の設定*/
     label_750:{};/*時間進行のための戻り点*/
     adv();/*時間進行*/
     ITRCNT=0;/*連続の式の収束判定の反復回数を表す変数*/
     label_760:{};/*速度・圧力・温度補正のための戻り点*/
     ITRCNT=ITRCNT+1;
     /*速度，圧力，温度の線形システム解法が収束したかどうかのパラメータ(IFLG)を初期化
       これは，クリロフ部分空間法を含む反復法で用いる
       IFLG -> 0:収束 1:発散(設定された許容回数NITR以下で解が得られない)*/
     IFLG=0;
     /*速度，圧力，温度を同時に補正しながら連続の式を満たすまで計算するが，
       この際，連続の式を満たすかどうかを判定する反復回数を意味する
       パラメータ(IFLGC)を初期化
       IFLGC-> 0:連続の式を満たす 1:発散(許容回数NCNT以下で解が得られない)*/
     IFLGC=0;
     /* --- 注意 ---
      本計算においては，ICYCLE=0のときの初期条件は速度場と圧力場をゼロと
      しているので，これらに関わる線形システムを解くのが困難となるため，
      ICYCLE=1のときだけは温度の計算にジャンプするようにする．
      計算条件や問題に応じて適宜削除あるいは変更*/
     if ( ICYCLE==1 ){
       goto label_770;
     }
     calvel();/*速度場の計算*/
     if ( IFLG==1 ) {/*速度計算の線形システム解法が収束していないとき*/
       printf(" NOT CONVERGE ! \n");
       prout(); /*データを出力して強制終了*/
       goto label_900;
     }
     press();/*圧力補正の計算と，圧力場と速度場の更新*/
     ITRALL = ITRALL + ITR;/*圧力計算の線形システム解法の総反復回数*/
     label_770:{};/* 速度場と圧力場がゼロとなるときのスキップ先*/
     if ( IFLG==1 ) {/*圧力計算の線形システム解法が収束していないとき*/
       printf(" NOT CONVERGE ! \n");
       prout(); /*データを出力して強制終了*/
       goto label_900;
     }
     if ( ITYPE==2 ){/*非等温場計算の場合*/
       caltem();/*温度場を計算*/
     }
     if ( IFLG==1 ) {/*温度計算の線形システム解法が収束していないとき*/
       printf(" NOT CONVERGE ! \n");
       prout(); /*データを出力して強制終了*/
       goto label_900;
     }
     if ( IFLGC==1 ){/*速度，圧力，温度場の計算が収束して連続の式をEPSP以下で満足していないとき(IFLGC=1)*/
       if ( ITRCNT < NCNT ) {/*予め設定された許容回数NCNTに達していないときは再計算*/
         goto label_760;
       }
       else{/*予め設定された許容回数NCNTに達したらデータを出力して強制終了*/
         printf(" NOT CONVERGE ! \n");
         prout(); /*データを出力して強制終了*/
         goto label_900;
       }
     }
     /*速度，圧力，温度場の計算が収束して連続の式をEPSP以下で満足し(IFLGC=0)
       時間進行カウンタ(ICYCLE)がNCYCLEより小さい時*/
     ICNTAL = ICNTAL + ITRCNT;
     if ( ICYCLE < NCYCLE ){/*時間進行カウンタ(ICYCLE)がNCYCLEより小さい時*/
       goto label_750;
     }
     else{/*時間進行カウンタがNCYCLEになったら計算終了*/
       prout();
     }
     label_900:{};
     tecplt(fname[11]);/*Tecplot用データの出力*/
     fclose (in_10); fclose (out_11); fclose (out_12);
     fclose (out_13); fclose (out_14); fclose(out_15);
     /* 総反復回数を出力 */
     printf ("  Total Iteration Number \n");
     printf ("  To solve Continuity Eq. %8d \n",ICNTAL);
     printf ("  To solve Linear system of PD %8d \n",ITRALL);
}
/* ２つの数の大きい方の値を計算する関数 */
double Max(double x,  double y)
{
    if (x < y){
        return y;
    }
    else {
        return x;
    }
}
/*初期設定*/
void cinit()
{
     int ix,iy,iz;
     if ( ICYCLE == 0 ){/*新規計算の場合*/
       for (ix=0;ix<=NX;ix++){/*Uの初期値設定*/
         for (iy=0;iy<=NY+1;iy++){
           for (iz=0;iz<=NZ+1;iz++){
             U[ix][iy][iz] = 0.0;
           }
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*Vの初期値設定*/
         for (iy=0;iy<=NY;iy++){
           for (iz=0;iz<=NZ+1;iz++){
             V[ix][iy][iz] = 0.0;
           }
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*Wの初期値設定*/
         for (iy=0;iy<=NY+1;iy++){
           for (iz=0;iz<=NZ;iz++){
             W[ix][iy][iz] = 0.0;
           }
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*PとP'の初期値設定*/
         for (iy=0;iy<=NY+1;iy++){
           for (iz=0;iz<=NZ+1;iz++){
             PS[ix][iy][iz] = 0.0; PD[ix][iy][iz] = 0.0;
           }
         }
       }
       for (ix=0;ix<=NX+1;ix++){/* Tの初期値設定(領域内は高温(+0.5)と低温(-0.5)の中間温度) */
         for (iy=0;iy<=NY+1;iy++){/*（注意）浮力項の計算で温度の配列を使用しているので*/
           for (iz=0;iz<=NZ+1;iz++){/*      等温場でもT=0として初期条件だけは設定する*/
             T[ix][iy][iz]=0.0;/*           必要がある．ゼロ以外の値を入れると浮力項が計算される可能性があるので注意．*/
           }
         }
       }
       for (iy=0;iy<=NY+1;iy++){/*Tの境界：右側壁（冷却）T=-0.5*/
         for (iz=0;iz<=NZ+1;iz++){
           T[NX+1][iy][iz] = 2.0 * ( -0.5 ) - T[NX][iy][iz];
         }
       }
       for (iy=0;iy<=NY+1;iy++){/*Tの境界：左側壁（加熱）T=+0.5*/
         for (iz=0;iz<=NZ+1;iz++){
           T[0][iy][iz] = 2.0 * ( +0.5 ) - T[1][iy][iz];
         }
       }
       for (ix=1;ix<=NX;ix++){/*Tの境界：上面（断熱）*/
         for (iz=1;iz<=NZ+1;iz++){
           T[ix][NY+1][iz] = T[ix][NY][iz];
         }
       }
       for (ix=1;ix<=NX;ix++){/*Tの境界：下面（断熱）*/
         for (iz=1;iz<=NZ+1;iz++){
           T[ix][0][iz] = T[ix][1][iz];
       }
     }
       for (ix=1;ix<=NX;ix++){/*Tの境界：後面（断熱）*/
         for (iy=1;iy<=NY;iy++){
           T[ix][iy][NZ+1] = T[ix][iy][NZ];
         }
       }
       for (ix=1;ix<=NX;ix++){/*Tの境界：前面（断熱）*/
         for (iy=1;iy<=NY;iy++){
           T[ix][iy][0] = T[ix][iy][1];
         }
       }
     }
     else{/*継続計算（すでにある計算結果からスタート）の場合*/
       /*Uデータファイルからの読み込み[Unit No.=16]*/
       fread(U, sizeof(double), NX1*NY2*NZ2, in_16);
       /*Vデータファイルからの読み込み[Unit No.=17]*/
       fread(V, sizeof(double), NX2*NY1*NZ2, in_17);
       /*Wデータファイルからの読み込み[Unit No.=18]*/
       fread(W, sizeof(double), NX2*NY2*NZ1, in_18);
       /*Pデータファイルからの読み込み[Unit No.=19]*/
       fread(PS, sizeof(double), NX2*NY2*NZ2, in_19);
       fread(PD, sizeof(double), NX2*NY2*NZ2, in_19);
       /*Tデータファイルからの読み込み[Unit No.=20]*/
       fread(T , sizeof(double), NX2*NY2*NZ2, in_20);
       fclose (in_16); fclose (in_17); fclose (in_18); fclose (in_19); fclose(in_20);
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
     /* U -> UO : 必要なら入れ替える前にUとUOから変動量を求める
        U       : 前の時間ステップにおいて最終的に得られた値
                  速度・圧力・温度同時補正の度に更新される
        UO      : 前の時間ステップの値を保存しておく．*/
     for (ix=0;ix<=NX;ix++){
       for (iy=0;iy<=NY+1;iy++){
         for (iz=0;iz<=NZ+1;iz++){
           UO[ix][iy][iz] = U[ix][iy][iz];
         }
       }
     }
     /* V -> VO : 必要なら入れ替える前にVとVOから変動量を求める
        V       : 前の時間ステップにおいて最終的に得られた値
                  速度・圧力・温度同時補正の度に更新される
        VO      : 前の時間ステップの値を保存しておく．*/
     for (ix=0;ix<=NX+1;ix++){
       for (iy=0;iy<=NY;iy++){
         for (iz=0;iz<=NZ+1;iz++){
           VO[ix][iy][iz] = V[ix][iy][iz];
         }
       }
     }
     /* W -> WO : 必要なら入れ替える前にWとWOから変動量を求める
        W       : 前の時間ステップにおいて最終的に得られた値
                  速度・圧力・温度同時補正の度に更新される
        WO      : 前の時間ステップの値を保存しておく．*/
     for (ix=0;ix<=NX+1;ix++){
       for (iy=0;iy<=NY+1;iy++){
         for (iz=0;iz<=NZ;iz++){
           WO[ix][iy][iz] = W[ix][iy][iz];
         }
       }
     }
     /* T -> TO : 必要なら入れ替える前にTとTOから変動量を求める
        T       : 前の時間ステップにおいて最終的に得られた値
                  速度・圧力・温度同時補正の繰り返しでTが更新される
        TO      : 前の時間ステップでの値を保存しておく．*/
     for (ix=0;ix<=NX+1;ix++){
       for (iy=0;iy<=NY+1;iy++){
         for (iz=0;iz<=NZ+1;iz++){
           TO[ix][iy][iz] = T[ix][iy][iz];
         }
       }
     }
}
/*仮の速度場 ( U*, V*, W* ) の計算*/
void calvel()
{
     double atp1[NEP1],atp2[NEP1],atp3[NEP1],atp4[NEP1],atp5[NEP1],atp6[NEP1],atp7[NEP1],bp[NEP1],xx[NEP1];
     int ix,iy,iz,i,k,nnx,nny,nnz,nnxy,nne;
     double du1,du2,du4,du5,du6,du7;
     double dv1,dv2,dv4,dv5,dv6,dv7;
     double dw1,dw2,dw4,dw5,dw6,dw7;
     double uu1,vu2,vu4,uu5,wu6,wu7;
     double uv1,vv2,vv4,uv5,wv6,wv7;
     double uw1,vw2,vw4,uw5,ww6,ww7;
     double fu1,fu2,fu4,fu5,fu6,fu7;
     double fv1,fv2,fv4,fv5,fv6,fv7;
     double fw1,fw2,fw4,fw5,fw6,fw7;
     double peu1,peu2,peu4,peu5,peu6,peu7;
     double pev1,pev2,pev4,pev5,pev6,pev7;
     double pew1,pew2,pew4,pew5,pew6,pew7;
     double apeu1,apeu2,apeu4,apeu5,apeu6,apeu7;
     double apev1,apev2,apev4,apev5,apev6,apev7;
     double apew1,apew2,apew4,apew5,apew6,apew7;
     double au1,au2,au4,au5,au6,au7;
     double av1,av2,av4,av5,av6,av7;
     double aw1,aw2,aw4,aw5,aw6,aw7;
     double scu,scv,scw,spu,spv,spw,tv;
     /* ****************** */
     /* U*(ix,iy,iz)の計算 */
     /* ****************** */
     i=1;
     for (iz=1;iz<=NZ;iz++){
       for (ix=1;ix<=NX-1;ix++){
         for (iy=1;iy<=NY;iy++){
         /* Uに関するDの計算 */
         du5=VIS*DY*DZ/DX; du1=VIS*DY*DZ/DX; du4=VIS*DX*DZ/DY; du2=VIS*DX*DZ/DY;
         du7=VIS*DX*DY/DZ; du6=VIS*DX*DY/DZ;
         /* 速度・圧力・温度場の繰り返し計算となっても線形システムの係数となる
            ここの速度は更新されないように(U,V,W)には前の時刻の収束値(UO,VO,WO)を用いる*/
         uu5=( UO[ix  ][iy  ][iz  ] +UO[ix+1][iy  ][iz  ] ) / 2.0;
         uu1=( UO[ix-1][iy  ][iz  ] +UO[ix  ][iy  ][iz  ] ) / 2.0;
         vu4=( VO[ix  ][iy  ][iz  ] +VO[ix+1][iy  ][iz  ] ) / 2.0;
         vu2=( VO[ix  ][iy-1][iz  ] +VO[ix+1][iy-1][iz  ] ) / 2.0;
         wu7=( WO[ix  ][iy  ][iz  ] +WO[ix+1][iy  ][iz  ] ) / 2.0;
         wu6=( WO[ix  ][iy  ][iz-1] +WO[ix+1][iy  ][iz-1] ) / 2.0;
         fu5=RO*uu5*DY*DZ; fu1=RO*uu1*DY*DZ; fu4=RO*vu4*DX*DZ; fu2=RO*vu2*DX*DZ;
         fu7=RO*wu7*DX*DY; fu6=RO*wu6*DX*DY;
         /* Uに関するPeの計算 */
         peu5=fu5/du5; peu1=fu1/du1; peu4=fu4/du4; peu2=fu2/du2;
         peu7=fu7/du7; peu6=fu6/du6;
         /* 関数A(|Pe|)の計算 */
         if ( NMT==1 ) {/* NMT=1なら 風上法 */
           apeu5=1.0; apeu1=1.0; apeu4=1.0; apeu2=1.0; apeu7=1.0; apeu6=1.0;
         }
         if ( NMT==2 ) {/* NMT=2なら べき乗法 */
           apeu5=Max( 0.0, pow( (1.0-0.1*fabs(peu5)),5 ) );
           apeu1=Max( 0.0, pow( (1.0-0.1*fabs(peu1)),5 ) );
           apeu4=Max( 0.0, pow( (1.0-0.1*fabs(peu4)),5 ) );
           apeu2=Max( 0.0, pow( (1.0-0.1*fabs(peu2)),5 ) );
           apeu7=Max( 0.0, pow( (1.0-0.1*fabs(peu7)),5 ) );
           apeu6=Max( 0.0, pow( (1.0-0.1*fabs(peu6)),5 ) );
         }
         /* U* の線形システムの係数の計算 */
         au5=( du5*apeu5 + Max(-fu5,0.0) );
         au1=( du1*apeu1 + Max( fu1,0.0) );
         au4=( du4*apeu4 + Max(-fu4,0.0) );
         au2=( du2*apeu2 + Max( fu2,0.0) );
         au7=( du7*apeu7 + Max(-fu7,0.0) );
         au6=( du6*apeu6 + Max( fu6,0.0) );
         scu=0.0; spu=0.0;/* 生成項 */
         /* U* の線形システムのBの計算 */
         bp[i]=scu*DX*DY*DZ+( PS[ix][iy][iz]-PS[ix+1][iy][iz] )*DY*DZ+RO*DX*DY*DZ/DT*UO[ix][iy][iz];
         /* 境界の処理
          1. 境界のU_BND についても変数とみなしてAU3[ix][iy][iz]を求める
          2. その後，境界の項 au(1,2,4,5,6,7)*U_BND をBに吸収させる
          3. 2に対応して(吸収の後)，境界のau(1,2,4,5,6,7)をゼロとしてリンクを断つ */
         AU3[ix][iy][iz] = au1+au2+au4+au5+au6+au7 + RO*DX*DY*DZ/DT - spu*DX*DY*DZ;
         if ( ix==NX-1 ) {/* U（右面）*/
           bp[i]=bp[i]+au5*U[NX][iy][iz]; au5=0.0;
         }
         if ( ix==1 ) {/* U（左面）*/
           bp[i]=bp[i]+au1*U[0][iy][iz]; au1=0.0;
         }
         if ( iy==NY ) {/* U（上面）*/
           bp[i]=bp[i]+au4*U[ix][NY+1][iz]; au4=0.0;
         }
         if ( iy==1 ) {/* U（下面）*/
           bp[i]=bp[i]+au2*U[ix][0][iz]; au2=0.0;
         }
         if ( iz==NZ ) {/* U（後面）*/
           bp[i]=bp[i]+au7*U[ix][iy][NZ+1]; au7=0.0;
         }
         if ( iz==1 ) {/* U（前面）*/
           bp[i]=bp[i]+au6*U[ix][iy][0]; au6=0.0;
         }
         /* 線形システム計算用に係数行列を計算
            ALPHAU : U* 計算のための緩和係数 */
         atp3[i]=AU3[ix][iy][iz]/ALPHAU;
         atp5[i]=-au5; atp1[i]=-au1; atp4[i]=-au4; atp2[i]=-au2;
         atp7[i]=-au7; atp6[i]=-au6;
         bp[i]=bp[i]+(1.0-ALPHAU)*atp3[i]*UO[ix][iy][iz];
         xx[i]=U[ix][iy][iz];/* 線形システム解法前のxx[i]の初期化(任意) */
         i=i+1;
         }
       }
     }
     /* --- SIMPLE --- 仮の速度場U*(US) に関する線形システムの解法 */
     nnx=NX-1; nny=NY; nnz=NZ; nne=NEU; nnxy=nnx*nny;
     /* 1. 直接法 : バンドマトリックスによるガウスの消去法 */
     if (METHOD==1) gb    (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 2. 反復法1 : point-SOR 法 */
     if (METHOD==2) psorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 3. 反復法2 : line-SOR 法 : OMG=1で十分．大きくしすぎると発散する */
     if (METHOD==3) lsorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 4. クリロフ部分空間法1 : 共役残差法 */
     if (METHOD==4) crb   (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 5. クリロフ部分空間法2 : BiCGSTAB */
     if (METHOD==5) bicgb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* METHOD=4,5における探索ベクトル計算時のゼロ除算対策 */
     if (IFLG==2)   psorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 線形システム解法で得られた1次元配列の解を3次元配列に置き換える */
     for (ix=1;ix<=NX-1;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
           k = iy + (ix-1)*nny + (iz-1)*(nnx*nny);
           US[ix][iy][iz] = xx[k];
         }
       }
     }
     ubnd (US);/* U*の境界条件の処理 */
     /* ****************** */
     /* V*(ix,iy,iz)の計算 */
     /* ****************** */
     i=1;
     for (iz=1;iz<=NZ;iz++){
       for (ix=1;ix<=NX;ix++){
         for (iy=1;iy<=NY-1;iy++){
         /* Vに関するDの計算 */
         dv5 = VIS*DY*DZ/DX; dv1 = VIS*DY*DZ/DX;
         dv4 = VIS*DX*DZ/DY; dv2 = VIS*DX*DZ/DY;
         dv7 = VIS*DX*DY/DZ; dv6 = VIS*DX*DY/DZ;
         /* 速度・圧力・温度場の繰り返し計算となっても線形システムの係数となる
            ここの速度は更新されないように(U,V,W)には前の時刻の収束値(UO,VO,WO)を用いる*/
         uv5 = ( UO[ix  ][iy  ][iz  ] +UO[ix  ][iy+1][iz  ] ) / 2.0;
         uv1 = ( UO[ix-1][iy  ][iz  ] +UO[ix-1][iy+1][iz  ] ) / 2.0;
         vv4 = ( VO[ix  ][iy  ][iz  ] +VO[ix  ][iy+1][iz  ] ) / 2.0;
         vv2 = ( VO[ix  ][iy-1][iz  ] +VO[ix  ][iy  ][iz  ] ) / 2.0;
         wv7 = ( WO[ix  ][iy  ][iz  ] +WO[ix  ][iy+1][iz  ] ) / 2.0;
         wv6 = ( WO[ix  ][iy  ][iz-1] +WO[ix  ][iy+1][iz-1] ) / 2.0;
         fv5 = RO*uv5*DY*DZ; fv1 = RO*uv1*DY*DZ;
         fv4 = RO*vv4*DX*DZ; fv2 = RO*vv2*DX*DZ;
         fv7 = RO*wv7*DX*DY; fv6 = RO*wv6*DX*DY;
         /* Vに関するPeの計算 */
         pev5 = fv5/dv5; pev1 = fv1/dv1;
         pev4 = fv4/dv4; pev2 = fv2/dv2;
         pev7 = fv7/dv7; pev6 = fv6/dv6;
         /* 関数A(|Pe|)の計算 */
         if ( NMT==1 ) {/* NMT=1なら 風上法 */
           apev5=1.0; apev1=1.0 ;apev4=1.0; apev2=1.0; apev7=1.0; apev6=1.0;
         }
         if ( NMT==2 ) {/* NMT=2なら べき乗法 */
           apev5 = Max( 0.0, pow( (1.0-0.1*fabs(pev5)),5 ) );
           apev1 = Max( 0.0, pow( (1.0-0.1*fabs(pev1)),5 ) );
           apev4 = Max( 0.0, pow( (1.0-0.1*fabs(pev4)),5 ) );
           apev2 = Max( 0.0, pow( (1.0-0.1*fabs(pev2)),5 ) );
           apev7 = Max( 0.0, pow( (1.0-0.1*fabs(pev7)),5 ) );
           apev6 = Max( 0.0, pow( (1.0-0.1*fabs(pev6)),5 ) );
         }
         /* V* の線形システムの係数の計算 */
         av5 = ( dv5*apev5 + Max(-fv5,0.0) );
         av1 = ( dv1*apev1 + Max( fv1,0.0) );
         av4 = ( dv4*apev4 + Max(-fv4,0.0) );
         av2 = ( dv2*apev2 + Max( fv2,0.0) );
         av7 = ( dv7*apev7 + Max(-fv7,0.0) );
         av6 = ( dv6*apev6 + Max( fv6,0.0) );
         /* 速度・圧力・温度場の繰り返し計算となった場合，更新されたT
            により浮力も更新される */
         tv = ( T[ix][iy][iz] + T[ix][iy+1][iz] ) / 2.0;
         scv = BUO * tv; spv=0.0;/* 生成項(浮力項)の計算 */
         /* V* の線形システムの B の計算 */
         bp[i] = scv*DX*DY*DZ+( PS[ix][iy][iz]-PS[ix][iy+1][iz] )*DX*DZ
               + RO*DX*DY*DZ/DT*VO[ix][iy][iz];
         /* 境界の処理
          1. 境界のV_BND についても変数とみなしてAV3[ix][iy][iz]を求める
          2. その後，境界の項 av(1,2,4,5,6,6)*V_BND をBに吸収させる
          3. 2に対応して(吸収の後)，境界のav(1,2,4,5,6,7)をゼロとしてリンクを断つ */
         AV3[ix][iy][iz] = av1+av2+av4+av5+av6+av7 + RO*DX*DY*DZ/DT - spv*DX*DY*DZ;
         if ( ix==NX ) {/* V（右面）*/
           bp[i] = bp[i]+av5*V[NX+1][iy][iz]; av5 = 0.0;
         }
         if ( ix==1 ) {/* V（左面）*/
           bp[i] = bp[i]+av1*V[0][iy][iz]; av1 = 0.0;
         }
         if ( iy==NY-1 ) {/* V（上面）*/
           bp[i] = bp[i] + av4*V[ix][NY][iz]; av4 = 0.0;
         }
         if ( iy==1 ) {/* V（下面）*/
           bp[i] = bp[i] + av2*V[ix][0][iz]; av2 = 0.0;
         }
         if ( iz==NZ ) {/* V（後面）*/
           bp[i] = bp[i] + av7*V[ix][iy][NZ+1]; av7 = 0.0;
         }
         if ( iz==1 ) {/* V（前面）*/
           bp[i] = bp[i] + av6*V[ix][iy][0]; av6 = 0.0;
         }
         /* 線形システム計算用に係数行列を計算
            ALPHAV : V* 計算のための緩和係数 */
         atp3[i] = AV3[ix][iy][iz] / ALPHAV;
         atp5[i] = -av5; atp1[i] = -av1;
         atp4[i] = -av4; atp2[i] = -av2;
         atp7[i] = -av7; atp6[i] = -av6;
         bp[i] = bp[i] + (1.0-ALPHAV)*atp3[i]*VO[ix][iy][iz];
         xx[i]=V[ix][iy][iz];/* 線形システム解法前のxx[i]の初期化 */
         i = i + 1;
         }
       }
     }
     /* --- SIMPLE --- 仮の速度場V*(VS) に関する線形システムの解法 */
     nnx = NX; nny = NY-1; nnz = NZ; nne = NEV; nnxy = nnx*nny;
     /* 1. 直接法 : バンドマトリックスによるガウスの消去法 */
     if (METHOD==1) gb    (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 2. 反復法1 : point-SOR 法 */
     if (METHOD==2) psorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 3. 反復法2 : line-SOR 法 : OMG=1で十分．大きくしすぎると発散する */
     if (METHOD==3) lsorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 4. クリロフ部分空間法1 : 共役残差法 */
     if (METHOD==4) crb   (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 5. クリロフ部分空間法2 : BiCGSTAB */
     if (METHOD==5) bicgb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* METHOD=4,5における探索ベクトル計算時のゼロ除算対策 */
     if (IFLG==2)   psorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 線形システム解法で得られた1次元配列の解を3次元配列に置き換える */
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY-1;iy++){
         for (iz=1;iz<=NZ;iz++){
           k = iy + (ix-1)*nny + (iz-1)*(nnx*nny);
           VS[ix][iy][iz] = xx[k];
         }
       }
     }
     vbnd (VS);/* V*の境界条件の処理 */
     /* ****************** */
     /* W*(ix,iy,iz)の計算 */
     /* ****************** */
     i=1;
     for (iz=1;iz<=NZ-1;iz++){
       for (ix=1;ix<=NX;ix++){
         for (iy=1;iy<=NY;iy++){
         /* Wに関するDの計算 */
         dw5 = VIS*DY*DZ/DX; dw1 = VIS*DY*DZ/DX;
         dw4 = VIS*DX*DZ/DY; dw2 = VIS*DX*DZ/DY;
         dw7 = VIS*DX*DY/DZ; dw6 = VIS*DX*DY/DZ;
         /* 速度・圧力・温度場の繰り返し計算となっても線形システムの係数となる
            ここの速度は更新されないように(U,V,W)には前の時刻の収束値(UO,VO,WO)を用いる*/
         uw5 = ( UO[ix  ][iy  ][iz] + UO[ix  ][iy  ][iz+1] ) / 2.0;
         uw1 = ( UO[ix-1][iy  ][iz] + UO[ix-1][iy  ][iz+1] ) / 2.0;
         vw4 = ( VO[ix  ][iy  ][iz] + VO[ix  ][iy  ][iz+1] ) / 2.0;
         vw2 = ( VO[ix  ][iy-1][iz] + VO[ix  ][iy-1][iz+1] ) / 2.0;
         ww7 = ( WO[ix  ][iy  ][iz] + WO[ix  ][iy  ][iz+1] ) / 2.0;
         ww6 = ( WO[ix  ][iy  ][iz] + WO[ix  ][iy  ][iz-1] ) / 2.0;
         fw5 = RO*uw5*DY*DZ; fw1 = RO*uw1*DY*DZ;
         fw4 = RO*vw4*DX*DZ; fw2 = RO*vw2*DX*DZ;
         fw7 = RO*ww7*DX*DY; fw6 = RO*ww6*DX*DY;
         /* Wに関するPeの計算 */
         pew5 = fw5/dw5; pew1 = fw1/dw1;
         pew4 = fw4/dw4; pew2 = fw2/dw2;
         pew7 = fw7/dw7; pew6 = fw6/dw6;
         /* 関数A(|Pe|)の計算 */
         if ( NMT==1 ) {/* NMT=1なら 風上法 */
           apew5=1.0; apew1=1.0 ;apew4=1.0; apew2=1.0; apew7=1.0; apew6=1.0;
         }
         if ( NMT==2 ) {/* NMT=2なら べき乗法 */
           apew5 = Max( 0.0, pow( (1.0-0.1*fabs(pew5)),5 ) );
           apew1 = Max( 0.0, pow( (1.0-0.1*fabs(pew1)),5 ) );
           apew4 = Max( 0.0, pow( (1.0-0.1*fabs(pew4)),5 ) );
           apew2 = Max( 0.0, pow( (1.0-0.1*fabs(pew2)),5 ) );
           apew7 = Max( 0.0, pow( (1.0-0.1*fabs(pew7)),5 ) );
           apew6 = Max( 0.0, pow( (1.0-0.1*fabs(pew6)),5 ) );
         }
         /* W* の線形システムの係数の計算 */
         aw5 = ( dw5*apew5 + Max(-fw5,0.0) );
         aw1 = ( dw1*apew1 + Max( fw1,0.0) );
         aw4 = ( dw4*apew4 + Max(-fw4,0.0) );
         aw2 = ( dw2*apew2 + Max( fw2,0.0) );
         aw7 = ( dw7*apew7 + Max(-fw7,0.0) );
         aw6 = ( dw6*apew6 + Max( fw6,0.0) );
         scw = 0.0; spw=0.0;/* 生成項(浮力項)の計算 */
         /* W* の線形システムの B の計算 */
         bp[i] = scw*DX*DY*DZ+( PS[ix][iy][iz]-PS[ix][iy][iz+1] )*DX*DY
               + RO*DX*DY*DZ/DT*WO[ix][iy][iz];
         /* 境界の処理
          1. 境界のW_BND についても変数とみなしてAW3[ix][iy][iz]を求める
          2. その後，境界の項 aw(1,2,4,5,6,7)*W_BND をBに吸収させる
          3. 2に対応して(吸収の後)，境界のaw(1,2,4,5,6,7)をゼロとしてリンクを断つ */
         AW3[ix][iy][iz] = aw1+aw2+aw4+aw5+aw6+aw7 + RO*DX*DY*DZ/DT - spw*DX*DY*DZ;
         if ( ix==NX ) {/* W（右面）*/
           bp[i] = bp[i]+aw5*W[NX+1][iy][iz]; aw5 = 0.0;
         }
         if ( ix==1 ) {/* W（左面）*/
           bp[i] = bp[i]+aw1*W[0][iy][iz]; aw1 = 0.0;
         }
         if ( iy==NY ) {/* W（上面）*/
           bp[i] = bp[i] + aw4*W[ix][NY+1][iz]; aw4 = 0.0;
         }
         if ( iy==1 ) {/* W（下面）*/
           bp[i] = bp[i] + aw2*W[ix][0][iz]; aw2 = 0.0;
         }
         if ( iz==NZ-1 ) {/* W（後面）*/
           bp[i] = bp[i] + aw7*W[ix][iy][NZ]; aw7 = 0.0;
         }
         if ( iz==1 ) {/* W（前面）*/
           bp[i] = bp[i] + aw6*W[ix][iy][0]; aw6 = 0.0;
         }
         /* 線形システム計算用に係数行列を計算
            ALPHAW : W* 計算のための緩和係数 */
         atp3[i] = AW3[ix][iy][iz] / ALPHAW;
         atp5[i] = -aw5; atp1[i] = -aw1;
         atp4[i] = -aw4; atp2[i] = -aw2;
         atp7[i] = -aw7; atp6[i] = -aw6;
         bp[i] = bp[i] + (1.0-ALPHAW)*atp3[i]*WO[ix][iy][iz];
         xx[i]=W[ix][iy][iz];/* 線形システム解法前のxx[i]の初期化 */
         i = i + 1;
         }
       }
     }
     /* --- SIMPLE --- 仮の速度場W*(WS) に関する線形システムの解法 */
     nnx = NX; nny = NY; nnz = NZ-1; nne = NEW; nnxy = nnx*nny;
     /* 1. 直接法 : バンドマトリックスによるガウスの消去法 */
     if (METHOD==1) gb    (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 2. 反復法1 : point-SOR 法 */
     if (METHOD==2) psorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 3. 反復法2 : line-SOR 法 : OMG=1で十分．大きくしすぎると発散する */
     if (METHOD==3) lsorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 4. クリロフ部分空間法1 : 共役残差法 */
     if (METHOD==4) crb   (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 5. クリロフ部分空間法2 : BiCGSTAB */
     if (METHOD==5) bicgb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* METHOD=4,5における探索ベクトル計算時のゼロ除算対策 */
     if (IFLG==2)   psorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 線形システム解法で得られた1次元配列の解を3次元配列に置き換える */
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ-1;iz++){
           k = iy + (ix-1)*nny + (iz-1)*(nnx*nny);
           WS[ix][iy][iz] = xx[k];
         }
       }
     }
     wbnd (WS);/* W*の境界条件の処理 */
}
/*圧力場の計算 Fortran プログラムのコメントを参照してください．*/
void press()
{
     double atp1[NEP1],atp2[NEP1],atp3[NEP1],atp4[NEP1],atp5[NEP1],atp6[NEP1],atp7[NEP1],bp[NEP1],xx[NEP1];
     int ix,iy,iz,k,nnx,nny,nnz,nne,nnxy,i;
     double dmax, du, dv, dw;
     double ap1,  ap2,  ap3, ap4,  ap5, ap6, ap7;
     double apu1, apv2, apv4, apu5, apw6, apw7;
     double dp1,  dp2,  dp4,  dp5, dp6, dp7;
     /* 速度・圧力・温度を連立させて繰り返し計算を行って，連続の式を満たす
        ( U*, V*, W* ) となったかどうかを判定する際の最大値の初期設定 */
     dmax = 0.0;
     i = 1; /* 線形システム解法のため，3次元配列を1次元配列にするための変数 */
     for (iz=1;iz<=NZ;iz++){
       for (ix=1;ix<=NX;ix++){
         for (iy=1;iy<=NY;iy++){
         /* 速度場(U*,V*,W*)が連続の式を満たすかどうかを計算し，満たしていれば，
            この(U*,V*,W*)をもとに，ここで計算するP'(PD)から得られる(U,V,W,P)を
            この時間ステップの収束解とする．なお，(U,V,W)は常に連続の式を満たす
            ことが保証される．したがって，(U*,V*,W*)が連続の式を満たすということは，
            前回の反復において得られた(U,V,W)を用いると，もはやこれを補正する必要が
            なくなるということを意味する．*/
         /* P'(PD) の線形システムのBの計算 */
         bp[i]=RO*( ( US[ix-1][iy  ][iz  ] - US[ix][iy][iz] )*DY*DZ
                  + ( VS[ix  ][iy-1][iz  ] - VS[ix][iy][iz] )*DX*DZ
                  + ( WS[ix  ][iy  ][iz-1] - WS[ix][iy][iz] )*DX*DY );
         /* (U*,V*,W*)による連続の式の収束性をチェックし，div \VEC{v}の
            最大値がDMAXとなるようにする */
         if ( fabs(bp[i])>dmax ){
           dmax = fabs(bp[i]);
         }
         if ( ix==NX ) {/* P'（右面）*/
           ap5=0.0;
         }
         else {
           apu5=AU3[ix][iy][iz]; dp5=(DY*DZ)/apu5; ap5=RO*dp5*DY*DZ;
         }
         if ( ix==1 ) {/* P'（左面）*/
           ap1=0.0;
         }
         else {
           apu1=AU3[ix-1][iy][iz]; dp1=(DY*DZ)/apu1; ap1=RO*dp1*DY*DZ;
         }
         if ( iy==NY ) {/* P'（上面）*/
           ap4=0.0;
         }
         else {
           apv4=AV3[ix][iy][iz]; dp4=(DX*DZ)/apv4; ap4=RO*dp4*DX*DZ;
         }
         if ( iy==1 ) {/* P'（下面）*/
           ap2=0.0;
         }
         else {
           apv2=AV3[ix][iy-1][iz]; dp2=(DX*DZ)/apv2; ap2=RO*dp2*DX*DZ;
         }
         if ( iz==NZ ) {/* P'（後面）*/
           ap7=0.0;
         }
         else {
           apw7=AW3[ix][iy][iz]; dp7=(DX*DY)/apw7; ap7=RO*dp7*DX*DY;
         }
         if ( iz==1 ) {/* P'（前面）*/
           ap6=0.0;
         }
         else {
           apw6=AW3[ix][iy][iz-1]; dp6=(DX*DY)/apw6; ap6=RO*dp6*DX*DY;
         }
         ap3=ap1+ap2+ap4+ap5+ap6+ap7;
         /* 線形システム計算用に係数行列を計算 */
         atp3[i]=ap3;
         atp5[i]=-ap5; atp1[i]=-ap1; atp4[i]=-ap4; atp2[i]=-ap2; atp7[i]=-ap7; atp6[i]=-ap6;
         /*１次独立な解を得るための処理 : IRELP=1 : 直接法においては必須
             (IX=1,IY=1,IZ=1 ---> K=1を基準点とし，常にP[1][1][1]=PD[1][1][1]=0とする) */
         if (IRELP==1 && ix==1 && iy==1 && iz==1) {
           atp4[1]=0.0; atp5[1]=0.0; atp7[1]=0.0; bp[1]=0.0;
           atp2[2]=0.0; /* K=2 の点の処理(K=1とのリンクを断つ) */
           atp1[1+NY]=0.0; /* K=1+NY の点の処理(K=1とのリンクを断つ) */
           atp6[1+NY+NX*NY]=0.0;/* K=1+NY+NX*NY の点の処理(K=1とのリンクを断つ)*/
         }
         xx[i]=PD[ix][iy][iz];/* 線形システム解法前のxx[i]の初期化 */
         i=i+1;
         }
       }
     }
     /* --- SIMPLE --- 圧力補正 P'(PD) に関する解法 */
     nnx = NX; nny = NY; nnz = NZ; nne = NEP; nnxy = nnx*nny;
     /* 1. 直接法 : バンドマトリックスによるガウスの消去法 */
     if (METHOD==1) gb    (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 2. 反復法1 : point-SOR 法 */
     if (METHOD==2) psorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 3. 反復法2 : line-SOR 法 : OMG=1で十分．大きくしすぎると発散する */
     if (METHOD==3) lsorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 4. クリロフ部分空間法1 : 共役残差法 */
     if (METHOD==4) crb   (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 5. クリロフ部分空間法2 : BiCGSTAB */
     if (METHOD==5) bicgb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* METHOD=4,5における探索ベクトル計算時のゼロ除算対策 */
     if (IFLG==2)   psorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
         k=iy+(ix-1)*NY+(iz-1)*(nnx*nny);
         /* 圧力の相対性の処理 : 基準値の設定
            1次独立な解を求める場合は圧力の基準点が強制的に設定される */
         if (IRELP != 2)  PD[ix][iy][iz] = xx[k]; /* 圧力の基準点を設けない場合 */
         /* 1次従属な解のうちの1つを求めた後，圧力基準を設ける場合 (IRELP=2) */
         if (IRELP == 2) PD[ix][iy][iz] = xx[k]-xx[1]; /* P'(1,1)=0 ---> P(1,1)=0 */
         }
       }
     }
     pdbnd (PD);/* 圧力補正の境界条件の処理 */
     for (ix=1;ix<=NX-1;ix++){/* Uの修正 */
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
            du = DY*DZ/AU3[ix][iy][iz];
             U[ix][iy][iz]=US[ix][iy][iz]+du*(PD[ix][iy][iz]-PD[ix+1][iy][iz]);
         }
       }
     }
     for (ix=1;ix<=NX;ix++){/* Vの修正 */
       for (iy=1;iy<=NY-1;iy++){
         for (iz=1;iz<=NZ;iz++){
           dv = DX*DZ/AV3[ix][iy][iz];
           V[ix][iy][iz]=VS[ix][iy][iz]+dv*(PD[ix][iy][iz]-PD[ix][iy+1][iz]);
         }
       }
     }
     for (ix=1;ix<=NX;ix++){/* Wの修正 */
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ-1;iz++){
           dw = DX*DY/AW3[ix][iy][iz];
           W[ix][iy][iz]=WS[ix][iy][iz]+dw*(PD[ix][iy][iz]-PD[ix][iy][iz+1]);
         }
       }
     }
     /* 新たに得られた速度を用いて境界条件を処理する */
     ubnd (U); vbnd (V); wbnd(W);
     for (ix=1;ix<=NX;ix++){/* P* の修正: P = P* + P': ただしPの値はPSに代入 */
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
           PS[ix][iy][iz] = PS[ix][iy][iz] + ALPHAP*PD[ix][iy][iz];
         }
       }
     }
     psbnd (PS);/* 新たに得られたP*を用いて境界条件を処理する */
     if ( dmax > EPSC ) IFLGC=1; /* 速度場 (U*,V*,W*) が EPSC 以下で連続の式を満たすかどうかを判定 */
}
/*温度場の計算*/
void caltem()
{
     double atp1[NEP1],atp2[NEP1],atp3[NEP1],atp4[NEP1],atp5[NEP1],atp6[NEP1],atp7[NEP1],bp[NEP1],xx[NEP1];
     int ix,iy,iz,i,k;
     int nnx,nny,nnz,nne,nnxy;
     double dt1,dt2,dt4,dt5,dt6,dt7;
     double ut1,vt2,vt4,ut5,wt6,wt7;
     double ft1,ft2,ft4,ft5,ft6,ft7;
     double pet1,pet2,pet4,pet5,pet6,pet7;
     double apet1,apet2,apet4,apet5,apet6,apet7;
     double at1,at2,at3,at4,at5,at6,at7;
     double sct,spt;
     i=1;
     for (iz=1;iz<=NZ;iz++){
       for (ix=1;ix<=NX;ix++){
         for (iy=1;iy<=NY;iy++){
         /* Tに関するDの計算 */
         dt5=ALP*DY*DZ/DX; dt1=ALP*DY*DZ/DX;
         dt4=ALP*DX*DZ/DY; dt2=ALP*DX*DZ/DY;
         dt7=ALP*DX*DY/DZ; dt6=ALP*DX*DY/DZ;
         /* 速度・圧力・温度場の繰り返し計算の際のU,V,Wを線形システムの係数に反映させる */
         ut5=U[ix][iy][iz]; ut1=U[ix-1][iy][iz]; vt4=V[ix][iy][iz]; vt2=V[ix][iy-1][iz];
         wt7=W[ix][iy][iz]; wt6=W[ix][iy][iz-1];
         ft5=RO*ut5*DY*DZ; ft1=RO*ut1*DY*DZ; ft4=RO*vt4*DX*DZ; ft2=RO*vt2*DX*DZ; ft7=RO*wt7*DX*DY; ft6=RO*wt6*DX*DY;
         /* Tに関するPeの計算 */
         pet5=ft5/dt5; pet1=ft1/dt1; pet4=ft4/dt4; pet2=ft2/dt2; pet7=ft7/dt7; pet6=ft6/dt6;
         /* 関数A(|Pe|)の計算 */
         if ( NMT==1 ) {/* NMT=1なら 風上法 */
           apet5=1.0; apet1=1.0; apet4=1.0; apet2=1.0; apet7=1.0; apet6=1.0;
         }
         if ( NMT==2 ) {/* NMT=2なら べき乗法 */
           apet5 = Max( 0.0, pow( (1.0-0.1*fabs(pet5)),5 ) );
           apet1 = Max( 0.0, pow( (1.0-0.1*fabs(pet1)),5 ) );
           apet4 = Max( 0.0, pow( (1.0-0.1*fabs(pet4)),5 ) );
           apet2 = Max( 0.0, pow( (1.0-0.1*fabs(pet2)),5 ) );
           apet7 = Max( 0.0, pow( (1.0-0.1*fabs(pet7)),5 ) );
           apet6 = Max( 0.0, pow( (1.0-0.1*fabs(pet6)),5 ) );
         }
         /* T* の線形システムの係数の計算 */
         at5 = ( dt5*apet5 + Max(-ft5,0.0) );
         at1 = ( dt1*apet1 + Max( ft1,0.0) );
         at4 = ( dt4*apet4 + Max(-ft4,0.0) );
         at2 = ( dt2*apet2 + Max( ft2,0.0) );
         at7 = ( dt7*apet7 + Max(-ft7,0.0) );
         at6 = ( dt6*apet6 + Max( ft6,0.0) );
         sct = 0.0; spt=0.0; /* 生成項 */
         /* T* の線形システムのBの計算 */
         bp[i] = sct*DX*DY*DZ + RO*DX*DY*DZ/DT*TO[ix][iy][iz];
         /* 境界の処理
          1. 境界のT_BND についても変数とみなしてAT3を求める
          2. その後，境界の項 at(1,2,4,5,6,7)*T_BND をBに吸収させる
          3. 2に対応して(吸収の後)，境界のat(1,2,4,5,6,7)をゼロとしてリンクを断つ */
         at3 = at1+at2+at4+at5+at6+at7 + RO*DX*DY*DZ/DT - spt*DX*DY*DZ;
         if ( ix==NX ) {/* T（右面）*/
           bp[i] = bp[i]+at5*T[NX+1][iy][iz]; at5=0.0;
         }
         if ( ix==1 ) {/* T（左面）*/
           bp[i] = bp[i]+at1*T[0][iy][iz]; at1=0.0;
         }
         if ( iy==NY ) {/* T（上面）*/
           bp[i] = bp[i] + at4*T[ix][NY+1][iz]; at4=0.0;
         }
         if ( iy==1 ) {/* T（下面）*/
           bp[i] = bp[i] + at2*T[ix][0][iz]; at2=0.0;
         }
         if ( iz==NZ ) {/* T（後面）*/
           bp[i] = bp[i] + at7*T[ix][iy][NZ+1]; at7=0.0;
         }
         if ( iz==1 ) {/* T（前面）*/
           bp[i] = bp[i] + at6*T[ix][iy][0]; at6=0.0;
         }
         /* 線形システム計算用に係数行列を計算
            ALPHAT : T 計算のための緩和係数 */
         atp3[i] = at3 / ALPHAT;
         atp5[i] = -at5; atp1[i] = -at1;
         atp4[i] = -at4; atp2[i] = -at2;
         atp7[i] = -at7; atp6[i] = -at6;
         bp[i] = bp[i] + (1.0-ALPHAT)*atp3[i]*TO[ix][iy][iz];
         xx[i]=T[ix][iy][iz];/* 線形システム解法前のxx[i]の初期化 */
         i = i + 1;
         }
       }
     }
     /* --- SIMPLE --- 仮の速度場U*(US) に関する線形システムの解法 */
     nnx=NX; nny=NY; nnz=NZ; nne=NEP, nnxy=nnx*nny;
     /* 1. 直接法 : バンドマトリックスによるガウスの消去法 */
     if (METHOD==1) gb    (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 2. 反復法1 : point-SOR 法 */
     if (METHOD==2) psorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 3. 反復法2 : line-SOR 法 : OMG=1で十分．大きくしすぎると発散する */
     if (METHOD==3) lsorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 4. クリロフ部分空間法1 : 共役残差法 */
     if (METHOD==4) crb   (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 5. クリロフ部分空間法2 : BiCGSTAB */
     if (METHOD==5) bicgb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* METHOD=4,5における探索ベクトル計算時のゼロ除算対策 */
     if (IFLG==2)   psorb (atp1,atp2,atp3,atp4,atp5,atp6,atp7,bp,xx,nnx,nny,nnz,nne,nnxy);
     /* 線形システム解法で得られた1次元配列の解を3次元配列に置き換える */
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
           k=iy+(ix-1)*nny+(iz-1)*(nnx*nny);
           T[ix][iy][iz] = xx[k];
         }
       }
     }
     tbnd (T);/* Tの境界条件の処理 */
}
/* 速度 ( U* と U ) の境界条件の処理 */
void ubnd(double XU[NX1][NY2][NZ2])
{
     int ix,iy,iz;
     for (iy=1;iy<=NY;iy++){/*U（左面）*/
       for (iz=1;iz<=NZ;iz++){
         XU[0][iy][iz] = 0.0;
       }
     }
     for (iy=1;iy<=NY;iy++){/*U（右面）*/
       for (iz=1;iz<=NZ;iz++){
         XU[NX][iy][iz] = 0.0;
       }
     }
     for (ix=0;ix<=NX;ix++){/*U（下面）*/
       for (iz=1;iz<=NZ;iz++){
         XU[ix][0][iz] = -XU[ix][1][iz];
       }
     }
     for (ix=0;ix<=NX;ix++){/*U（上面）*/
       for (iz=1;iz<=NZ;iz++){
         XU[ix][NY+1][iz] = -XU[ix][NY][iz];
       }
     }
     for (ix=0;ix<=NX;ix++){/*U（前面）*/
       for (iy=0;iy<=NY+1;iy++){
         XU[ix][iy][0] = -XU[ix][iy][1];
       }
     }
     for (ix=0;ix<=NX;ix++){/*U（後面）*/
       for (iy=0;iy<=NY+1;iy++){
         XU[ix][iy][NZ+1] = -XU[ix][iy][NZ];
       }
     }
}
/* 速度 ( V* と V ) の境界条件の処理 */
void vbnd(double XV[NX2][NY1][NZ2])
{
     int ix,iy,iz;
     for (iy=1;iy<=NY-1;iy++){/*V（左面）*/
       for (iz=1;iz<=NZ;iz++){
         XV[0][iy][iz] = -XV[1][iy][iz];
       }
     }
     for (iy=1;iy<=NY-1;iy++){/*V（右面）*/
       for (iz=1;iz<=NZ;iz++){
         XV[NX+1][iy][iz] = -XV[NX][iy][iz];
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*V（上面）*/
       for (iz=1;iz<=NZ;iz++){
         XV[ix][NY][iz] = 0.0;
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*V（下面）*/
       for (iz=1;iz<=NZ;iz++){
         XV[ix][0][iz] = 0.0;
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*V（前面）*/
       for (iy=0;iy<=NY;iy++){
         XV[ix][iy][0] = -XV[ix][iy][1];
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*V（後面）*/
       for (iy=0;iy<=NY;iy++){
         XV[ix][iy][NZ+1] = -XV[ix][iy][NZ];
       }
     }
}
/* 速度 ( W* と W ) の境界条件の処理 */
void wbnd(double XW[NX2][NY2][NZ1])
{
     int ix,iy,iz;
     for (iy=1;iy<=NY;iy++){/*W（左面）*/
       for (iz=1;iz<=NZ-1;iz++){
         XW[0][iy][iz] = -XW[1][iy][iz];
       }
     }
     for (iy=1;iy<=NY;iy++){/*W（右面）*/
       for (iz=1;iz<=NZ-1;iz++){
         XW[NX+1][iy][iz] = -XW[NX][iy][iz];
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*W（下面）*/
       for (iz=1;iz<=NZ-1;iz++){
         XW[ix][0][iz] = -XW[ix][1][iz];
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*W（上面）*/
       for (iz=1;iz<=NZ-1;iz++){
         XW[ix][NY+1][iz] = -XW[ix][NY][iz];
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*W（前面）*/
       for (iy=0;iy<=NY+1;iy++){
         XW[ix][iy][0] = 0.0;
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*W（後面）*/
       for (iy=0;iy<=NY+1;iy++){
         XW[ix][iy][NZ] = 0.0;
       }
     }
}
/*温度の境界条件の処理*/
void tbnd(double XT[NX2][NY2][NZ2])
{
     int ix,iy,iz;
      for (iy=1;iy<=NY;iy++){/* 左面 */
        for (iz=1;iz<=NZ;iz++){
          XT[0][iy][iz] = 2.0 * ( +0.5 ) - XT[1][iy][iz];
        }
      }
      for (iy=1;iy<=NY;iy++){/* 右面 */
        for (iz=1;iz<=NZ;iz++){
          XT[NX+1][iy][iz] = 2.0 * ( -0.5 ) - XT[NX][iy][iz];
        }
      }
      for (ix=0;ix<=NX+1;ix++){/* 下面 */
        for (iz=1;iz<=NZ;iz++){
          XT[ix][0][iz] = XT[ix][1][iz];
        }
      }
      for (ix=0;ix<=NX+1;ix++){/* 上面 */
        for (iz=1;iz<=NZ;iz++){
          XT[ix][NY+1][iz] = XT[ix][NY][iz];
        }
      }
      for (ix=0;ix<=NX+1;ix++){/* 前面 */
        for (iy=0;iy<=NY+1;iy++){
          XT[ix][iy][0] = XT[ix][iy][1];
        }
      }
      for (ix=0;ix<=NX+1;ix++){/* 後面 */
        for (iy=0;iy<=NY+1;iy++){
          XT[ix][iy][NZ+1] = XT[ix][iy][NZ];
        }
      }
}
/*圧力補正の境界条件の処理*/
void pdbnd(double XP[NX2][NY2][NZ2])
{
     int ix,iy,iz;
      for (iy=1;iy<=NY;iy++){/* 左面 */
        for (iz=1;iz<=NZ;iz++){
          XP[0][iy][iz] = XP[1][iy][iz];
        }
      }
      for (iy=1;iy<=NY;iy++){/* 右面 */
        for (iz=1;iz<=NZ;iz++){
          XP[NX+1][iy][iz] = XP[NX][iy][iz];
        }
      }
      for (ix=0;ix<=NX+1;ix++){/* 下面 */
        for (iz=1;iz<=NZ;iz++){
          XP[ix][0][iz] = XP[ix][1][iz];
        }
      }
      for (ix=0;ix<=NX+1;ix++){/* 上面 */
        for (iz=1;iz<=NZ;iz++){
          XP[ix][NY+1][iz] = XP[ix][NY][iz];
        }
      }
      for (ix=0;ix<=NX+1;ix++){/* 前面 */
        for (iy=0;iy<=NY+1;iy++){
          XP[ix][iy][0] = XP[ix][iy][1];
        }
      }
      for (ix=0;ix<=NX+1;ix++){/* 後面 */
        for (iy=0;iy<=NY+1;iy++){
          XP[ix][iy][NZ+1] = XP[ix][iy][NZ];
        }
      }
}
/*圧力の境界条件の処理*/
void psbnd(double XP[NX2][NY2][NZ2])
{
     int ix,iy,iz;
      for (iy=1;iy<=NY;iy++){/* 左面 */
        for (iz=1;iz<=NZ;iz++){
          XP[0][iy][iz] = XP[1][iy][iz];
        }
      }
      for (iy=1;iy<=NY;iy++){/* 右面 */
        for (iz=1;iz<=NZ;iz++){
          XP[NX+1][iy][iz] = XP[NX][iy][iz];
        }
      }
      for (ix=0;ix<=NX+1;ix++){/* 下面 */
        for (iz=1;iz<=NZ;iz++){
          XP[ix][0][iz] = XP[ix][1][iz];
        }
      }
      for (ix=0;ix<=NX+1;ix++){/* 上面 */
        for (iz=1;iz<=NZ;iz++){
          XP[ix][NY+1][iz] = XP[ix][NY][iz];
        }
      }
      for (ix=0;ix<=NX+1;ix++){/* 前面 */
        for (iy=0;iy<=NY+1;iy++){
          XP[ix][iy][0] = XP[ix][iy][1];
        }
      }
      for (ix=0;ix<=NX+1;ix++){/* 後面 */
        for (iy=0;iy<=NY+1;iy++){
          XP[ix][iy][NZ+1] = XP[ix][iy][NZ];
        }
      }
}
void prout()/*データ出力用*/
{
    fwrite(U,  sizeof(double), NX1*NY2*NZ2, out_11);
    fwrite(V,  sizeof(double), NX2*NY1*NZ2, out_12);
    fwrite(W,  sizeof(double), NX2*NY2*NZ1, out_13);
    fwrite(PS, sizeof(double), NX2*NY2*NZ2, out_14);
    fwrite(PD, sizeof(double), NX2*NY2*NZ2, out_14);
    fwrite(T,  sizeof(double), NX2*NY2*NZ2, out_15);
}
/* Tecplot用データの出力 */
void tecplt ( char *fname_tec )
{
     int ix,iy,iz;
     double x,y,z,uu,vv,ww,tt;

     out_21=fopen(fname_tec,"wt");
     fprintf (out_21," VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"T\"\n");
     fprintf (out_21," ZONE I=%3d, J=%3d, K=%3d, F=POINT \n",NX1,NY1,NZ1);

     for (iz=0;iz<=NZ;iz++){
       for (iy=0;iy<=NY;iy++){
         for (ix=0;ix<=NX;ix++){
           x = DX * (double)ix;
           y = DY * (double)iy;
           z = DZ * (double)iz;
           uu = (U[ix][iy][iz]+U[ix][iy+1][iz]+U[ix][iy][iz+1]+U[ix][iy+1][iz+1])/4.0;
           vv = (V[ix][iy][iz]+V[ix+1][iy][iz]+V[ix][iy][iz+1]+V[ix+1][iy][iz+1])/4.0;
           ww = (W[ix][iy][iz]+W[ix+1][iy][iz]+W[ix][iy+1][iz]+W[ix+1][iy+1][iz])/4.0;
           tt = (T[ix][iy][iz]+  T[ix+1][iy][iz]+  T[ix][iy+1][iz]+  T[ix+1][iy+1][iz]
                +T[ix][iy][iz+1]+T[ix+1][iy][iz+1]+T[ix][iy+1][iz+1]+T[ix+1][iy+1][iz+1] )/8.0;
         fprintf(out_21,"%12.3lE %12.3lE %12.3lE %12.3lE %12.3lE %12.3lE %12.3lE\n",x,y,z,uu,vv,ww,tt);
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
void gb (double a1[NEP1], double a2[NEP1], double a3[NEP1],
         double a4[NEP1], double a5[NEP1], double a6[NEP1], double a7[NEP1],
         double b[NEP1],  double x[NEP1],
         int nx, int ny, int nz, int ne, int nxy)
{
     double a[NXYNXY][NEP1];
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
        IFLG = 1; return;/* 線形システムの計算を発散に設定 */
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
     ITR = 1;
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
void psorb (double a1[NEP1], double a2[NEP1], double a3[NEP1], double a4[NEP1], 
            double a5[NEP1], double a6[NEP1], double a7[NEP1],
            double b[NEP1],  double x[NEP1],
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
void lsorb (double at1[NEP1] ,double at2[NEP1], double at3[NEP1],
            double at4[NEP1], double at5[NEP1], double at6[NEP1], double at7[NEP1],
            double bx[NEP1],  double xn[NEP1],
            int nx, int ny, int nz, int ne, int nxy)
{
     double x1[NEP1],xo[NEP1],xold[NEP1],a[NEP1],b[NEP1],c[NEP1],d[NEP1],u[NEP1],y[NEP1];
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
void crb (double a1[NEP1], double a2[NEP1], double a3[NEP1],
          double a4[NEP1], double a5[NEP1], double a6[NEP1], double a7[NEP1],
          double b[NEP1],  double xp[NEP1],
          int nx, int ny, int nz, int ne, int nxy)
{
     double r[NEP1], p[NEP1], ap[NEP1], ar[NEP1], x[NEP1], xold[NEP1];/* 作業用配列 */
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
double provv (double a[NEP1], double b[NEP1], int ne)
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
void promv (double a1[NEP1], double a2[NEP1], double a3[NEP1], double a4[NEP1],
            double a5[NEP1], double a6[NEP1], double a7[NEP1],
            double b[NEP1],  double c[NEP1],
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
void bicgb ( double a1[NEP1], double a2[NEP1], double a3[NEP1], double a4[NEP1], double a5[NEP1],
             double a6[NEP1], double a7[NEP1], double b[NEP1],  double x[NEP1],
             int nx, int ny, int nz, int ne, int nxy)
{
     double r[NEP1], ap[NEP1], at[NEP1], p[NEP1], s[NEP1], t[NEP1], xold[NEP1];
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
