*********************************************************************
* ファイル名  ：MAT3D.FOR                                           *
* タイトル    ：線形システム解法プログラム                          *
* 製作者      ：平野　博之                                          *
* 所属        ：岡山理科大学 工学部 応用化学科                      *
* 製作日      ：2003.12.25                                          *
* 言語        ：FORTRAN                                             *
*********************************************************************
*［内容］                                                           *
*    本プログラムは本書の15.6節にある2次元の例題で作成された線形    *
*    システム解法の1例である．                                      *
*    解法は，直接法，反復法，そしてクリロフ部分空間法を用いている． *
*    解法アルゴリズムとそのサブルーチン名は以下に示す通りである．   *
*    -------------------------------------------------------------  *
*    直接法                              サブルーチン名             *
*    -------------------------------------------------------------  *
*    LU分解法(軸選択付)               ---> DLUP                     *
*    Gaussの消去法(軸選択付)          ---> DFGP                     *
*    Gaussの消去法(軸選択無)          ---> DFG                      *
*    Gauss-Jordan法(軸選択付)         ---> DGJP                     *
*    Gaussの消去法(2次元(5点)差分     ---> DGBND                    *
*    バンドマトリックス用: 軸選択無)                                *
*    -------------------------------------------------------------  *
*    反復法                              サブルーチン名             *
*    -------------------------------------------------------------  *
*    JACOBI法                         ---> HJACOB                   *
*    point-SOR法                      ---> HSOR                     *
*    point-SOR法(バンドマトリックス用)---> SORBND                   *
*    line-SOR法(バンドマトリックス用) ---> LBLBND                   *
*    注意: SOR法で，ω=1とするとGauss-Seidel法となる．              *
*    -------------------------------------------------------------  *
*    クリロフ部分空間法                  サブルーチン名             *
*    -------------------------------------------------------------  *
*    共役残差法                       ---> KCR                      *
*    共役残差法(バンドマトリックス用) ---> KCRBND                   *
*    Bi-CGSTAB                        ---> KBICG                    *
*    Bi-CGSTAB(バンドマトリックス用)  ---> KBIBND                   *
*［  例  題  ］                                                     *
*    3次元偏微分方程式 ∇^{2}f + ∇_{x}f + ∇_{y}f + ∇_{z}f= b     *
*  について，2次精度中心差分近似を用いて離散化せよ．そして以下に示す*
*  計算領域を考え，通し番号をつけられた計算格子の番号そのものが解と *
*  なるような線形システム(AX=B)を作り，各種の解法を用いて解け．     *
*    以下に示すような3次元の計算領域について，これをx,y,z方向にそれ *
*  ぞれ5等分ずつしこれに通し番号を(1次元的に)m=1から125まで付ける． *
*  そして，X_{i,j,k}=X_{m}=mが解となるようにマトリックス(A,B)を設定 *
*  する．詳細は2次元のものと同様．                                  *
*                         ------------------------                  *
*                       / 105| 110| 115| 120| 125/ |                *
*                      /-------------------------  |                *
*                     /                        /   |                *
*                    /                        /    |                *
*                   +------------------------+     |                *
*            (NY)5  |  5 | 10 | 15 | 20 | 25 |     |                *
*                   +------------------------+     |                *
*                4  |  4 |  9 | 14 | 19 | 24 |     |                *
*                   +------------------------+     |                *
*                3  |  3 |  8 | 13 | 18 | 23 |     /                *
*                   +------------------------+    /                 *
*                2  |  2 |  7 | 12 | 17 | 22 |   /k=1,5(NZ)         *
*                   +------------------------+  /                   *
*                1  |  1 |  6 | 11 | 16 | 21 | /                    *
*                   +------------------------+/                     *
*                j   i=1     2    3    4    5(NX)                   *
*   f の定義点は各格子の中心であるとする．f の定義点間距離，および  *
*   各計算格子の幅はいずれもΔ=1とする．また，境界条件に関しては，  *
*   計算領域外の定義点の値をゼロとして扱う．                        *
*   通し番号:m=(k-1)*(NX*NY) + (i-1)*NY + j                         *
* [　解　答　例　]                                                  *
* ( 1 )　離　散　化                                                 *
*   ∇^{2}f + ∇_{x}f + ∇_{y}f ∇_{z}f = b                         *
*   ===> ( f_{i-1,j,k}-2f_{i,j,k}+f_{i+1,j,k} ) / Δ^{2}            *
*       +( f_{i,j-1,k}-2f_{i,j,k}+f_{i,j+1,k} ) / Δ^{2}            *
*       +( f_{i,j,k-1}-2f_{i,j,k}+f_{i,j,k+1} ) / Δ^{2}            *
*       +( f_{i+1,j,k}-f_{i-1,j,k} ) / (2Δ)                        *
*       +( f_{i,j+1,k}-f_{i,j-1,k} ) / (2Δ)                        *
*       +( f_{i,j,k+1}-f_{i,j,k-1} ) / (2Δ)                        *
*       = b_{m}   where m=(k-1)*(NX*NY) + (i-1)*NY + j              *
*   ===> (  1/Δ^{2} - 1/(2Δ)  ) f_{i-1,j,k}                       *
*        ~~~~~~~~~~~~~~~~~~~~~~~~->AT1                              *
*       +( -2/Δ^{2} - 2/Δ^{2} - 2/Δ^{2} ) f_{i,j,k}              *
*        ~~~~~~~~~~~~~~~~~~~~~~~~->AT3                              *
*       +(  1/Δ^{2} + 1/(2Δ)  ) f_{i+1,j,k}                       *
*        ~~~~~~~~~~~~~~~~~~~~~~~~->AT5                              *
*       +(  1/Δ^{2} - 1/(2Δ) )  f_{i,j-1,k}                       *
*        ~~~~~~~~~~~~~~~~~~~~~~~~->AT2                              *
*       +(  1/Δ^{2} + 1/(2Δ) )  f_{i,j+1,k}                       *
*        ~~~~~~~~~~~~~~~~~~~~~~~~->AT4                              *
*       +(  1/Δ^{2} - 1/(2Δ) )  f_{i,j,k-1}                       *
*        ~~~~~~~~~~~~~~~~~~~~~~~~->AT6                              *
*       +(  1/Δ^{2} + 1/(2Δ) )  f_{i,j,k+1}                       *
*        ~~~~~~~~~~~~~~~~~~~~~~~~->AT7                              *
*       = b_{m}   where m=(k-1)*(NX*NY) + (i-1)*NY + j              *
*       注意：Δ=1とする                                            *
*   ===> Af = B  ===> AX = B                                        *
* ( 2 )　マトリックスの設定                                         *
*  2-1 Aの設定                                                      *
*   Aは78行から91行にあるように設定する．A1,A2,A3,A4,A5,A6,A7 は    *
*   バンドマトリックス用の解法のためのものである．これ以外の，ゼロの*
*   要素を含む全(フル)マトリックスを扱う場合は，AFを用いる．        *
*  2-2 Bの設定                                                      *
*   78行から91行において，解となるfの値を代入したものをBとすればよ  *
*   い．具体的には，以下のようにして順にマトリックスを求めてゆく．  *
*   m= 1 (i=1,j=1,k=1)                                              *
*   ---> f_{i-1,j  ,k  }:領域外となるので0とする                    *
*   ---> f_{i  ,j  ,k  }:(k-1)*(NX*NY)+(i-1)*NY+jを代入             *
*   ---> f_{i+1,j  ,k  }:(k-1)*(NX*NY)+((i+1)-1)*NY+jを代入         *
*   ---> f_{i  ,j-1,k  }:領域外となるので0とする                    *
*   ---> f_{i  ,j+1,k  }:(k-1)*(NX*NY)+(i-1)*NY+(j+1)を代入         *
*   ---> f_{i  ,j  ,k-1}:領域外となるので0とする                    *
*   ---> f_{i  ,j  ,k+1}:領域外となるので0とする                    *
*   b_{m}:以上のfの値とΔ(=1)から求まる                             *
*   m=2 (i=1,j=2,k=1), m=3 (i=1,j=3,k=1), m=4 (i=1,j=4,k=1),        *
*   ................., m=NX*NY=25 (i=NX=5,j=NY=5,k=1),              *
*   m=NX*NY+1=26 (i=1,j=1,k=2), ..............,                     *
*   m=NX*NY*NZ=125 (i=NX5,j=NY=5,k=NZ=5)                            *
*   係数行列Aを配列AT1,AT2,AT3,AT4,AT5,AT6,AT7(バンドマトリックス用)*
*   およびA(ゼロを含むフルマトリックス用)に格納する．               *
*   i=3, j=3, k=3 (m=59)における，AT1,AT2,AT3,AT4,AT5,AT6,AT7とA の *
*   関係                     k=1                                    *
*        +------------------/-----+                                 *
* (NY)5  |    |    |    |  / |    |                                 *
*        +---------------AT6------+                                 *
*     4  |    |    |AT4 /    |    |    m=(k-1)*(NX*NY)+(i-1)*NY+j   *
*        +-------------/----------+                                 *
*     3  |    | AT1|AT3 |AT5 |    |    AT1(m) -> A(m,m-NY)          *
*        +---------/--------------+    AT2(m) -> A(m,m- 1)          *
*     2  |    |   /|AT2 |    |    |    AT3(m) -> A(m,m   )          *
*        +-----AT7----------------+    AT4(m) -> A(m,m+ 1)          *
*     1  |      /  |    |    |    |    AT5(m) -> A(m,m+NY)          *
*        +-----/------------------+    AT6(m) -> A(m,m-(NX*NY))     *
*     j   i=1 /   2    3    4    5(NX) AT7(m) -> A(m,m+(NX*NY))     *
*            /                         B(m)   -> B(m)               *
*         k=5(NZ)                                                   *
*********************************************************************
      PROGRAM MAT3D
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
* パラメータ変数(NX0と0が付いているのはサブルーチン引数にはPARAMETER
* 変数であるNX0を直接渡せないため．後でNX=NX0と代入し直してこれを渡し
* ている．) NX:x方向格子数; NY:y方向格子数; NZ:z方向格子数; NXY:NX*NY;
* NE:全格子数=NX*NY*NZ
      PARAMETER ( NX0=5, NY0=5, NZ0=5, NXY0=25, NE0=125 )
      DIMENSION A(NE0,NE0),B(NE0),X1(NE0),X3(NX0,NY0,NZ0)
*バンドマトリックスの配列の定義
      DIMENSION AT1(NE0),AT2(NE0),AT3(NE0),AT4(NE0),AT5(NE0),
     &          AT6(NE0),AT7(NE0)
* 作業用配列
      DIMENSION W1(NE0),W2(NE0),W3(NE0),W4(NE0),W5(NE0),W6(NE0),
     $          W7(NE0),W8(NE0),AGB(-NXY0:NXY0,NE0)
      DIMENSION IPV(NE0)
*上記離散化方程式のΔ=1を設定
      DELTA = 1.0D0
*PARAMETER変数の値をSUBROUTINE引数に渡すために代入
      NX  = NX0
      NY  = NY0
      NZ  = NZ0
      NXY = NXY0
      NE  = NE0
* 戻り点 : 解法の変更毎の戻り点
  700 CONTINUE
* A,B を設定するためにあらかじめX3の値を格子番号と同じになるようにする
      I = 1
      DO 10 IZ = 1,NZ
        DO 20 IX = 1,NX
          DO 30 IY = 1,NY
            X3(IX,IY,IZ) = DBLE(I)
            I =I + 1
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
* X3 以外の配列のゼロクリア
      DO 40 I = 1,NE
        AT1(I) = 0.0D0
        AT2(I) = 0.0D0
        AT3(I) = 0.0D0
        AT4(I) = 0.0D0
        AT5(I) = 0.0D0
        AT6(I) = 0.0D0
        AT7(I) = 0.0D0
        B(I) = 0.0D0
        X1(I) = 0.0D0
        W1(I)=0.0D0
        W2(I)=0.0D0
        W3(I)=0.0D0
        W4(I)=0.0D0
        W5(I)=0.0D0
        W6(I)=0.0D0
        W7(I)=0.0D0
        W8(I)=0.0D0
        DO 50 II = 1,NE
          A(I,II) = 0.0D0
   50   CONTINUE
        DO 55 III = -NXY,NXY
          AGB(III,I) = 0.0D0
   55   CONTINUE
   40 CONTINUE
* X_{i,j,k}=m=(k-1)*(NX*NY)+(i-1)*NY+j となるように A,B を求める
      I = 1
      DO 60 IZ = 1,NZ
        DO 70 IX = 1,NX
          DO 80 IY = 1,NY
*           f(i,j)は常に計算領域内
            AT3(I) = -2.0D0/DELTA**2 - 2.0D0/DELTA**2 -2.0D0/DELTA**2
            A(I,I) = AT3(I)
*           AT3*f(i,j)をB3とする
            B3 = X3(IX,IY,IZ)
     $         * ( -2.0D0/DELTA**2 - 2.0D0/DELTA**2 - 2.0D0/DELTA**2 )
*           f(i-1,j,k)が計算領域外（左側）の境界条件の処理
            IF (IX.EQ.1) THEN
*             f(0,j,k)=0として処理する
              AT1(I) = 0.0D0
*             AT1*f(0,j,k)をB1とする
              B1 = 0.0D0
*           f(i-1,j,k)が計算領域内の場合
            ELSE
              AT1(I) = 1.0D0/DELTA**2 - 1.0D0/(2.0D0*DELTA)
              A(I,I-NY) = AT1(I)
              B1=X3(IX-1,IY,IZ)*( 1.0D0/DELTA**2-1.0D0/(2.0D0*DELTA) )
            END IF
*           f(i+1,j,k)が計算領域外（右側）の境界条件の処理
            IF (IX.EQ.NX) THEN
*             f(NX,j,k)=0として処理
              AT5(I) = 0.0D0
*             AT5*f(NX,j,k)=B5とする
              B5 = 0.0D0
*           f(i+1,j,k)が計算領域内の場合
            ELSE
              AT5(I) = 1.0D0/DELTA**2 + 1.0D0/(2.0D0*DELTA)
              A(I,I+NY) = AT5(I)
              B5 = X3(IX+1,IY,IZ)*( 1.0D0/DELTA**2+1.0D0/(2.0D0*DELTA) )
            END IF
*           f(i,j-1,k)が計算領域外（下側）の境界条件の処理
            IF (IY.EQ.1) THEN
*             f(i,0,k)=0として処理
              AT2(I) = 0.0D0
*             AT2*f(i,0,k)=B2とする
              B2 = 0.0D0
*           f(i,j-1,k)が計算領域内の場合
            ELSE
              AT2(I) = 1.0D0/DELTA**2 - 1.0D0/(2.0D0*DELTA)
              A(I,I-1) = AT2(I)
              B2=X3(IX,IY-1,IZ)*( 1.0D0/DELTA**2-1.0D0/(2.0D0*DELTA) )
            END IF
*           f(i,j+1,k)が計算領域外（上側）の境界条件の処理
            IF (IY.EQ.NY) THEN
*             f(i,NY,k)=0として処理
              AT4(I) = 0.0D0
*             AT4*f(i,NY,k)=B4とする
              B4 = 0.0D0
*           f(i,j+1,k)が計算領域内の場合
            ELSE
              AT4(I) = 1.0D0/DELTA**2 + 1.0D0/(2.0D0*DELTA)
              A(I,I+1) = AT4(I)
              B4=X3(IX,IY+1,IZ)*( 1.0D0/DELTA**2+1.0D0/(2.0D0*DELTA) )
            END IF
*           f(i,j,k-1)が計算領域外（後面）の境界条件の処理
            IF (IZ.EQ.1) THEN
*             f(i,j,0)=0として処理
              AT6(I) = 0.0D0
*             AT6*f(i,j,0)=B6とする
              B6 = 0.0D0
*           f(i,j,k-1)が計算領域内の場合
            ELSE
              AT6(I) = 1.0D0/DELTA**2 - 1.0D0/(2.0D0*DELTA)
              A(I,I-NXY) = AT6(I)
              B6=X3(IX,IY,IZ-1)*( 1.0D0/DELTA**2-1.0D0/(2.0D0*DELTA) )
            END IF
*           f(i,j,k+1)が計算領域外（前面）の境界条件の処理
            IF (IZ.EQ.NZ) THEN
*             f(i,j,NZ)=0として処理
              AT7(I) = 0.0D0
*             AT7*f(i,j,NZ)=B7とする
              B7 = 0.0D0
*           f(i,j,k+1)が計算領域内の場合
            ELSE
              AT7(I) = 1.0D0/DELTA**2 + 1.0D0/(2.0D0*DELTA)
              A(I,I+NXY) = AT7(I)
              B7=X3(IX,IY,IZ+1)*( 1.0D0/DELTA**2 + 1.0D0/(2.0D0*DELTA) )
            END IF
*           Bの計算
            B(I) = B1 + B2 + B3 + B4 + B5 + B6 + B7
            I = I + 1
   80     CONTINUE
   70   CONTINUE
   60 CONTINUE
* 戻り点 : 解法の選択が無効なときの戻り点
  710 CONTINUE
* 各サブルーチンを用いて計算する
      WRITE (6,2000)
 2000 FORMAT ( '  1 : LU Decomposition '/,
     $         '  2 : Gauss Elimination '/,
     $         '  3 : Gauss Elimination without Pivotting '/,
     $         '  4 : Gauss Jordan '/,
     $         '  5 : Gauss Elimination for Band Matrix '/,
     $         '  6 : Jacobi '/,
     $         '  7 : point-SOR '/,
     $         '  8 : point-SOR for Band Matrix '/,
     $         '  9 : line-SOR for Band Matrix '/,
     $         ' 10 : Conjugate Residual '/,
     $         ' 11 : Conjugate Residual for Band Matrix '/,
     $         ' 12 : Bi-CGSTAB '/,
     $         ' 13 : Bi-CGSTAB for Band Matrix '/,
     $         ' 14 : END '/,
     $         ' Input Number ---> ')
      READ (5,*) IJ
*   直接法
      IF (IJ.EQ.1) THEN
        CALL DLU (A,B,X1,NE,
     $            W1,IPV)
      ELSE IF (IJ.EQ.2) THEN
        CALL DFGP (A,B,X1,NE,
     $             IPV)
      ELSE IF (IJ.EQ.3) THEN
        CALL DFG (A,B,X1,NE)
      ELSE IF (IJ.EQ.4) THEN
        CALL DGJP (A,B,X1,NE,
     $             IPV)
      ELSE IF (IJ.EQ.5) THEN
        CALL DGBND (AT1,AT2,AT3,AT4,AT5,AT6,AT7,B,X1,NY,NXY,NE,
     $              AGB)
*   反復法
      ELSE IF (IJ.EQ.6) THEN
        CALL HJACOB (A,B,X1,NE,
     $               W1)
      ELSE IF (IJ.EQ.7) THEN
        CALL HSOR (A,B,X1,NE)
      ELSE IF (IJ.EQ.8) THEN
        CALL SORBND (AT1,AT2,AT3,AT4,AT5,AT6,AT7,B,X1,NY,NXY,NE)
      ELSE IF (IJ.EQ.9) THEN
        CALL LBLBND (AT1,AT2,AT3,AT4,AT5,AT6,AT7,B,X1,NX,NY,NZ,NE,NXY,
     $               W1,W2,W3,W4,W5,W6,W7,W8)
*   クリロフ部分空間法
      ELSE IF (IJ.EQ.10) THEN
        CALL KCR (A,B,X1,NE,
     $            W1,W2,W3,W4)
      ELSE IF (IJ.EQ.11) THEN
        CALL KCRBND (AT1,AT2,AT3,AT4,AT5,AT6,AT7,B,X1,NY,NXY,NE,
     $               W1,W2,W3,W4)
      ELSE IF (IJ.EQ.12) THEN
        CALL KBICG (A,B,X1,NE,
     $              W1,W2,W3,W4,W5,W6)
      ELSE IF (IJ.EQ.13) THEN
        CALL KBIBND (AT1,AT2,AT3,AT4,AT5,AT6,AT7,B,X1,NY,NXY,NE,
     $               W1,W2,W3,W4,W5,W6)
      ELSE IF (IJ.EQ.14) THEN
        GO TO 900
      ELSE
        GO TO 710
      END IF
      GO TO 700
  900 CONTINUE
* 計算終了
      STOP
      END
*********************************************************************
*                  直接法 : 解法 1 - 5
*********************************************************************
*    LU分解による非対称行列 A を含む線形システム解法サブルーチン    *
*    (部分軸選択付)        AX=B                                     *
*    LU分解では次のような2段階で解を求める                          *
*    第1段階 : LY = B から Y を求める                               *
*    第2段階 : UX = Y から X を求める                               *
*    NX : x方向格子分割数; NY : y方向格子分割数; NE : 総格子点数    *
*    A(NE,NE) : AX=B のマトリックス A                               *
*    B(NE) : 線形システム AX=B のマトリックスB                      *
*    X(NE) : 線形システム AX=B のマトリックスX                      *
*            これは第2段階の UX=Y の解 X                            *
*            (このサブルーチンでこれを求める)                       *
*    Y(NE) : LY=B の解 Y                                            *
*    IPV(NE) : 部分軸選択用の配列                                   *
*             全部で NE 回の消去操作を行うが，各消去段階で行う除算  *
*             について分母が最も大きくなるような列の番号が入る      *
*********************************************************************
      SUBROUTINE DLU (A,B,X,NE,
     $                Y,IPV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE),Y(NE)
      DIMENSION IPV(NE)
* 配列のゼロクリア
      DO 10 I = 1,NE
        Y(I)  = 0.0D0
   10 CONTINUE
* 部分軸選択用配列の初期設定
      DO 20 I = 1, NE
        IPV(I) = I
   20 CONTINUE
* 部分軸選択による特異性の判定値
      EPS = 1.0D-50
      DO 30 K = 1, NE
        L = K
        APV = ABS( A(IPV(L),K) )
        DO 40 I = K+1, NE
*         部分軸選択
          IF ( ABS( A(IPV(I),K) ). GT. APV ) THEN
            L = I
            APV = ABS( A(IPV(L),K) )
          END IF
   40   CONTINUE
*       部分軸選択を行った方がよいと判断したら入れ替えを行う
        IF (L.NE.K) THEN
          IPVEX = IPV(K)
          IPV(K) = IPV(L)
          IPV(L) = IPVEX
        END IF
*       部分軸選択を行っても計算不可能なとき : 行列は特異(singular)
        IF ( ABS( A(IPV(K),K) ). LE. EPS ) THEN
          WRITE (*,2000) K
 2000     FORMAT (' Matrix is singular at K = ', I4)
          STOP
        END IF
* U を求める
        A(IPV(K),K) = 1.0D0 / A(IPV(K),K)
        DO 50 I = K+1, NE
          A(IPV(I),K) = A(IPV(I),K) * A(IPV(K),K)
          DO 60 J = K+1, NE
            A(IPV(I),J) = A(IPV(I),J) - A(IPV(I),K)*A(IPV(K),J)
   60     CONTINUE
   50   CONTINUE
   30 CONTINUE
* 第1段階 LY = B を解く
      Y(1) = B(IPV(1))
      DO 70 I = 2, NE
        T = B(IPV(I))
        DO 80 J = 1, I-1
          T = T - A(IPV(I),J) * Y(J)
   80   CONTINUE
        Y(I) = T
   70 CONTINUE
* 第2段階 UX = Y を解く
      X(NE) = Y(NE) * A(IPV(NE),NE)
      DO 90 I = NE-1, 1, -1
        T = Y(I)
        DO 100 J = I+1, NE
          T = T - A(IPV(I),J) * X(J)
  100   CONTINUE
        X(I) = T * A(IPV(I),I)
   90 CONTINUE
* 結果出力
      DO 110 I = 1,NE
        WRITE(6,2010) I,X(I)
 2010   FORMAT(' X_1_( ',I3,' ) = ',1PD13.5)
  110 CONTINUE
      RETURN
      END
*********************************************************************
*    Gaussの消去法(部分軸選択付)による非対称行列Aを含む線形システム *
*    解法サブルーチン    AX=B                                       *
*    NX : x方向格子分割数; NY : y方向格子分割数; NE : 総格子点数    *
*    A(NE,NE) : AX=B のマトリックス A                               *
*    B(NE) : 線形システム AX=B のマトリックスB                      *
*    X(NE) : 線形システム AX=B のマトリックスX                      *
*            (このサブルーチンでこれを求める)                       *
*    IPV(NE) : 部分軸選択用の配列                                   *
*             全部で NE 回の消去操作を行うが，各消去段階で行う除算  *
*             について分母が最も大きくなるような列の番号が入る      *
*********************************************************************
      SUBROUTINE DFGP (A,B,X,NE,
     $                 IPV)
      IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE)
      DIMENSION IPV(NE)
* 部分軸選択による特異性の判定値
      EPS = 1.0D-50
* 部分軸選択用配列の初期設定
      DO 10 I = 1, NE
        IPV(I) = I
   10 CONTINUE
      DO 20 K = 1, NE
        L = K
        APV = ABS( A(IPV(L),K) )
*       部分軸選択
        DO 30 I = K+1, NE
          IF ( ABS( A(IPV(I),K) ). GT. APV ) THEN
            L = I
            APV = ABS( A(IPV(L),K) )
          END IF
   30   CONTINUE
*       部分軸選択を行った方がよいと判断したら入れ替えを行う
        IF (L.NE.K) THEN
          IPVEX = IPV(K)
          IPV(K) = IPV(L)
          IPV(L) = IPVEX
        END IF
*       部分軸選択を行っても計算不可能なとき : 行列は特異
        IF ( ABS( A(IPV(K),K) ). LE. EPS ) THEN
          WRITE (*,2000) K
 2000     FORMAT (' Matrix is singular at K= ', I4)
          STOP
        END IF
* 前進消去
        DO 40 I = K+1, NE
          A(IPV(K),I)=A(IPV(K),I)/A(IPV(K),K)
   40   CONTINUE
        B(IPV(K)) = B(IPV(K))/A(IPV(K),K)
        DO 50 I = K+1, NE
          DO 60 J = K+1, NE
            A(IPV(I),J)=A(IPV(I),J)-A(IPV(I),K)*A(IPV(K),J)
   60     CONTINUE
          B(IPV(I))=B(IPV(I))-A(IPV(I),K)*B(IPV(K)) 
   50   CONTINUE
   20 CONTINUE
* 後退代入
      X(NE) = B(IPV(NE))
      DO 70 I = NE-1, 1, -1
        T = B(IPV(I))
        DO 80 J = I+1, NE
          T = T - A( IPV(I), J ) * X(J)
   80   CONTINUE
        X(I) = T
   70 CONTINUE
* 結果出力
      DO 90 I = 1,NE
        WRITE(6,2010) I,X(I)
 2010   FORMAT(' X_2_( ',I3,' ) = ',1PD13.5)
   90 CONTINUE
      RETURN
      END
*********************************************************************
*    Gaussの消去法(部分軸選択無)による非対称行列Aを含む線形システム *
*    解法サブルーチン    AX=B                                       *
*    NX : x方向格子分割数; NY : y方向格子分割数; NE : 総格子点数    *
*    A(NE,NE) : AX=B のマトリックス A                               *
*    B(NE) : 線形システム AX=B のマトリックスB                      *
*    X(NE) : 線形システム AX=B のマトリックスX ---> これを求める    *
*********************************************************************
      SUBROUTINE DFG (A,B,X,NE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE)
* 特異性の判定値
      EPS = 1.0D-50
* 前進消去
      DO 10 K = 1, NE
*       除算を行う係数が小さく計算不可能なとき : 行列は特異
        IF ( ABS( A(K, K) ). LE. EPS ) THEN
          WRITE (*,2000) K
 2000     FORMAT (' Matrix is singular at K = ', I4)
          STOP
        END IF
        DO 20 I = K+1, NE
          A(K,I)=A(K,I)/A(K,K)
   20   CONTINUE
        B(K) = B(K)/A(K,K)
        DO 30 I = K+1, NE
          DO 40 J = K+1, NE
            A(I,J)=A(I,J)-A(I,K)*A(K,J)
   40     CONTINUE
          B(I)=B(I)-A(I,K)*B(K) 
   30   CONTINUE
   10 CONTINUE
* 後退代入
      X(NE) = B(NE)
      DO 50 I = NE-1, 1, -1
        T = B(I)
        DO 60 J = I+1, NE
          T = T - A( I, J ) * X(J)
   60   CONTINUE
        X(I) = T
   50 CONTINUE
* 結果出力
      DO 70 I = 1,NE
        WRITE(6,2010) I,X(I)
 2010   FORMAT(' X_3_( ',I3,' ) = ',1PD13.5)
   70 CONTINUE
      RETURN
      END
*********************************************************************
*    Gauss-Jordan法(部分軸選択付)による非対称行列Aを含む線形システム*
*    解法サブルーチン    AX=B                                       *
*    NX : x方向格子分割数; NY : y方向格子分割数; NE : 総格子点数    *
*    A(NE,NE) : AX=B のマトリックス A                               *
*    B(NE) : 線形システム AX=B のマトリックスB                      *
*    X(NE) : 線形システム AX=B のマトリックスX ---> これを求める    *
*    IPV(NE) : 部分軸選択用の配列                                   *
*             全部で NE 回の消去操作を行うが，各消去段階で行う除算  *
*             について分母が最も大きくなるような列の番号が入る      *
*********************************************************************
      SUBROUTINE DGJP (A,B,X,NE,
     $                 IPV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE)
      DIMENSION IPV(NE)
* 部分軸選択による特異性の判定値
      EPS = 1.0D-50
* 部分軸選択用配列の初期設定
      DO 10 I = 1, NE
        IPV(I) = I
   10 CONTINUE
      DO 20 K = 1, NE
        L = K
        AIPV = ABS( A(IPV(L),K) )
*       部分軸選択
        DO 30 I = K+1, NE
          IF ( ABS( A(IPV(I),K) ). GT. AIPV ) THEN
            L = I
            AIPV = ABS( A(IPV(L),K) )
          END IF
   30   CONTINUE
*       部分軸選択を行った方がよいと判断したら入れ替えを行う
        IF (L.NE.K) THEN
          IPVEX = IPV(K)
          IPV(K) = IPV(L)
          IPV(L) = IPVEX
        END IF
*       部分軸選択を行っても計算不可能なとき
        IF ( ABS( A(IPV(K),K) ). LE. EPS ) THEN
          WRITE (*,2000) K
 2000     FORMAT (' Matrix is singular at K= ', I4)
          STOP
        END IF
* 消去過程
        DO 40 I = K+1, NE
          A(IPV(K),I)=A(IPV(K),I)/A(IPV(K),K)
   40   CONTINUE
        B(IPV(K)) = B(IPV(K))/A(IPV(K),K)
        DO 50 I = 1, NE
          IF (I.EQ.K) GO TO 50
          DO 60 J = K+1, NE
            A(IPV(I),J)=A(IPV(I),J)-A(IPV(I),K)*A(IPV(K),J)
   60     CONTINUE
          B(IPV(I))=B(IPV(I))-A(IPV(I),K)*B(IPV(K)) 
   50   CONTINUE
   20 CONTINUE
* 上の消去過程が終了すると，IX=Bとなり解はすぐ求まる
      DO 70 I = 1, NE
        X(I) = B(IPV(I))
   70 CONTINUE
* 結果出力
      DO 80 I = 1,NE
        WRITE(6,2010) I,X(I)
 2010   FORMAT(' X_4_( ',I3,' ) = ',1PD13.5)
   80 CONTINUE
      RETURN
      END
*********************************************************************
*    3次元(7点)差分用バンドマトリックス解法サブルーチン AX=B        *
*    3次元(7点)差分にて得られた規則的非対称行列Aを含んだ線形システム*
*    をGaussの消去法を用いて解くサブルーチン．(部分軸選択無)        *
*    高速解法のためにバンドマトリックス用にしてある．               *
*［変数の説明］                                                     *
*    NX:x方向格子分割数; NY:y方向格子分割数; NZ:z方向格子分割数     *
*    NXY:NX*NY; NE : 総格子点数=NX*NY*NZ                            *
*［配列の説明］                                                     *
*    AT1(NE),AT2(NE),AT3(NE),AT4(NE),AT5(NE),AT6(NE),AT7(NE)        *
*     i=3, j=3, k=3 (m=59)における，AT1,AT2,AT3,AT4,AT5,AT6,AT7     *
*                            k=1                                    *
*        +------------------/-----+                                 *
* (NY)5  |    |    |    |  / |    |                                 *
*        +---------------AT6------+                                 *
*     4  |    |    |AT4 /    |    |    m=(k-1)*(NX*NY)+(i-1)*NY+j   *
*        +-------------/----------+                                 *
*     3  |    | AT1|AT3 |AT5 |    |    AT1(m) -> A(-NY ,m)          *
*        +---------/--------------+    AT2(m) -> A(-1  ,m)          *
*     2  |    |   /|AT2 |    |    |    AT3(m) -> A( 0  ,m)          *
*        +-----AT6----------------+    AT4(m) -> A(+1  ,m)          *
*     1  |      /  |    |    |    |    AT5(m) -> A(+NY ,m)          *
*        +-----/------------------+    AT6(m) -> A(-NXY,m)          *
*     j   i=1 /   2    3    4    5(NX) AT7(m) -> A(+NXY,m)          *
*            /                         B(m)                         *
*         k=5(NZ)                      NXY = NX * NY                *
*                                                                   *
*    A(-NXY:NXY,NE) : 3次元差分近似による規則的非対称行列           *
*      -NXY:NXY->m=59のときAT1からAT7までの通し番号 m は，          *
*      AT1は m"-NY"; AT2は m"-1 "; AT3は m"+0 "; AT4は m"+1 "       *
*      AT5は m"+NY"; AT6は m"-(NX*NY)"; AT7は m"+(NX*NY)"           *
*                  となる．この"と"で囲まれた値が-NXY:NXYである．   *
*                  なお，これ以外のATはゼロである．                 *
*      NE->上述のようなkは全部で1からNE                             *
*    B(NE) : 線形システム AX=B のマトリックスB                      *
*    X(NE) : 線形システム AX=B のマトリックスX  ---> これを求める   *
*********************************************************************
      SUBROUTINE DGBND (AT1,AT2,AT3,AT4,AT5,AT6,AT7,B,X,NY,NXY,NE,
     $                  A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AT1(NE),AT2(NE),AT3(NE),AT4(NE),AT5(NE),AT6(NE),AT7(NE)
      DIMENSION B(NE),X(NE),A(-NXY:NXY,NE)
* マトリックスAのゼロクリア
      DO 10 INE = 1,NE
        DO 20 I = -NXY,NXY
          A(I,INE) = 0.0D0
   20   CONTINUE
   10 CONTINUE
* AT1からAT7までをAに格納
      DO 30 INE = 1,NE
        A(- NY,INE) = AT1(INE)
        A(  -1,INE) = AT2(INE)
        A(   0,INE) = AT3(INE)
        A(   1,INE) = AT4(INE)
        A(  NY,INE) = AT5(INE)
        A(-NXY,INE) = AT6(INE)
        A( NXY,INE) = AT7(INE)
   30 CONTINUE
*前進消去
      DO 40 I = 1,NE-1
        IF ( I.LE.NE-NXY ) THEN
          DO 50 J = 1,NXY
            AA = A(-J,I+J)/A(0,I)
            B(I+J) = B(I+J) - B(I)*AA
            N = 1
            DO 60 K = -J+1,NXY-J
              A(K,I+J) = A(K,I+J)-A(N,I)*AA
              N = N + 1
   60       CONTINUE
   50     CONTINUE
        ELSE
          DO 70 J = 1,NE-I
            AA = A(-J,I+J)/A(0,I)
            B(I+J) = B(I+J) - B(I)*AA
            N = 1
            DO 80 K = -J+1,NE-I-J
              A(K,I+J) = A(K,I+J)-A(N,I)*AA
              N = N + 1
   80       CONTINUE
   70     CONTINUE
        END IF
   40 CONTINUE
*係数行列の特異性を判定
      IF ( DABS(A(0,NE)).LE.1.0D-20 ) THEN
        WRITE (6,*) ' Matrix is singular : |A(0,NE)| < 1E-20 '
      END IF
*後退代入
      X(NE) = B(NE) / A(0,NE)
      DO 90 I = NE-1,1,-1
        S = 0.0D0
        IF ( I.GT.NE-NXY ) THEN
          DO 100 N = 1,NE-I
            S = S + A(N,I)* X(I+N)
  100     CONTINUE
          X(I) = ( B(I)-S ) / A(0,I)
        ELSE
          DO 110 N = 1,NXY
            S = S + A(N,I)* X(I+N)
  110     CONTINUE
          X(I) = ( B(I)-S ) / A(0,I)
        END IF
   90 CONTINUE
* 結果出力
      DO 120 I = 1,NE
        WRITE(6,2000) I,X(I)
 2000   FORMAT(' X_5_( ',I3,' ) = ',1PD13.5)
  120 CONTINUE
      RETURN
      END
*********************************************************************
*                  反復法 : 解法 6 - 9
*********************************************************************
*    Jacobi法による非対称行列 A を含む線形システム解法サブルーチン  *
*                        AX=B                                       *
*    NX : x方向格子分割数; NY : y方向格子分割数; NE : 総格子点数    *
*    A(NE,NE) : AX=B のマトリックス A                               *
*    B(NE) : 線形システム AX=B のマトリックスB                      *
*    X(NE) : 旧値 ; XN(NE) : 新値 ---> これを求める                 *
*********************************************************************
      SUBROUTINE HJACOB (A,B,XN,NE,
     $                   X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE),XN(NE)
* 最大繰り返し回数(NITR)と収束判定値(EITR)の設定
      NITR = 200
      EITR = 1.0D-9
* X のゼロクリア
      DO 10 I = 1,NE
        X(I) = 0.0D0
   10 CONTINUE
* Bの2乗ノルムの計算
      BNORM = 0.0D0
      DO 20 I = 1,NE
        BNORM = BNORM + B(I)**2
   20 CONTINUE
* ITR : 繰り返し回数のカウンタ
      ITR = 0
* 繰り返しのための戻り点
  700 CONTINUE
* 先に得られた値を旧値に設定し直す
      ITR = ITR + 1
      DO 30 I = 1,NE
        X(I) = XN(I)
   30 CONTINUE
* 新値の計算
      DO 40 I = 1, NE
        SUM = 0.0D0
        DO 50 J = 1, NE
          IF (I.EQ.J) GO TO 50
          SUM = SUM + A(I,J)*X(J)
   50   CONTINUE
        XN(I) = ( B(I)-SUM )/A(I,I)
   40 CONTINUE
* 収束判定のためのノルムの計算
*     RNORM : 残差ベクトル(R=AX-B)の2乗ノルム
      RNORM = 0.0D0
      DO 60 I = 1,NE
        SUM = 0.0D0
        DO 70 J = 1,NE
          SUM = SUM + A(I,J)*XN(J)
   70   CONTINUE
        RNORM = RNORM + (B(I)-SUM)**2
   60 CONTINUE
*     残差の計算 ZANSA = || R || / || B ||
      ZANSA = DSQRT(RNORM/BNORM)
* 収束判定
      IF (ITR.LE.NITR) THEN
        IF (ZANSA.GT.EITR) THEN
          GO TO 700
        ELSE
          WRITE (6,*) 'Converged : Total ITR = ',ITR
        END IF
      ELSE
        WRITE (6,*) ' Not converged ! '
      END IF
* 結果出力
      DO 80 I = 1,NE
        WRITE(6,2000) I,XN(I)
 2000   FORMAT(' X_6_( ',I3,' ) = ',1PD13.5)
   80 CONTINUE
      RETURN
      END
*********************************************************************
*   point-SOR法による非対称行列 A を含む線形システム解法サブルーチン*
*    (緩和係数OMGを1とするとGauss-Seidel法)  AX=B                   *
*    NX : x方向格子分割数; NY : y方向格子分割数; NE : 総格子点数    *
*    A(NE,NE) : AX=B のマトリックス A                               *
*    B(NE) : 線形システム AX=B のマトリックスB                      *
*    X(NE) : マトリックス X の1次元配列                             *
*********************************************************************
      SUBROUTINE HSOR (A,B,X,NE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE)
* 最大繰り返し回数(NITR)と収束判定値(EITR)の設定
      NITR = 200
      EITR = 1.0D-9
* 緩和係数の設定 ( OMG=1.0D0ならGauss-Seidel法 )
      WRITE (6,*) ' Input OMG ( if Gauss-Seidel, OMG=1.0D0 ) ---> '
      READ (5,*) OMG
* Bの2乗ノルムの計算
      BNORM = 0.0D0
      DO 10 I = 1,NE
        BNORM = BNORM + B(I)**2
   10 CONTINUE
* ITR : 繰り返し回数のカウンタ
      ITR = 0
* 繰り返しのための戻り点
  700 CONTINUE
      ITR = ITR + 1
      DO 20 I = 1, NE
        XOLD = X(I)
        SUM = 0.0D0
        DO 30 J = 1, NE
          IF (I.EQ.J) GO TO 30
          SUM = SUM + A(I,J)*X(J)
   30   CONTINUE
        XNEW = ( B(I)-SUM )/A(I,I)
        X(I) = XOLD + OMG * ( XNEW - XOLD )
   20 CONTINUE
* 収束判定のためのノルムの計算
*     RNORM : 残差ベクトル(R=AX-B)の2乗ノルム
      RNORM = 0.0D0
      DO 40 I = 1,NE
        SUM = 0.0D0
        DO 50 J = 1,NE
          SUM = SUM + A(I,J)*X(J)
   50   CONTINUE
        RNORM = RNORM + (B(I)-SUM)**2
   40 CONTINUE
*     残差の計算 ZANSA = || R || / || B ||
      ZANSA = DSQRT(RNORM/BNORM)
* 収束判定
      IF (ITR.LE.NITR) THEN
        IF (ZANSA.GT.EITR) THEN
          GO TO 700
        ELSE
          WRITE (6,*) 'Converged : Total ITR = ',ITR
        END IF
      ELSE
        WRITE (6,*) ' Not converged ! '
      END IF
* 結果出力
      DO 60 I = 1,NE
        WRITE(6,2000) I,X(I)
 2000   FORMAT(' X_7_( ',I3,' ) = ',1PD13.5)
   60 CONTINUE
      RETURN
      END
*********************************************************************
*  point-SOR 法による非対称行列 A を含む線形システム解法サブルーチン*
*  (バンドマトリックス用)                                           *
*    線形システム --->  A_{i,j} X_{i} = B_{i}                       *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE)                *
*    B : 既知ベクトル                                               *
*    X : 未知ベクトル ---> これを求める                             *
*［変数の説明］                                                     *
*    NX:x方向格子分割数; NY:y方向格子分割数; NE:総格子点数= NX*NY   *
*    NITR : 最大反復回数 ; EITR : 収束判定値                        *
*********************************************************************
      SUBROUTINE SORBND (A1,A2,A3,A4,A5,A6,A7,B,X,NY,NXY,NE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE),
     $          B(NE),X(NE)
* 最大繰り返し回数(NITR)と収束判定値(EITR)の設定
      NITR = 200
      EITR = 1.0D-9
* 緩和係数の設定 ( OMG=1.0D0ならGauss-Seidel法 )
      WRITE (6,*) ' Input OMG ( if Gauss-Seidel, OMG=1.0D0 ) ---> '
      READ (5,*) OMG
* Bの2乗ノルムの計算
      BNORM = 0.0D0
      DO 10 I = 1,NE
        BNORM = BNORM + B(I)**2
   10 CONTINUE
* ITR : 繰り返し回数のカウンタ
      DO 20 ITR = 1,NITR
*       RNORM : 残差ベクトル(R=AX-B)の2乗ノルム
        RNORM = 0.0D0
*       A3,A4,A5,A7の範囲
        I=1
          XOLD = X(I)
          SUM =  A4(I)*X(I+1)+A5(I)*X(I+NY)+A7(I)*X(I+NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
*       A2,A3,A4,A5,A7の範囲
        DO 30 I=2,NY
          XOLD = X(I)
          SUM =  A2(I)*X(I-1)+A4(I)*X(I+1)+A5(I)*X(I+NY)+A7(I)*X(I+NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
   30   CONTINUE
*       A1 - A5, A7 の範囲
        DO 40 I=NY+1,NXY
          XOLD = X(I)
          SUM = A1(I)*X(I-NY)+A2(I)*X(I-1)+A4(I)*X(I+1)+A5(I)*X(I+NY)
     $         +A7(I)*X(I+NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
   40   CONTINUE
*       A6, A1 - A5, A7 の範囲
        DO 50 I=NXY+1,NE-NXY
          XOLD = X(I)
          SUM = A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A4(I)*X(I+1)
     $         +A5(I)*X(I+NY)+A7(I)*X(I+NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
   50   CONTINUE
*       A6, A1 - A5 の範囲
        DO 60 I=NE-NXY+1,NE-NY
          XOLD = X(I)
          SUM = A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A4(I)*X(I+1)
     $         +A5(I)*X(I+NY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
   60   CONTINUE
*       A6, A1 - A4 の範囲
        DO 70 I=NE-NY+1,NE-1
          XOLD = X(I)
          SUM = A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A4(I)*X(I+1)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
   70   CONTINUE
*       A6, A1 - A3 の範囲
        I=NE
          XOLD = X(I)
          SUM = A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
*       RNORM : 残差ベクトル(R=AX-B)の2乗ノルム
        RNORM= 0.0D0
        DO 80 I = 1,NE
          SUM = 0.0D0
          IF (I.EQ.1) THEN
            SUM=A3(I)*X(I)+A4(I)*X(I+1)+A5(I)*X(I+NY)+A7(I)*X(I+NXY)
          ELSE IF (I.GE.2. AND. I.LE.NY) THEN
            SUM=A2(I)*X(I-1)+A3(I)*X(I)+A4(I)*X(I+1)+A5(I)*X(I+NY)
     $         +A7(I)*X(I+NXY)
          ELSE IF (I.GE.NY+1. AND. I.LE.NXY) THEN
            SUM=A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)+A5(I)*X(I+NY)+A7(I)*X(I+NXY)
          ELSE IF (I.GE.NXY+1. AND. I.LE.NE-NXY) THEN
            SUM=A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)+A5(I)*X(I+NY)+A7(I)*X(I+NXY)
          ELSE IF (I.GE.NE-NXY+1. AND. I.LE.NE-NY) THEN
            SUM=A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)+A5(I)*X(I+NY)
          ELSE IF (I.GE.NE-NY+1. AND. I.LE.NE-1) THEN
            SUM=A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)
          ELSE IF (I.EQ.NE) THEN
            SUM=A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
          END IF
          RNORM= RNORM + (B(I)-SUM)**2
   80   CONTINUE
*       残差の計算 ZANSA = || R || / || B ||
        ZANSA = DSQRT(RNORM/BNORM)
*       収束判定
        IF (ZANSA.LE.EITR) THEN
          WRITE (6,*) 'Converged : Total ITR = ',ITR
          GO TO 700
        END IF
   20 CONTINUE
* NITRまで計算しても収束せず
      WRITE (6,*) ' Not converged ! '
* 収束と判定されたときの分岐点
  700 CONTINUE
* 結果出力
      DO 90 I = 1,NE
        WRITE(6,2000) I,X(I)
 2000   FORMAT(' X_8_( ',I3,' ) = ',1PD13.5)
   90 CONTINUE
      RETURN
      END
*********************************************************************
*  line-SOR 法による非対称行列 A を含む線形システム解法サブルーチン *
*  (バンドマトリックス用)                                           *
*    線形システム --->  A_{i,j} X_{i} = B_{i}                       *
* A_{i,j}--->AT1(NE),AT2(NE),AT3(NE),AT4(NE),AT5(NE),AT6(NE),AT7(NE)*
*    B -> BX(NE) : 既知ベクトル                                     *
*    X -> XN(NE) : 未知ベクトル ---> これを求める                   *
* [トーマス法のための係数行列]                                      *
*    x,y,z方向に陰的に離散化された結果を以下のように表す．          *
*  |B(1) C(1) 0    0 ...             |X(1)   |  |D(1)   |           *
*  |A(2) B(2) C(2) 0 ...             |X(2)   |  |D(2)   |           *
*  |0    A(3) B(3) C(3) 0 ...        |X(3)   |= |D(3)   |           *
*  |              ...                |...    |  |...    |           *
*  |0   ...   A(NE-1) B(NE-1) C(NE-1)|X(NE-1)|  |D(NE-1)|           *
*  |0    0    0 ....  A(NE  ) B(NE)  |X(NE)  |  |D(NE)  |           *
*［変数の説明］                                                     *
*    NX:x方向格子分割数; NY:y方向格子分割数; NZ:z方向格子分割数     *
*    NE:総格子点数=NX*NY; NXY=NX*NY                                 *
*    NITR : 最大反復回数; EITR : 収束判定値                         *
*    OMG : 緩和係数．1.0で十分．                                    *
*           注意：Point-SORと異なり，あまり大きくしすぎると発散する *
* [配列の説明]                                                      *
* XN...各方向への掃引後のX(番号付けは不変)                          *
*      はじめにこのサブルーチンへ渡されるXでもある                  *
* X1...各方向への掃引後のX(番号付けは軸方向に異なる)                *
* XO...各方向への掃引前のX(番号付けは不変)                          *
*********************************************************************
      SUBROUTINE LBLBND (AT1,AT2,AT3,AT4,AT5,AT6,AT7,
     $                 BX,XN,NX,NY,NZ,NE,NXY,
     $                 X1,XO,A,B,C,D,U,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AT1(NE),AT2(NE),AT3(NE),AT4(NE),AT5(NE),
     $          AT6(NE),AT7(NE)
      DIMENSION XN(NE),X1(NE),BX(NE),XO(NE)
      DIMENSION A(NE),B(NE),C(NE),D(NE),U(NE),Y(NE)
* 最大繰り返し回数(NITR)と収束判定値(EITR)の設定
      NITR = 200
      EITR = 1.0D-9
* 緩和係数の設定 ( OMG=1.0D0ならGauss-Seidel法 )
      WRITE (6,*) ' Input OMG ( if Gauss-Seidel, OMG=1.0D0 ) ---> '
      READ (5,*) OMG
* Bの2乗ノルムの計算
      BNORM = 0.0D0
      DO 10 I = 1,NE
        BNORM = BNORM + BX(I)**2
   10 CONTINUE
      DO 20 K=1,NITR
*       x 軸方向への掃引 : トーマス法による
        INX = 1
        DO 100 IZ = 1,NZ
          DO 110 IY = 1,NY
            DO 120 IX = 1,NX
            INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*           トーマス法のための係数A,B,C,Dの設定 : XNは最新のX
            A(INX) = AT1(INZ)
            B(INX) = AT3(INZ)
            C(INX) = AT5(INZ)
            D(INX) = BX(INZ)
            IF (INZ-1.GE.1) THEN
              D(INX)=D(INX)-AT2(INZ)*XN(INZ-1)
            END IF
            IF (INZ+1.LE.NE) THEN
              D(INX)=D(INX)-AT4(INZ)*XN(INZ+1)
            END IF
            IF (INZ-NXY.GE.1) THEN
              D(INX)=D(INX)-AT6(INZ)*XN(INZ-NXY)
            END IF
            IF (INZ+NXY.LE.NE) THEN
              D(INX)=D(INX)-AT7(INZ)*XN(INZ+NXY)
            END IF
*           トーマス法で答えを求める前のXNをXOに保存
            XO(INZ) = XN(INZ)
            INX = INX + 1
  120       CONTINUE
  110     CONTINUE
  100   CONTINUE
*       Ly=b を解く
        U(1) = C(1) / B(1)
        DO 130 J = 2,NE-1
          U(J) = C(J) / ( B(J)-A(J)*U(J-1) )
  130   CONTINUE
        Y(1) = D(1) / B(1)
        DO 140 J = 2,NE
          Y(J) = ( D(J)-A(J)*Y(J-1) ) / ( B(J)-A(J)*U(J-1) )
  140   CONTINUE
*       Ux=y を解く
        X1(NE) = Y(NE)
        DO 150 J = NE-1,1,-1
          X1(J) = Y(J) - U(J)*X1(J+1)
  150   CONTINUE
        INX = 1
        DO 160 IZ = 1,NZ
          DO 170 IY = 1,NY
            DO 180 IX = 1,NX
              INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*             得られたX1と反復前のXOにより最新のXNを緩和
              XN(INZ)=(1.0D0-OMG)*XO(INZ)+OMG*X1(INX)
              INX = INX + 1
  180       CONTINUE
  170     CONTINUE
  160   CONTINUE
*       y 軸方向への掃引 : トーマス法による
        INY = 1
        DO 200 IZ = 1,NZ
          DO 210 IX = 1,NX
            DO 220 IY = 1,NY
            INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*           トーマス法のための係数A,B,C,Dの設定 : XNは最新のX
            A(INY) = AT2(INZ)
            B(INY) = AT3(INZ)
            C(INY) = AT4(INZ)
            D(INY) = BX(INZ)
            IF (INZ-NY.GE.1) THEN
              D(INY)=D(INY)-AT1(INZ)*XN(INZ-NY)
            END IF
            IF (INZ+NY.LE.NE) THEN
              D(INY)=D(INY)-AT5(INZ)*XN(INZ+NY)
            END IF
            IF (INZ-NXY.GE.1) THEN
              D(INY)=D(INY)-AT6(INZ)*XN(INZ-NXY)
            END IF
            IF (INZ+NXY.LE.NE) THEN
              D(INY)=D(INY)-AT7(INZ)*XN(INZ+NXY)
            END IF
*           トーマス法で答えを求める前のXNをXOに保存
            XO(INZ) = XN(INZ)
            INY = INY + 1
  220       CONTINUE
  210     CONTINUE
  200   CONTINUE
*       Ly=b を解く
        U(1) = C(1) / B(1)
        DO 230 J = 2,NE-1
          U(J) = C(J)/( B(J)-A(J)*U(J-1) )
  230   CONTINUE
        Y(1) = D(1) / B(1)
        DO 240 J = 2,NE
          Y(J) = ( D(J)-A(J)*Y(J-1) ) / ( B(J)-A(J)*U(J-1) )
  240   CONTINUE
*       Ux=y を解く
        X1(NE) = Y(NE)
        DO 250 J = NE-1,1,-1
          X1(J) = Y(J) - U(J)*X1(J+1)
  250   CONTINUE
        INY = 1
        DO 260 IZ = 1,NZ
          DO 270 IX = 1,NX
            DO 280 IY = 1,NY
              INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*             得られたX1と反復前のXOにより最新のXNを緩和
              XN(INZ)=(1.0D0-OMG)*XO(INZ)+OMG*X1(INY)
              INY = INY + 1
  280       CONTINUE
  270     CONTINUE
  260   CONTINUE
*       z 軸方向への掃引 : トーマス法による
        INZZ = 1
        DO 300 IY = 1,NY
          DO 310 IX = 1,NX
            DO 320 IZ = 1,NZ
            INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*           トーマス法のための係数A,B,C,Dの設定 : XNは最新のX
            A(INZZ) = AT6(INZ)
            B(INZZ) = AT3(INZ)
            C(INZZ) = AT7(INZ)
            D(INZZ) = BX(INZ)
            IF (INZ-NY.GE.1) THEN
              D(INZZ)=D(INZZ)-AT1(INZ)*XN(INZ-NY)
            END IF
            IF (INZ+NY.LE.NE) THEN
              D(INZZ)=D(INZZ)-AT5(INZ)*XN(INZ+NY)
            END IF
            IF (INZ-1.GE.1) THEN
              D(INZZ)=D(INZZ)-AT2(INZ)*XN(INZ-1)
            END IF
            IF (INZ+1.LE.NE) THEN
              D(INZZ)=D(INZZ)-AT4(INZ)*XN(INZ+1)
            END IF
*           トーマス法で答えを求める前のXNをXOに保存
            XO(INZ ) = XN(INZ)
            INZZ = INZZ + 1
  320       CONTINUE
  310     CONTINUE
  300   CONTINUE
*       Ly=b を解く
        U(1) = C(1) / B(1)
        DO 330 J = 2,NE-1
          U(J) = C(J)/( B(J)-A(J)*U(J-1) )
  330   CONTINUE
        Y(1) = D(1) / B(1)
        DO 340 J = 2,NE
          Y(J) = ( D(J)-A(J)*Y(J-1) ) / ( B(J)-A(J)*U(J-1) )
  340   CONTINUE
*       Ux=y を解く
        X1(NE) = Y(NE)
        DO 350 J = NE-1,1,-1
          X1(J) = Y(J) - U(J)*X1(J+1)
  350   CONTINUE
        INZZ = 1
        DO 360 IY = 1,NY
          DO 370 IX = 1,NX
            DO 380 IZ = 1,NZ
              INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*             得られたX1と反復前のXOにより最新のXNを緩和
              XN(INZ)=(1.0D0-OMG)*XO(INZ)+OMG*X1(INZZ)
              INZZ = INZZ + 1
  380       CONTINUE
  370     CONTINUE
  360   CONTINUE
*       RNORM : 残差ベクトル(R=AX-B)の2乗ノルム : サブルーチンSORBNDの場合分けを参照
        RNORM= 0.0D0
        DO 40 I = 1,NE
          SUM = 0.0D0
          IF (I.EQ.1) THEN
            SUM=AT3(I)*XN(I)+AT4(I)*XN(I+1)+AT5(I)*XN(I+NY)
     $         +AT7(I)*XN(I+NXY)
          ELSE IF (I.GE.2. AND. I.LE.NY) THEN
            SUM=AT2(I)*XN(I-1)+AT3(I)*XN(I)+AT4(I)*XN(I+1)
     $         +AT5(I)*XN(I+NY)+AT7(I)*XN(I+NXY)
          ELSE IF (I.GE.NY+1. AND. I.LE.NXY) THEN
            SUM=AT1(I)*XN(I-NY)+AT2(I)*XN(I-1)+AT3(I)*XN(I)
     $         +AT4(I)*XN(I+1)+AT5(I)*XN(I+NY)+AT7(I)*XN(I+NXY)
          ELSE IF (I.GE.NXY+1. AND. I.LE.NE-NXY) THEN
            SUM=AT6(I)*XN(I-NXY)+AT1(I)*XN(I-NY)+AT2(I)*XN(I-1)
     $         +AT3(I)*XN(I)+AT4(I)*XN(I+1)+AT5(I)*XN(I+NY)
     $         +AT7(I)*XN(I+NXY)
          ELSE IF (I.GE.NE-NXY+1. AND. I.LE.NE-NY) THEN
            SUM=AT6(I)*XN(I-NXY)+AT1(I)*XN(I-NY)+AT2(I)*XN(I-1)
     $         +AT3(I)*XN(I)+AT4(I)*XN(I+1)+AT5(I)*XN(I+NY)
          ELSE IF (I.GE.NE-NY+1. AND. I.LE.NE-1) THEN
            SUM=AT6(I)*XN(I-NXY)+AT1(I)*XN(I-NY)+AT2(I)*XN(I-1)
     $         +AT3(I)*XN(I)+AT4(I)*XN(I+1)
          ELSE IF (I.EQ.NE) THEN
            SUM=AT6(I)*XN(I-NXY)+AT1(I)*XN(I-NY)+AT2(I)*XN(I-1)
     $         +AT3(I)*XN(I)
          END IF
          RNORM= RNORM + (BX(I)-SUM)**2
   40   CONTINUE
*       残差の計算 ZANSA = || R || / || B ||
        ZANSA = DSQRT(RNORM/BNORM)
*       収束判定
        IF (ZANSA.LE.EITR) THEN
          WRITE (6,*) 'Converged : Total ITR = ',K
          GO TO 700
        END IF
   20 CONTINUE
* NITRまで計算しても収束せず
      WRITE (6,*) ' Not converged ! '
* 収束と判定されたときの分岐点
  700 CONTINUE
* 結果出力
      DO 50 I = 1,NE
        WRITE(6,2000) I,XN(I)
 2000   FORMAT(' X_9_( ',I3,' ) = ',1PD13.5)
   50 CONTINUE
      RETURN
      END
*********************************************************************
*              クリロフ部分空間法 : 解法 10 - 13
*********************************************************************
*    共役残差(Conjugate Residual)法による非対称行列 A を含む        *
*  線形システム解法サブルーチン AX=B                                *
*［変数の説明］                                                     *
*    NX:x方向格子分割数; NY:y方向格子分割数; NE:総格子点数=NX*NY    *
*    NITR : 最大反復回数; EPSP : 収束判定値                         *
*［配列の説明］                                                     *
*    A(NE,NE) : AX=B のフルマトリックス A                           *
*    B(NE) : 線形システム AX=B のベクトル B                         *
*    X(NE) : 線形システム AX=B のベクトル X ---> これを求める       *
*    R(NE) : r_{k} = B - A x_{k}                                    *
*    P(NE) :  p_{k+1} = r_{k+1} + β_{k} p_{k}, p_{0} = r_{0}       *
*    AP(NE) : A * P , AR(NE) : A * R                                *
*********************************************************************
      SUBROUTINE KCR (A,B,X,NE,
     $                R,P,AP,AR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE)
* 作業用配列
      DIMENSION R(NE),P(NE),AP(NE),AR(NE)
* 最大繰り返し回数の設定
      NITR = 200
* 収束判定値(変数 ZANSA がこの値以下になれば NITR 以下でも収束と判定する)
      EITR = 1.0D-9
* 配列のゼロクリア
      DO 10 J=1, NE
        R(J) = 0.0D0
        P(J) = 0.0D0
        AP(J) = 0.0D0
        AR(J) = 0.0D0
   10 CONTINUE
* R に AX を代入
      CALL PROFMV (A,X,R,NE)
* Bのノルムの計算
      BNORM = 0.0D0
      DO 20 I = 1,NE
        BNORM = BNORM + B(I)**2
   20 CONTINUE
* r_{0}とp_{0}(初期値)の設定
      DO 30 I = 1, NE
*       r_{0} = B - AX の計算
        R(I) = B(I) - R(I)
*       p_{0} = r_{0}
        P(I) = R(I)
   30 CONTINUE
* APに A p_{0} を代入(以降のAPは 式(A) で求める)
      CALL PROFMV (A,P,AP,NE)
* 繰り返し計算
      DO 40 K = 1,NITR
*       ( r_{k}, A p_{k} )の計算 => RAP
        CALL PROVV (R,AP,RAP,NE)
*       ( A p_{k}, A p_{k} )の計算 => APAP
        CALL PROVV (AP,AP,APAP,NE)
*       α_{k} = ( r_{k}, A p_{k} ) / ( A p_{k}, A p_{k} )
        ALP = RAP / APAP
        DO 50 I = 1,NE
*         x_{k+1}=x_{k}+α_{k}p_{k}
          X(I) = X(I) + ALP*P(I)
*         r_{k+1}=r_{k}-α_{k}Ap_{k}
          R(I) = R(I) - ALP*AP(I)
   50   CONTINUE
*       RNORM : 残差ベクトル(R=AX-B)の2乗ノルム
        RNORM = 0.0D0
        DO 60 I = 1,NE
          SUM = 0.0D0
          DO 70 J = 1,NE
            SUM = SUM + A(I,J)*X(J)
   70     CONTINUE
          RNORM = RNORM + (B(I)-SUM)**2
   60   CONTINUE
*       残差の計算 ZANSA = || R || / || B ||
        ZANSA = DSQRT(RNORM/BNORM)
*       ZANSA が EITR 以下なら収束とみなして結果を表示して 900 へ
        IF (ZANSA.LE.EITR) THEN
          WRITE (6,*) ' Converged : Total ITR = ',K
          GO TO 900
        ELSE
*         A r_{k+1} の計算 => AR(NE)
          CALL PROFMV (A,R,AR,NE)
*         ( A r_{k+1}, A p_{k} )の計算 => ARAP
          CALL PROVV( AR,AP,ARAP,NE)
*         β_{k} = - ( A r_{k+1}, A p_{k} ) / ( A p_{k}, A p_{k} )
          BETA = - ARAP / APAP
          DO 80 I = 1, NE
*           p_{k+1} = r_{k+1} + β_{k} p_{k}
            P(I) = R(I) + BETA*P(I)
*           A p_{k+1} = A r_{k+1} + β_{k}A p_{k}<-------式(A)
            AP(I) = AR(I) + BETA*AP(I)
   80     CONTINUE
        END IF
   40 CONTINUE
*     NITRまで計算しても収束せず
      WRITE (6,*) ' Not Converged ! '
*     収束と判定されたときの分岐点
  900 CONTINUE
*結果出力
      DO 90 I = 1,NE
        WRITE(6,2000) I,X(I)
 2000   FORMAT(' X_10_( ',I3,' ) = ',1PD13.5)
   90 CONTINUE
      RETURN
      END
*********************************************************************
*    フルマトリックス A とベクトル B の積の計算サブルーチン AB=C    *
*    NE : 総格子点数                                                *
*********************************************************************
      SUBROUTINE PROFMV (A,B,C,NE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),C(NE)
      DO 10 I = 1, NE
        S = 0.0D0
        DO 20 J = 1, NE
          S = S + A(I,J)*B(J)
   20   CONTINUE
        C(I) = S
   10 CONTINUE
      RETURN
      END
*********************************************************************
*    ベクトル A とベクトル B の積の計算サブルーチン                 *
*                        AB=C                                       *
*［変数の説明］                                                     *
*    NE : 総格子点数(ベクトル A,B のサイズ)                         *
*    C  : A と B の積(スカラー)                                     *
*********************************************************************
      SUBROUTINE PROVV(A,B,C,NE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  A(NE),B(NE)
      C = 0.0D0
      DO 10 I=1,NE
        C = C + A(I)*B(I)
   10 CONTINUE
      RETURN
      END
*********************************************************************
*    共役残差(Conjugate Residual)法による非対称行列 A を含む        *
*  線形システム解法サブルーチン   (バンドマトリックス用)            *
*                        AX=B                                       *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE)                *
*  B : 既知ベクトル                                                 *
*  X : 未知ベクトル ---> これを求める                               *
*［変数の説明］                                                     *
*    NX:x方向格子分割数; NY:y方向格子分割数; NE:総格子点数=NX*NY    *
*    NITR : 最大反復回数; EPSP : 収束判定値                         *
*［配列の説明］                                                     *
*    R(NE) : r_{k} = B - A x_{k}                                    *
*    P(NE) :  p_{k+1} = r_{k+1} + β_{k} p_{k}, p_{0} = r_{0}       *
*    AP(NE) : A * P                                                 *
*    AR(NE) : A * R                                                 *
*********************************************************************
      SUBROUTINE KCRBND (A1,A2,A3,A4,A5,A6,A7,B,X,NY,NXY,NE,
     $                   R,P,AP,AR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE),
     $          B(NE),X(NE)
* 作業用配列
      DIMENSION R(NE),P(NE),AP(NE),AR(NE)
* 最大繰り返し回数の設定
      NITR = 200
* 収束判定条件(変数 ZANSA がこの値以下になれば NITR 以下でも収束と判定する)
      EITR = 1.0D-9
* 配列のゼロクリア
      DO 10 J=1, NE
        R(J) = 0.0D0
        P(J) = 0.0D0
        AP(J) = 0.0D0
        AR(J) = 0.0D0
   10 CONTINUE
* R に AX を代入
      CALL PROBMV (A1,A2,A3,A4,A5,A6,A7,X,R,NY,NXY,NE)
* Bのノルムの計算
      BNORM = 0.0D0
      DO 20 I = 1,NE
        BNORM = BNORM + B(I)**2
   20 CONTINUE
* r_{0}とp_{0}(初期値)の設定
      DO 30 I = 1, NE
*       r_{0} = B - AX の計算
        R(I) = B(I) - R(I)
*       p_{0} = r_{0}
        P(I) = R(I)
   30 CONTINUE
* APに A p_{0} を代入(以降のAPは 式(A) で求める)
      CALL PROBMV (A1,A2,A3,A4,A5,A6,A7,P,AP,NY,NXY,NE)
* 繰り返し計算
      DO 40 K = 1,NITR
*       ( r_{k}, A p_{k} )の計算 => RAP
        CALL PROVV (R,AP,RAP,NE)
*       ( A p_{k}, A p_{k} )の計算 => APAP
        CALL PROVV (AP,AP,APAP,NE)
*       α_{k} = ( r_{k}, A p_{k} ) / ( A p_{k}, A p_{k} )
        ALP = RAP / APAP
        DO 50 I = 1,NE
*         x_{k+1}=x_{k}+α_{k}p_{k}
          X(I) = X(I) + ALP*P(I)
*         r_{k+1}=r_{k}-α_{k}Ap_{k}
          R(I) = R(I) - ALP*AP(I)
   50   CONTINUE
*       RNORM : 残差ベクトル(R=AX-B)の2乗ノルム : サブルーチンSORBNDの場合分けを参照
        RNORM= 0.0D0
        DO 60 I = 1,NE
          SUM = 0.0D0
          IF (I.EQ.1) THEN
            SUM=A3(I)*X(I)+A4(I)*X(I+1)+A5(I)*X(I+NY)+A7(I)*X(I+NXY)
          ELSE IF (I.GE.2. AND. I.LE.NY) THEN
            SUM=A2(I)*X(I-1)+A3(I)*X(I)+A4(I)*X(I+1)+A5(I)*X(I+NY)
     $         +A7(I)*X(I+NXY)
          ELSE IF (I.GE.NY+1. AND. I.LE.NXY) THEN
            SUM=A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)+A5(I)*X(I+NY)+A7(I)*X(I+NXY)
          ELSE IF (I.GE.NXY+1. AND. I.LE.NE-NXY) THEN
            SUM=A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)+A5(I)*X(I+NY)+A7(I)*X(I+NXY)
          ELSE IF (I.GE.NE-NXY+1. AND. I.LE.NE-NY) THEN
            SUM=A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)+A5(I)*X(I+NY)
          ELSE IF (I.GE.NE-NY+1. AND. I.LE.NE-1) THEN
            SUM=A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)
          ELSE IF (I.EQ.NE) THEN
            SUM=A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
          END IF
          RNORM= RNORM + (B(I)-SUM)**2
   60   CONTINUE
*       残差の計算 ZANSA = || R || / || B ||
        ZANSA = DSQRT(RNORM/BNORM)
*       ZANSA が EITR 以下なら収束とみなして 700 へ
        IF (ZANSA.LE.EITR) THEN
          WRITE (6,*) ' Converged : Total ITR = ',K
          GO TO 700
        END IF
*       収束せずの場合
*       A r_{k+1} の計算 => AR(NE)
        CALL PROBMV (A1,A2,A3,A4,A5,A6,A7,R,AR,NY,NXY,NE)
*       ( A r_{k+1}, A p_{k} )の計算 => ARAP
        CALL PROVV( AR,AP,ARAP,NE)
*       β_{k} = - ( A r_{k+1}, A p_{k} ) / ( A p_{k}, A p_{k} )
        BETA = - ARAP / APAP
        DO 70 I = 1, NE
*         p_{k+1} = r_{k+1} + β_{k} p_{k}
          P(I) = R(I) + BETA*P(I)
*         A p_{k+1} = A r_{k+1} + β_{k}A p_{k}<-------式(A)
          AP(I) = AR(I) + BETA*AP(I)
   70   CONTINUE
   40 CONTINUE
* NITR まで計算しても収束せず
      WRITE (6,*) ' Not Converged ! '
* 収束と判定されたときの分岐点
  700 CONTINUE
* 結果出力
      DO 80 I = 1,NE
        WRITE(6,2000) I,X(I)
 2000   FORMAT(' X_11_( ',I3,' ) = ',1PD13.5)
   80 CONTINUE
      RETURN
      END
*********************************************************************
*    バンドマトリックス A とベクトル B の積の計算サブルーチン       *
*                        AB=C                                       *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE)                *
*［変数の説明］                                                     *
*    NY : y方向格子分割数                                           *
*                 ---> バンドマトリックスとベクトルの積の計算で使用 *
*    NE : 総格子点数                                                *
*    C  : A と B の積(ベクトル)                                     *
*********************************************************************
      SUBROUTINE PROBMV (A1,A2,A3,A4,A5,A6,A7,B,C,NY,NXY,NE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE),
     $          B(NE),C(NE)
* サブルーチンSORBNDの場合分けを参照
      I=1
        C(I) = A3(I)*B(I)+A4(I)*B(I+1)+A5(I)*B(I+NY)+A7(I)*B(I+NXY)
      DO 10 I=2,NY
        C(I) = A2(I)*B(I-1)+A3(I)*B(I)+A4(I)*B(I+1)+A5(I)*B(I+NY)
     $        +A7(I)*B(I+NXY)
   10 CONTINUE
      DO 20 I=NY+1,NXY
        C(I) = A1(I)*B(I-NY)+A2(I)*B(I-1)+A3(I)*B(I)
     $      +A4(I)*B(I+1)+A5(I)*B(I+NY)+A7(I)*B(I+NXY)
   20 CONTINUE
      DO 30 I=NXY+1,NE-NXY
        C(I) = +A6(I)*B(I-NXY)+A1(I)*B(I-NY)+A2(I)*B(I-1)+A3(I)*B(I)
     $      +A4(I)*B(I+1)+A5(I)*B(I+NY)+A7(I)*B(I+NXY)
   30 CONTINUE
      DO 40 I=NE-NXY+1,NE-NY
        C(I) = +A6(I)*B(I-NXY)+A1(I)*B(I-NY)+A2(I)*B(I-1)+A3(I)*B(I)
     $      +A4(I)*B(I+1)+A5(I)*B(I+NY)
   40 CONTINUE
      DO 50 I=NE-NY+1,NE-1
        C(I) = +A6(I)*B(I-NXY)+A1(I)*B(I-NY)+A2(I)*B(I-1)+A3(I)*B(I)
     $      +A4(I)*B(I+1)
   50 CONTINUE
      I=NE
        C(I) = +A6(I)*B(I-NXY)+A1(I)*B(I-NY)+A2(I)*B(I-1)+A3(I)*B(I)
      RETURN
      END
*********************************************************************
*  Bi-CGSTAB 法による非対称行列 A を含む線形システム解法サブルーチン*
*                        AX=B                                       *
*    B : 既知ベクトル                                               *
*    X : 未知ベクトル ---> これを求める                             *
*［変数の説明］                                                     *
*    NX : x方向格子分割数                                           *
*    NY : y方向格子分割数                                           *
*    NE : 総格子点数 = NX * NY                                      *
*    NITR : 最大反復回数                                            *
*    EPSP : 収束判定条件                                            *
*［配列の説明］                                                     *
*    T(NE) : t_{k} = r_{k} - α_{k} A p_{k}                         *
*    X(NE) : x_{k+1} = x_{k} + α_{k} p_{k} + ξ_{K} t_{k}          *
*    R(NE) : r_{k+1} = t_{k} - ξ_{k} A t_{k}                       *
*            r_{0} = B - A x_{0}                                    *
*    P(NE) : p_{k+1} = r_{k+1} + β_{k} ( p_{k}-ξ_{k} A p_{k} )    *
*            p_{0} = r_{0}                                          *
*    AP(NE) : A * P                                                 *
*    AR(NE) : A * T                                                 *
*********************************************************************
      SUBROUTINE KBICG (A,B,X,NE,
     $                  R,AP,AT,P,S,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE)
* 作業用配列
      DIMENSION R(NE),AP(NE),AT(NE),P(NE),S(NE),T(NE)
* 最大繰り返し回数の設定
      NITR = 200
* 収束判定条件(変数 ZANSA がこの値以下になれば NITR 以下でも収束と判定する)
      EITR = 1.0D-9
* Bのノルムの計算
      BNORM = 0.0D0
      DO 10 I = 1,NE
        BNORM = BNORM + B(I)**2
   10 CONTINUE
* 配列のゼロクリア
      DO 20 J=1, NE
        R(J) = 0.0D0
        AP(J) = 0.0D0
        AT(J) = 0.0D0
        P(J) = 0.0D0
        S(J) = 0.0D0
        T(J) = 0.0D0
   20 CONTINUE
* R に AX を代入
      CALL PROFMV (A,X,R,NE)
*r_{0}とp_{0}(初期値)，そして s=r_{0} の設定
      DO 30 I=1,NE
        R(I) = B(I) - R(I)
        P(I) = R(I)
        S(I) = R(I)
   30 CONTINUE
* 繰り返し計算
      DO 40 J =1,NITR
*       ( s, r_{k} ) の計算 => SR1
        CALL PROVV (S,R,SR1,NE)
*       A p_{k} の計算 => AP(NE)
        CALL PROFMV (A,P,AP,NE)
*       ( s, A p_{k} ) の計算 => SAP
        CALL PROVV (S,AP,SAP,NE)
*       α_{k} = ( s, r_{k} ) / ( s, A p_{k} )
        ALPHA = SR1/SAP
          DO 50 I=1,NE
*           t_{k} = r_{k} - α_{k} A p_{k}
            T(I) = R(I) - ALPHA*AP(I)
   50     CONTINUE
*       A t_{k} の計算 => AT(NE)
        CALL PROFMV (A,T,AT,NE)
*       ( A t_{k}, t_{k} ) の計算 => ATT
        CALL PROVV (AT,T,ATT,NE)
*       ( A t_{k}, A t_{k} ) の計算 => ATAT
        CALL PROVV (AT,AT,ATAT,NE)
*       ξ_{k} = ( A t_{k}, t_{k} ) / ( A t_{k}, A t_{k} )
        XI = ATT/ATAT
        DO 60 I=1,NE
*         x_{k+1} = x_{k} + α_{k} p_{k} + ξ_{k} t_{k}
          X(I) = X(I) + ALPHA*P(I) + XI*T(I)
*         r_{k+1} = t_{k} - ξ_{k} A t_{k}
          R(I) = T(I) - XI*AT(I)
   60   CONTINUE
* RNORM : 残差ベクトル(R=AX-B)の2乗ノルム
      RNORM = 0.0D0
      DO 70 I = 1,NE
        SUM = 0.0D0
        DO 80 M = 1,NE
          SUM = SUM + A(I,M)*X(M)
   80   CONTINUE
        RNORM = RNORM + (B(I)-SUM)**2
   70 CONTINUE
*     残差の計算 ZANSA = || R || / || B ||
        ZANSA = DSQRT(RNORM/BNORM)
*       ZANSA が EITR 以下なら収束とみなして 900 へ
        IF (ZANSA.LE.EITR) THEN
          WRITE (6,*) ' Converged : Total ITR = ',J
          GO TO 900
        END IF
*       収束せずの場合 : β_{k}と p_{k+1} を求めて繰り返し計算
*       ( s, r_{k+1} ) の計算 => SR2
        CALL PROVV (S,R,SR2,NE)
*       β_{k} = ( α_{k} / ξ_{k} ) * ( s, r_{k+1} )/( s, r_{k} )
        BETA = (ALPHA / XI) * (SR2 / SR1)
        DO 90 I=1,NE
*         p_{k+1} = r_{k+1} + β_{k}( p_{k}-ξ_{k} A p_{k} )
          P(I) = R(I) + BETA * ( P(I) - XI*AP(I) )
   90   CONTINUE
   40 CONTINUE
*     NITR まで計算しても収束せず
      WRITE (6,*) ' Not Converged ! '
* 収束と判定されたときの分岐点
  900 CONTINUE
* 結果出力
      DO 100 I = 1,NE
        WRITE(6,2000) I,X(I)
 2000   FORMAT(' X_12_( ',I3,' ) = ',1PD13.5)
  100 CONTINUE
      RETURN
      END
*********************************************************************
*  Bi-CGSTAB 法による非対称行列 A を含む線形システム解法サブルーチン*
*  (バンドマトリックス用)                                           *
*                        AX=B                                       *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE)                *
*    B : 既知ベクトル                                               *
*    X : 未知ベクトル ---> これを求める                             *
*［変数の説明］                                                     *
*    NX : x方向格子分割数                                           *
*    NY : y方向格子分割数                                           *
*    NE : 総格子点数 = NX * NY                                      *
*    NITR : 許容反復回数(in2d.mac)にて設定                          *
*    EPSP : 収束判定条件で用いる値(in2d.mac)にて設定                *
*［配列の説明］                                                     *
*    T(NE) : t_{k} = r_{k} - α_{k} A p_{k}                         *
*    X(NE) : x_{k+1} = x_{k} + α_{k} p_{k} + ξ_{K} t_{k}          *
*    R(NE) : r_{k+1} = t_{k} - ξ_{k} A t_{k}                       *
*            r_{0} = B - A x_{0}                                    *
*    P(NE) : p_{k+1} = r_{k+1} + β_{k} ( p_{k}-ξ_{k} A p_{k} )    *
*            p_{0} = r_{0}                                          *
*    AP(NE) : A * P                                                 *
*    AT(NE) : A * T                                                 *
*********************************************************************
      SUBROUTINE KBIBND (A1,A2,A3,A4,A5,A6,A7,B,X,NY,NXY,NE,
     $                  R,AP,AT,P,S,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE),
     $          B(NE),X(NE)
* 作業用配列
      DIMENSION R(NE),AP(NE),AT(NE),P(NE),S(NE),T(NE)
* 最大繰り返し回数の設定
      NITR = 200
* 収束判定条件(変数 ZANSA がこの値以下になれば NITR 以下でも収束と判定する)
      EITR = 1.0D-9
* Bのノルムの計算
      BNORM = 0.0D0
      DO 10 I = 1,NE
        BNORM = BNORM + B(I)**2
   10 CONTINUE
* 配列のゼロクリア
      DO 20 J=1, NE
        R(J) = 0.0D0
        AP(J) = 0.0D0
        AT(J) = 0.0D0
        P(J) = 0.0D0
        S(J) = 0.0D0
        T(J) = 0.0D0
   20 CONTINUE
* R に AX を代入
      CALL PROBMV (A1,A2,A3,A4,A5,A6,A7,X,R,NY,NXY,NE)
* r_{0}とp_{0}(初期値)，そして s=r_{0} の設定
      DO 30 I=1,NE
        R(I) = B(I) - R(I)
        P(I) = R(I)
        S(I) = R(I)
   30 CONTINUE
* 繰り返し計算
      DO 40 J =1,NITR
*       ( s, r_{k} ) の計算 => SR1
        CALL PROVV (S,R,SR1,NE)
*       A p_{k} の計算 => AP(NE)
        CALL PROBMV (A1,A2,A3,A4,A5,A6,A7,P,AP,NY,NXY,NE)
*       ( s, A p_{k} ) の計算 => SAP
        CALL PROVV (S,AP,SAP,NE)
*       α_{k} = ( s, r_{k} ) / ( s, A p_{k} )
        ALPHA = SR1/SAP
        DO 50 I=1,NE
*         t_{k} = r_{k} - α_{k} A p_{k}
          T(I) = R(I) - ALPHA*AP(I)
   50   CONTINUE
*       A t_{k} の計算 => AT(NE)
        CALL PROBMV (A1,A2,A3,A4,A5,A6,A7,T,AT,NY,NXY,NE)
*       ( A t_{k}, t_{k} ) の計算 => ATT
        CALL PROVV (AT,T,ATT,NE)
*       ( A t_{k}, A t_{k} ) の計算 => ATAT
        CALL PROVV (AT,AT,ATAT,NE)
*       ξ_{k} = ( A t_{k}, t_{k} ) / ( A t_{k}, A t_{k} )
        XI = ATT/ATAT
        DO 60 I=1,NE
*         x_{k+1} = x_{k} + α_{k} p_{k} + ξ_{k} t_{k}
          X(I) = X(I) + ALPHA*P(I) + XI*T(I)
*         r_{k+1} = t_{k} - ξ_{k} A t_{k}
          R(I) = T(I) - XI*AT(I)
   60   CONTINUE
*       RNORM : 残差ベクトル(R=AX-B)の2乗ノルム : サブルーチンSORBNDの場合分けを参照
        RNORM= 0.0D0
        DO 70 I = 1,NE
          SUM = 0.0D0
          IF (I.EQ.1) THEN
            SUM=A3(I)*X(I)+A4(I)*X(I+1)+A5(I)*X(I+NY)+A7(I)*X(I+NXY)
          ELSE IF (I.GE.2. AND. I.LE.NY) THEN
            SUM=A2(I)*X(I-1)+A3(I)*X(I)+A4(I)*X(I+1)+A5(I)*X(I+NY)
     $         +A7(I)*X(I+NXY)
          ELSE IF (I.GE.NY+1. AND. I.LE.NXY) THEN
            SUM=A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)+A5(I)*X(I+NY)+A7(I)*X(I+NXY)
          ELSE IF (I.GE.NXY+1. AND. I.LE.NE-NXY) THEN
            SUM=A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)+A5(I)*X(I+NY)+A7(I)*X(I+NXY)
          ELSE IF (I.GE.NE-NXY+1. AND. I.LE.NE-NY) THEN
            SUM=A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)+A5(I)*X(I+NY)
          ELSE IF (I.GE.NE-NY+1. AND. I.LE.NE-1) THEN
            SUM=A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)
          ELSE IF (I.EQ.NE) THEN
            SUM=A6(I)*X(I-NXY)+A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
          END IF
          RNORM= RNORM + (B(I)-SUM)**2
   70   CONTINUE
*       残差の計算 ZANSA = || R || / || B ||
        ZANSA = DSQRT(RNORM/BNORM)
*       ZANSA が EITR 以下なら収束とみなして 900 へ
        IF (ZANSA.LE.EITR) THEN
          WRITE (6,*) ' Converged : Total ITR = ',J
          GO TO 900
        END IF
*       収束せずの場合
*       ( s, r_{k+1} ) の計算 => SR2
        CALL PROVV (S,R,SR2,NE)
*       β_{k} = ( α_{k} / ξ_{k} ) * ( s, r_{k+1} )/( s, r_{k} )
        BETA = (ALPHA / XI) * (SR2 / SR1)
        DO 90 I=1,NE
*         p_{k+1} = r_{k+1} + β_{k}( p_{k}-ξ_{k} A p_{k} )
          P(I) = R(I) + BETA * ( P(I) - XI*AP(I) )
   90   CONTINUE
   40 CONTINUE
* NITR まで計算しても収束せず
      WRITE (6,*) ' Not Converged ! '
* 収束と判定されたときの分岐点
  900 CONTINUE
* 結果出力
      DO 100 I = 1,NE
        WRITE(6,2000) I,X(I)
 2000   FORMAT(' X_13_( ',I3,' ) = ',1PD13.5)
  100 CONTINUE
      RETURN
      END
