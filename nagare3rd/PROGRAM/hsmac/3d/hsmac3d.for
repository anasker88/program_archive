*********************************************************************
* ファイル名  ：HSMAC3D.FOR                                         *
* タイトル    ：HSMAC法による3次元熱流動解析プログラム              *
* 製作者      ：平野　博之                                          *
* 所属        ：岡山理科大学 工学部 応用化学科                      *
* 製作日      ：2003.12.25                                          *
* 言語        ：FORTRAN (FORTRAN77でも実行可能)                     *
*********************************************************************
*                                                                   *
*    本プログラムでは，対流項は非保存形で1次精度風上差分を，拡散項  *
*  は2次精度中心差分を用いて離散化している．さらに高精度の近似を行う*
*  場合は，適宜変更のこと．                                         *
*  格子分割数を変更するときはPARAMETER文中のNX0,NY0,NZ0をすべて変更.*
*                                                                   *
*  ●変数                                                           *
*         →                                                        *
*    速度 Ｖ=(U,V,W),   圧力 P,   温度 T                            *
*                                                                   *
*  ●基礎方程式について                                             *
*        →                                                         *
*    ▽・Ｖ ＝０                                                    *
*      →                                                           *
*    ＤＶ                    ２→                                   *
*    －－＝－▽Ｐ＋ＶＩＳ（▽　Ｖ）＋ＢＵＯ×Ｔ(0,1,0)              *
*    Ｄｔ          ~~~~~~            ~~~~~~                         *
*    ＤＴ            ２                                             *
*    －－＝ＡＬＰ（▽  Ｔ）                                         *
*    Ｄｔ  ~~~~~~                                                   *
*    方程式に応じて，"BIS","BUO","ALP"を定義して与える．            *
*                                                                   *
*  ●格子分割について (NX=3, NY=3 の例)                             *
*       仮想セル 左側境界                  右側境界 仮想セル        *
*             ↓     ↓                         ↓ ↓               *
*           +--------＋-------+--------+--------＋-------+          *
*           |P(0,NY+1)        |        |        ∥P(NX+1,NY+1)      *
* 仮想セル→|   ・   →  ・   |   ・   |   ・   →  ・   |←仮想セル*
*           |      U(0,NY+1)  |        |      U(NX,NY+1) |          *
* 上方境界→+===↑===＋=======+========+========＋==↑===+←上方境界*
* (IY=NY)   |V(0,NY) ∥       |        |        ∥V(NX+1,NY)        *
*           |   ・   ∥  ・   |   ・   |   ・   ∥  ・   |          *
*           |        ∥       |        |        ∥       |          *
*           +--------＋-------+--------+--------＋-------+          *
*           |        ∥       |        |        ∥       |          *
*           |   ・   ∥  ・   |   ・   |   ・   ∥  ・   |          *
*           |        ∥       |        |        ∥       |          *
*           +--------＋-------+--------+--------＋-------+          *
*           |        ∥       |        |        ∥       |          *
*           |   ・   ∥  ・   |   ・   |   ・   ∥  ・   |          *
*           |  V(0,0)∥       |        |        ∥V(NX+1,0)         *
* 下方境界→+===↑===＋=======+========+========＋==↑===+←下方境界*
* (IY=0)    | P(0,0) ∥       |        |        ∥P(NX+1,0)         *
* 仮想セル→|   ・   →  ・   |   ・   |   ・   →  ・   |←仮想セル*
*           |      U(0,0)     |        |      U(NX,0)    |          *
*           +--------＋-------+--------+--------＋-------+          *
*              ↑    ↑                         ↑    ↑            *
*         仮想セル  左側境界                  右側境界 仮想セル     *
*                   (IX=0)                    (IX=NX)               *
*                                                                   *
*         Y                                                         *
*         ↑                                                        *
*         ｜                                                        *
*      ↓ ｜                                                        *
*      ｇ ｜                                                        *
*         ＋－－→X                                                 *
*        /                                                          *
*       /                                                           *
*      /                                                            *
*     Z                                                             *
*                                                                   *
*                                                                   *
*  ●スタッガードメッシュについて                                   *
*                               ←    DX     →                     *
*    +------------+-------------+-------------+                     *
*    |            |             |             | ↑                  *
*    |            |    P(i,j+1) |             |                     *
*    |            |      ・     |             | DY                  *
*    |            |             |             |                     *
*    |            |    V(i,j)   |             | ↓                  *
*    +------------+------↑-----+-------------+                     *
*    |            |             |             |                     *
*    |   P(i-1,j) |    P(i,j)   |   P(i+1,j)  | → U 定義点         *
*    |     ・     →     ・     →     ・     | ↑ V 定義点         *
*    |         U(i-1,j)       U(i,j)          | ・ P,T 定義         *
*    |            |    V(i,j-1) |             | (注）               *
*    +------------+------↑-----+-------------+ プログラム中のTは， *
*    |            |             |             | 本文中ではΘとなって*
*    |            |    P(i,j-1) |             | いる．              *
*    |            |      ・     |             |                     *
*    |            |             |             |                     *
*    |            |             |             |                     *
*    +------------+-------------+-------------+                     *
*                                                                   *
*             Y                                                     *
*             +--------------------------------------- +            *
*            /|         V(IX,IY,IZ)                   /|            *
*           / |             ↑                       / |            *
*          /  |             ｜                      /  |            *
*         /   |             ｜ W(IX,IY,IZ-1)       /   |            *
*        /    |             ｜ /                  /    |            *
*       /     |             －/                  /     |            *
*      /      |              /                  /      |            *
*     /       |             －                 /       |            *
*    /        |       P,T(IX,IY,IZ)           /        |            *
*   +----------------------------------------+         |            *
*   |   |---> |            ○                |  |--->U(IX,IY,IZ)    *
*   | U(IX-1,IY,IZ)--------------------------|---------+ X          *
*   |        /0          /                   |        /             *
*   |       /           /   ↑               |       /              *
*   |      /           /    ｜               |      /               *
*   |     /           －    ｜               |     /                *
*   |    /   W(IX,IY,IZ)    ｜               |    /                 *
*   |   /                   －               |   /                  *
*   |  /                   V(IX,IY-1,IZ)     |  /                   *
*   | /                                      | /                    *
*   |/                                       |/                     *
*   +----------------------------------------+                      *
*  Z                                                                *
*                                                                   *
*                                                                   *
*  ●パラメーターファイル（ｉｎ３ｄ．ｍａｃ）について               *
*    プログラムを実行すると，”ｉｎ３ｄ．ｍａｃ”というパラメータ   *
*    ファイルを読みに行くので，あらかじめ作成しておく．             *
*                                                                   *
*  "in3d.mac"のリスト                                               *
* 1: U.NEW ....Ｕの計算結果の出力ファイル名                         *
* 2: V.NEW ....Ｖの計算結果の出力ファイル名                         *
* 3: W.NEW ....Ｗの計算結果の出力ファイル名                         *
* 4: P.NEW ....Ｐの計算結果の出力ファイル名                         *
* 5: T.NEW ....Ｔの計算結果の出力ファイル名                         *
* 6: U.OLD ....継続計算のＵの入力データ                             *
* 7: V.OLD ....継続計算のＶの入力データ                             *
* 8: W.OLD ....継続計算のＷの入力データ                             *
* 9: P.OLD ....継続計算のＰの入力データ                             *
*10: T.OLD ....継続計算のＴの入力データ                             *
*11: UVWT.NEW ....Tecplot用データ                                   *
*12: +============================+                                 *
*13: | ITYPE   1===>   isothermal |                                 *
*14: |         2===>nonisothermal |                                 *
*15: +============================+                                 *
*16: -----------------------------------------------------------    *
*17: ITYPE     ICYCLE   NITR    NCYCLE                              *
*18: 2         0        10000   10000                               *
*19: -----------------------------------------------------------    *
*20: EPSP      OMG                                                  *
*21: 1.0e-3    1.7e+0                                               *
*22: -----------------------------------------------------------    *
*23: DT        RE       PR       GR                                 *
*24: 1.0e-4    0.0e+0   7.1e-1   1.0e+5                             *
*25: -----------------------------------------------------------    *
*26: DLX       DLY      DLZ      IRELP   METHOD                     *
*27: 1.0e+0    1.0e+0   1.0e0    0       5                          *
*                                                                   *
*    ITYPE.....１３，１４行を参照                                   *
*    ICYCLE.....計算開始のサイクル数（時間t=ICYCLE*DT）             *
*    NITR.....圧力計算のための，１サイクルあたりの最大反復回数      *
*    NCYCLE.....計算終了サイクル数                                  *
*                             →                                    *
*    EPSP.....収束判定値（▽・Ｖ ≦EPSPを満足するまで，反復計算によ *
*             って，圧力場を計算する．）                            *
*    OMG.....圧力計算のための緩和係数 , DT.....時間刻み             *
*    RE.....レイノルズ数 , PR.....プラントル数 , GR.....グラスホフ数*
*    DLX.....解析領域の横幅(DLX=NX*DX)                              *
*    DLY.....解析領域の高さ(DLY=NY*DY)                              *
*    DLZ.....解析領域の高さ(DLZ=NZ*DZ)                              *
*            格子幅（等間隔）DX,DY,DZはプログラムの中で求める．     *
*    IRELP...圧力の基準値の設定(0:行わない; 1:行う)                 *
*    METHOD..SMAC法においてのみ有効                                 *
*    パラメータファイルの数値について，                             *
*    FORTRAN プログラム.....できれば倍精度実数で与える"1.0D0,1.0d0" *
*            Cプログラムと共用させて"1.0E0,1.0e0"としても問題はない *
*    C プログラム....."1.0E0,1.0e0"として与える                     *
*    TECPLOT用データを除いた入出力ファイルは書式なし形式で，使用    *
*    するコンパイラーに依存する．コンパイラーに依存しない形式にする *
*    には，容量は増えるが書式付き形式に変更すればよい．             *
*                                                                   *
*  ●圧力の相対性について                                           *
*    圧力の相対性設定の変数名:IRELP                                 *
*    0:圧力の基準を設けない．                                       *
*    1:圧力の基準を設ける．                                         *
*      本プログラムではPO(1,1,1)=0となるように設定してある．        *
*      サブルーチン PRESS を参照．                                  *
*                                                                   *
*  ●変数・配列の説明                                               *
*    ICYCLE........タイムステップのためのカウンタ                   *
*    ITR...........圧力計算のための反復回数のカウンタ               *
*    IX,IY,IZ......スカラー量の定義点のｘ，ｙ座標の添字             *
*    UO,VO,W0,TO......圧力の反復計算を行う前の値                    *
*    UN,VN,WN,TN......収束した新たな圧力を用いて計算された速度      *
*                                                                   *
*********************************************************************
      PROGRAM HSMAC3D
*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
      CHARACTER FNAME(11)*20
*x方向の格子分割数
      NX  = NX0
*y方向の格子分割数
      NY  = NY0
*Z方向の格子分割数
      NZ  = NZ0
*パラメータファイルのオープン
      OPEN (10,FILE='IN3D.MAC',STATUS='OLD')
*出力ファイル名の読み込み
      DO 10 I = 1,11
        READ (10,'(A20)') FNAME(I)
   10 CONTINUE
*Uの計算結果出力用ファイルオープン(書式なし形式)
*書式なし形式はコンパイラーに依存するので注意
      OPEN (11, FILE=FNAME(1), STATUS='NEW', FORM='UNFORMATTED')
*Vの計算結果出力用ファイルオープン(書式なし形式)
      OPEN (12, FILE=FNAME(2), STATUS='NEW', FORM='UNFORMATTED')
*Wの計算結果出力用ファイルオープン(書式なし形式)
      OPEN (13, FILE=FNAME(3), STATUS='NEW', FORM='UNFORMATTED')
*Pの計算結果出力用ファイルオープン(書式なし形式)
      OPEN (14, FILE=FNAME(4), STATUS='NEW', FORM='UNFORMATTED')
*Tの計算結果出力用ファイルオープン(書式なし形式)
      OPEN (15, FILE=FNAME(5), STATUS='NEW', FORM='UNFORMATTED')
*in3d.mac中のコメント行(12-17行目)のスキップ
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) ITYPE, ICYCLE, NITR, NCYCLE
*in3d.mac中のコメント行(19-20行目)のスキップ
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) EPSP, OMG
*in3d.mac中のコメント行(22-23行目)のスキップ
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) DT, RE, PR, GR
*in3d.mac中のコメント行(25-26行目)のスキップ
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) DLX, DLY, DLZ, IRELP, METHOD
*継続の計算の場合
      IF (ICYCLE.NE.0) THEN
*Uデータファイルのオープン(書式なし形式)
        OPEN (16, FILE=FNAME(6), STATUS='OLD', FORM='UNFORMATTED')
*Vデータファイルのオープン(書式なし形式)
        OPEN (17, FILE=FNAME(7), STATUS='OLD', FORM='UNFORMATTED')
*Wデータファイルのオープン(書式なし形式)
        OPEN (18, FILE=FNAME(8), STATUS='OLD', FORM='UNFORMATTED')
*Pデータファイルのオープン(書式なし形式)
        OPEN (19, FILE=FNAME(9), STATUS='OLD', FORM='UNFORMATTED')
*Tデータファイルのオープン(等温場でもT=0.0のデータを読み込む)(書式なし形式)
        OPEN (20, FILE=FNAME(10), STATUS='OLD', FORM='UNFORMATTED')
      END IF
*x方向の格子幅
      DX  = DLX / FLOAT(NX)
*y方向の格子幅
      DY  = DLY / FLOAT(NY)
*Z方向の格子幅
      DZ  = DLZ / FLOAT(NZ)
*運動方程式中の拡散項の係数(ここではPr)
      VIS = PR
*エネルギー方程式中の拡散項の係数(ここでは1)
      ALP = 1.0D+0
*浮力項の係数(ここでは Gr * Pr**2)
      BUO = GR * PR**2
*等温場なら浮力項の係数はゼロ
      IF (ITYPE.EQ.1) BUO = 0.0D0
*初期値の設定
      CALL CINITI
*時間進行のための戻り点
  700 CONTINUE
*時間進行
      CALL ADV
*速度場の計算
      CALL CALVEL
*圧力計算の反復回数を1に初期化
      ITR = 1
*圧力反復のための戻り点
  710 CONTINUE
*圧力計算が収束したかどうかのパラメータIFLGを初期化
*IFLG -> 0:収束 1:発散(設定された許容回数NITR以下で解が得られない)
      IFLG = 0
*圧力場の計算
      CALL PRESS
*     Newton法による圧力場の計算が収束したとき
      IF ( IFLG. EQ. 0 ) THEN
*       非等温場計算の場合
        IF ( ITYPE.EQ.2 ) THEN
*         温度場を計算
          CALL CALTEM
        END IF
*     圧力場の計算が収束していないとき
      ELSE IF ( IFLG. EQ. 1 ) THEN
*       圧力計算の反復回数があらかじめ設定された最大値NITRより小さいとき
        IF ( ITR. LT. NITR ) THEN
*         さらに反復を繰り返す
          ITR = ITR + 1
          GO TO 710
*       圧力計算の反復回数がNITR以上のとき発散とみなして計算終了
        ELSE
          WRITE (6,*) ' NOT CONVERGE ! '
C         データを出力して強制終了
          CALL PROUT
          GO TO 900
        END IF
      END IF
*     時間進行カウンタ(ICYCLE)がNCYCLEより小さい時
      IF ( ICYCLE. LT. NCYCLE ) THEN
        GO TO 700
*     時間進行カウンタがNCYCLEになったら計算終了
      ELSE
        CALL PROUT
      END IF
  900 CONTINUE
*Tecplot用データの出力
      CALL TECPLT (FNAME(11))
      CLOSE (10)
      CLOSE (11)
      CLOSE (12)
      CLOSE (13)
      CLOSE (14)
      CLOSE (15)
      STOP
      END
*********************************************************************
*                          初期設定
*********************************************************************
      SUBROUTINE CINITI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
*新規計算の場合
      IF ( ICYCLE. EQ. 0 ) THEN
*       Uの初期値設定
        DO 10 IX = 0,NX
          DO 20 IY = 0,NY+1
            DO 30 IZ = 0,NZ+1
              UN(IX,IY,IZ) = 0.0D0
   30       CONTINUE
   20     CONTINUE
   10   CONTINUE
*       Vの初期値設定
        DO 40 IX = 0,NX+1
          DO 50 IY = 0,NY
            DO 60 IZ = 0,NZ+1
              VN(IX,IY,IZ) = 0.0D0
   60     CONTINUE
   50     CONTINUE
   40   CONTINUE
*       Wの初期値設定
        DO 70 IX = 0,NX+1
          DO 80 IY = 0,NY+1
            DO 90 IZ = 0,NZ
              WN(IX,IY,IZ) = 0.0D0
   90       CONTINUE
   80     CONTINUE
   70   CONTINUE
*       Pの初期値設定
        DO 100 IX = 0,NX+1
          DO 110 IY = 0,NY+1
            DO 120 IZ = 0,NZ+1
              PO(IX,IY,IZ) = 0.0D0
  120       CONTINUE
  110     CONTINUE
  100   CONTINUE
*-----------------------------------------------------------------------*
* （注意）浮力項の計算で温度の配列を使用しているので等温場でもT=0として *
* 初期条件だけは設定する必要がある．ゼロ以外の値を入れると浮力項が計算  *
* される可能性があるので注意．                                          *
*-----------------------------------------------------------------------*
*       Tの初期値設定(領域内は高温(+0.5)と低温(-0.5)の中間温度)
        DO 130 IX = 0,NX+1
          DO 140 IY = 0,NY+1
            DO 150 IZ = 0,NZ+1
              TN(IX,IY,IZ) = 0.0D0
  150     CONTINUE
  140     CONTINUE
  130   CONTINUE
*       Tの境界：右側壁（冷却）T=-0.5
        DO 160 IY = 0,NY+1
          DO 170 IZ = 0,NZ+1
            TN(NX+1,IY,IZ) = 2.0D0 * ( -0.5D0 ) - TN(NX,IY,IZ)
  170     CONTINUE
  160   CONTINUE
*       Tの境界：左側壁（加熱）T=+0.5
        DO 180 IY = 0,NY+1
          DO 190 IZ = 0,NZ+1
            TN(0,IY,IZ) = 2.0D0 * ( +0.5D0 ) - TN(1,IY,IZ)
  190   CONTINUE
  180   CONTINUE
*       Tの境界：上面（断熱）
        DO 200 IX = 1,NX
          DO 210 IZ = 0,NZ+1
            TN(IX,NY+1,IZ) = TN(IX,NY,IZ)
  210     CONTINUE
  200   CONTINUE
*       Tの境界：下面（断熱）
        DO 220 IX = 1,NX
          DO 230 IZ = 0,NZ+1
            TN(IX,0,IZ) = TN(IX,1,IZ)
  230     CONTINUE
  220   CONTINUE
*       Tの境界：後面（断熱）
        DO 240 IX = 1,NX
          DO 250 IY = 1,NY
            TN(IX,IY,NZ+1) = TN(IX,IY,NZ)
  250     CONTINUE
  240   CONTINUE
*       Tの境界：前面（断熱）
        DO 260 IX = 1,NX
          DO 270 IY = 1,NY
            TN(IX,IY,0) = TN(IX,IY,1)
  270     CONTINUE
  260   CONTINUE
*継続計算（すでにある計算結果からスタート）の場合
      ELSE
*       Uデータファイルからの読み込み[Unit No.=16](書式なし形式)
        READ (16) UN
*       Vデータファイルからの読み込み[Unit No.=17](書式なし形式)
        READ (17) VN
*       Wデータファイルからの読み込み[Unit No.=18](書式なし形式)
        READ (18) WN
*       Pデータファイルからの読み込み[Unit No.=19](書式なし形式)
        READ (19) PO
*---------------------------------------------------------------------*
*    (注意) 等温場の計算でもT(=0)のファイルを読み込む必要がある       *
*---------------------------------------------------------------------*
*       Tデータファイルからの読み込み[Unit No.=20](書式なし形式)
        READ (20) TN
        CLOSE (16)
        CLOSE (17)
        CLOSE (18)
        CLOSE (19)
        CLOSE (20)
      END IF
      RETURN
      END
*********************************************************************
*                         時間進行
*********************************************************************
      SUBROUTINE ADV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
      TIME   = DT*FLOAT(ICYCLE)
      ICYCLE = ICYCLE + 1
*ICYCLEを100回毎に表示
      IF (MOD(ICYCLE,100).EQ.0) THEN
        WRITE (6,2000) ICYCLE
 2000   FORMAT ('  CYC = ',I8)
      END IF
*--------------------------------------------------------------------
* UN -> UO : 必要なら入れ替える前にUNとUOから変動量を求める
* UN : 前の時間ステップにおいて最終的に得られた値，圧力補正の度に更新される
* UO : 新しい時間ステップでの初期値．UNを保存．
*--------------------------------------------------------------------
      DO 130 IX = 0,NX
        DO 140 IY = 0,NY+1
          DO 150 IZ = 0,NZ+1
            UO(IX,IY,IZ) = UN(IX,IY,IZ)
  150   CONTINUE
  140   CONTINUE
  130 CONTINUE
*--------------------------------------------------------------------
* VN -> VO : 必要なら入れ替える前にVNとVOから変動量を求める
* VN : 前の時間ステップにおいて最終的に得られた値，圧力補正の度に更新される
* VO : 新しい時間ステップでの初期値．VNを保存．
*--------------------------------------------------------------------
      DO 160 IX = 0,NX+1
        DO 170 IY = 0,NY
          DO 180 IZ = 0,NZ+1
            VO(IX,IY,IZ) = VN(IX,IY,IZ)
  180     CONTINUE
  170   CONTINUE
  160 CONTINUE
*--------------------------------------------------------------------
* WN -> WO : 必要なら入れ替える前にWNとWOから変動量を求める
* WN : 前の時間ステップにおいて最終的に得られた値，圧力補正の度に更新される
* WO : 新しい時間ステップでの初期値．WNを保存．
*--------------------------------------------------------------------
      DO 190 IX = 0,NX+1
        DO 200 IY = 0,NY+1
          DO 210 IZ = 0,NZ
            WO(IX,IY,IZ) = WN(IX,IY,IZ)
  210     CONTINUE
  200   CONTINUE
  190 CONTINUE
*--------------------------------------------------------------------
* TN -> TO : 必要なら入れ替える前にTNとTOから変動量を求める
* TN : 前の時間ステップでの計算値
* TO : 新しい時間ステップでの初期値．TNを保存．
*--------------------------------------------------------------------
      DO 220 IX = 0,NX+1
        DO 230 IY = 0,NY+1
          DO 240 IZ = 0,NZ+1
            TO(IX,IY,IZ) = TN(IX,IY,IZ)
  240     CONTINUE
  230   CONTINUE
  220 CONTINUE
      RETURN
      END
*********************************************************************
*                     速度場の計算
*********************************************************************
      SUBROUTINE CALVEL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
*--------------------------------------------------------------------
*                       U(IX,IY,IZ)の計算
*--------------------------------------------------------------------
      DO 10 IX = 1,NX-1
        DO 20 IY = 1,NY
          DO 30 IZ = 1,NZ
*       VVはU(IX,IY,IZ)におけるVの補間値
        VV = (  VO(IX,IY  ,IZ)+VO(IX+1,IY  ,IZ)
     $         +VO(IX,IY-1,IZ)+VO(IX+1,IY-1,IZ) )/4.0D0
*       WWはU(IX,IY,IZ)におけるWの補間値
        WW = (  WO(IX,IY,IZ  )+WO(IX+1,IY,IZ  )
     $         +WO(IX,IY,IZ-1)+WO(IX+1,IY,IZ-1) )/4.0D0
*       対流項(CNVUX,CNVUY,CNVUZ)を１次精度風上差分にて計算
        IF ( UO(IX,IY,IZ). GE. 0.0D0 ) THEN
          CNVUX = UO(IX,IY,IZ)*( UO(IX,IY,IZ) - UO(IX-1,IY,IZ) ) / DX
        ELSE IF ( UO(IX,IY,IZ). LT. 0.0D0 ) THEN
          CNVUX = UO(IX,IY,IZ)*( UO(IX+1,IY,IZ) - UO(IX,IY,IZ) ) / DX
        END IF
        IF ( VV. GE. 0.0D0 ) THEN
          CNVUY = VV*( UO(IX,IY,IZ) - UO(IX,IY-1,IZ) ) / DY
        ELSE IF ( VV. LT. 0.0D0 ) THEN
          CNVUY = VV*( UO(IX,IY+1,IZ) - UO(IX,IY,IZ) ) / DY
        END IF
        IF ( WW. GE. 0.0D0 ) THEN
          CNVUZ = WW*( UO(IX,IY,IZ) - UO(IX,IY,IZ-1) ) / DZ
        ELSE IF ( WW. LT. 0.0D0 ) THEN
          CNVUZ = WW*( UO(IX,IY,IZ+1) - UO(IX,IY,IZ) ) / DZ
        END IF
*       x方向の浮力項(BUOU)はゼロ
        TU = 0.0D0
        BUOU = BUO * TU
*       拡散項(DIFU)の計算
        DIFU = VIS*(
     $  +( UO(IX-1,IY,IZ)-2.0D0*UO(IX,IY,IZ)+UO(IX+1,IY,IZ) )/DX**2
     $  +( UO(IX,IY-1,IZ)-2.0D0*UO(IX,IY,IZ)+UO(IX,IY+1,IZ) )/DY**2
     $  +( UO(IX,IY,IZ-1)-2.0D0*UO(IX,IY,IZ)+UO(IX,IY,IZ+1) )/DZ**2
     $              )
*       仮の速度(U)の計算
        UN(IX,IY,IZ) = UO(IX,IY,IZ)
     $  + DT*( -CNVUX-CNVUY-CNVUZ+DIFU+BUOU
     $         +(PO(IX,IY,IZ)-PO(IX+1,IY,IZ))/DX )
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
*--------------------------------------------------------------------
*                    V(IX,IY,IZ)の計算
*--------------------------------------------------------------------
      DO 40 IX = 1,NX
        DO 50 IY = 1,NY-1
          DO 60 IZ = 1,NZ
*       UUはV(IX,IY,IZ)におけるUの補間値
        UU = (  UO(IX-1,IY  ,IZ)+UO(IX,IY  ,IZ)
     $         +UO(IX-1,IY+1,IZ)+UO(IX,IY+1,IZ) )/4.0D0
*       WWはV(IX,IY,IZ)におけるWの補間値
        WW = (  WO(IX,IY  ,IZ-1)+WO(IX,IY  ,IZ)
     $         +WO(IX,IY+1,IZ-1)+WO(IX,IY+1,IZ) )/4.0D0
*       対流項(CNVVX,CNVVY,CNVVZ)を１次精度風上差分にて計算
        IF ( UU. GE. 0.0D0 ) THEN
          CNVVX = UU*( VO(IX,IY,IZ) - VO(IX-1,IY,IZ) ) / DX
        ELSE IF ( UU. LT. 0.0D0 ) THEN
          CNVVX = UU*( VO(IX+1,IY,IZ) - VO(IX,IY,IZ) ) / DX
        END IF
        IF ( VO(IX,IY,IZ). GE. 0.0D0 ) THEN
          CNVVY = VO(IX,IY,IZ)*( VO(IX,IY,IZ) - VO(IX,IY-1,IZ) ) / DY
        ELSE IF ( VO(IX,IY,IZ). LT. 0.0D0 ) THEN
          CNVVY = VO(IX,IY,IZ)*( VO(IX,IY+1,IZ) - VO(IX,IY,IZ) ) / DY
        END IF
        IF ( WW. GE. 0.0D0 ) THEN
          CNVVZ = WO(IX,IY,IZ)*( VO(IX,IY,IZ) - VO(IX,IY,IZ-1) ) / DZ
        ELSE IF ( WW. LT. 0.0D0 ) THEN
          CNVVZ = WO(IX,IY,IZ)*( VO(IX,IY,IZ+1) - VO(IX,IY,IZ) ) / DZ
        END IF
*       浮力項(BUOV)の計算
        TV = ( TO(IX,IY,IZ) + TO(IX,IY+1,IZ) )/2.0D0
        BUOV = BUO*TV
*       拡散項(DIFV)の計算
        DIFV = VIS*(
     $  +(VO(IX-1,IY,IZ)-2.0D0*VO(IX,IY,IZ)+VO(IX+1,IY,IZ))/DX**2
     $  +(VO(IX,IY-1,IZ)-2.0D0*VO(IX,IY,IZ)+VO(IX,IY+1,IZ))/DY**2
     $  +(VO(IX,IY,IZ-1)-2.0D0*VO(IX,IY,IZ)+VO(IX,IY,IZ+1))/DZ**2
     $             )
*       仮の速度(V)の計算
        VN(IX,IY,IZ) = VO(IX,IY,IZ)
     $  + DT*( -CNVVX-CNVVY-CNVVZ+DIFV+BUOV
     $         +(PO(IX,IY,IZ)-PO(IX,IY+1,IZ))/DY )
   60     CONTINUE
   50   CONTINUE
   40 CONTINUE
*--------------------------------------------------------------------
*                     W(IX,IY,IZ)の計算
*--------------------------------------------------------------------
      DO 70 IX = 1,NX
        DO 80 IY = 1,NY
          DO 90 IZ = 1,NZ-1
*       UUはW(IX,IY,IZ)におけるUの補間値
        UU = (  UO(IX-1,IY,IZ  )+UO(IX,IY,IZ)
     $         +UO(IX-1,IY,IZ+1)+UO(IX,IY,IZ+1) )/4.0D0
*       VVはW(IX,IY,IZ)におけるVの補間値
        VV = (  VO(IX,IY-1,IZ  )+VO(IX,IY,IZ)
     $         +VO(IX,IY-1,IZ+1)+VO(IX,IY,IZ+1) )/4.0D0
*       対流項(CNVWX,CNVWY,CNVWZ)を１次精度風上差分にて計算
        IF ( UU. GE. 0.0D0 ) THEN
          CNVWX = UU*( WO(IX,IY,IZ) - WO(IX-1,IY,IZ) ) / DX
        ELSE IF ( UU. LT. 0.0D0 ) THEN
          CNVWX = UU*( WO(IX+1,IY,IZ) - WO(IX,IY,IZ) ) / DX
        END IF
        IF ( VV. GE. 0.0D0 ) THEN
          CNVWY = VV*( WO(IX,IY,IZ) - WO(IX,IY-1,IZ) ) / DY
        ELSE IF ( VV. LT. 0.0D0 ) THEN
          CNVWY = VV*( WO(IX,IY+1,IZ) - WO(IX,IY,IZ) ) / DY
        END IF
        IF ( WO(IX,IY,IZ). GE. 0.0D0 ) THEN
          CNVWZ = WO(IX,IY,IZ)*( WO(IX,IY,IZ) - WO(IX,IY,IZ-1) ) / DZ
        ELSE IF ( WO(IX,IY,IZ). LT. 0.0D0 ) THEN
          CNVWZ = WO(IX,IY,IZ)*( WO(IX,IY,IZ+1) - WO(IX,IY,IZ) ) / DZ
        END IF
*       浮力項(BUOW)の計算
        TW = 0.0D0
        BUOW = BUO * TW
*       拡散項(DIFW)の計算
        DIFW = VIS*(
     $  +(WO(IX-1,IY,IZ)-2.0D0*WO(IX,IY,IZ)+WO(IX+1,IY,IZ))/DX**2
     $  +(WO(IX,IY-1,IZ)-2.0D0*WO(IX,IY,IZ)+WO(IX,IY+1,IZ))/DY**2
     $  +(WO(IX,IY,IZ-1)-2.0D0*WO(IX,IY,IZ)+WO(IX,IY,IZ+1))/DZ**2
     $             )
*       仮の速度(W)の計算
        WN(IX,IY,IZ) = WO(IX,IY,IZ)
     $  + DT*( -CNVWX-CNVWY-CNVWZ+DIFW+BUOW
     $         +(PO(IX,IY,IZ)-PO(IX,IY,IZ+1))/DZ )
   90     CONTINUE
   80   CONTINUE
   70 CONTINUE
*速度の境界条件の処理
      CALL VELBND
      RETURN
      END
*********************************************************************
*                    圧力場の計算
*********************************************************************
      SUBROUTINE PRESS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
      IXMAX = 0
      IYMAX = 0
      IZMAX = 0
      DMAX = 0.0D0
*P(IX,IY,IZ)の計算
      DO 10 IX = 1,NX
        DO 20 IY = 1,NY
          DO 30 IZ = 1,NZ
          DEL = DT*( 2.0D0/DX**2 + 2.0D0/DY**2 + 2.0D0/DZ**2 )
          DIV = ( UN(IX,IY,IZ) - UN(IX-1,IY  ,IZ  ) )/DX
     $        + ( VN(IX,IY,IZ) - VN(IX  ,IY-1,IZ  ) )/DY
     $        + ( WN(IX,IY,IZ) - WN(IX  ,IY  ,IZ-1) )/DZ
          IF ( DABS(DIV). GE. DABS(DMAX) ) THEN
            IXMAX = IX
            IYMAX = IY
            IZMAX = IZ
            DMAX = DIV
          END IF
          DELP = - OMG * DIV / DEL
          PO(IX,  IY  ,IZ  ) = PO(IX,  IY,  IZ  ) + DELP
          UN(IX,  IY  ,IZ  ) = UN(IX,  IY,  IZ  ) + DT/DX*DELP
          UN(IX-1,IY  ,IZ  ) = UN(IX-1,IY,  IZ  ) - DT/DX*DELP
          VN(IX,  IY  ,IZ  ) = VN(IX,  IY,  IZ  ) + DT/DY*DELP
          VN(IX,  IY-1,IZ  ) = VN(IX,  IY-1,IZ  ) - DT/DY*DELP
          WN(IX,  IY  ,IZ  ) = WN(IX,  IY,  IZ  ) + DT/DZ*DELP
          WN(IX,  IY  ,IZ-1) = WN(IX,  IY,  IZ-1) - DT/DZ*DELP
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
*
* 圧力の相対性に関する処理(IRELP=1なら以下の処理を行う)
      IF (IRELP.EQ.1) THEN
        POSTN = PO(1,1,1)
        DO 40 IX = 1,NX
          DO 50 IY = 1,NY
            DO 60 IZ = 1,NZ
              PO(IX,IY,IZ) = PO(IX,IY,IZ)-POSTN
   60       CONTINUE
   50     CONTINUE
   40   CONTINUE
      END IF
*
* IFLG=1なら連続の式を満たしていないと判定し再び圧力計算を行う．
      IF ( DABS(DMAX). GE. EPSP ) IFLG = 1
* 圧力計算の回数を100回ごとに表示
      IF (MOD(ITR,100).EQ.0) THEN
        WRITE (6,2000) ITR, IXMAX, IYMAX, IZMAX, DMAX
 2000   FORMAT (' Iteration=',I8,'   Div(max)(',3I6,')=',1PE13.5)
      END IF
*新たに得られた速度を用いて境界条件を処理する
      CALL VELBND
      RETURN
      END
*********************************************************************
*                     温度場の計算
*********************************************************************
      SUBROUTINE CALTEM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
*T(IX,IY,IZ)の計算
      DO 10 IX = 1,NX
        DO 20 IY = 1,NY
          DO 30 IZ = 1,NZ
*         UUT,VVT,WWTはそれぞれT(IX,IY,IZ)におけるU,Vの補間値
          UUT = ( UO(IX,IY,IZ) + UO(IX-1,IY,  IZ  ) ) / 2.0D0
          VVT = ( VO(IX,IY,IZ) + VO(IX  ,IY-1,IZ  ) ) / 2.0D0
          WWT = ( WO(IX,IY,IZ) + WO(IX  ,IY,  IZ-1) ) / 2.0D0
*         対流項(CNVTX,CNVTY,CNVTZ)を１次精度風上差分にて計算
          IF ( UUT. GE. 0.0D0 ) THEN
            CNVTX = UUT*( TO(IX,IY,IZ) - TO(IX-1,IY,IZ) ) / DX
          ELSE IF ( UUT. LT. 0.0D0 ) THEN
            CNVTX = UUT*( TO(IX+1,IY,IZ) - TO(IX,IY,IZ) ) / DX
          END IF
          IF ( VVT. GE. 0.0D0 ) THEN
            CNVTY = VVT*( TO(IX,IY,IZ) - TO(IX,IY-1,IZ) ) / DY
          ELSE IF ( VVT. LT. 0.0D0 ) THEN
            CNVTY = VVT*( TO(IX,IY+1,IZ) - TO(IX,IY,IZ) ) / DY
          END IF
          IF ( WWT. GE. 0.0D0 ) THEN
            CNVTZ = WWT*( TO(IX,IY,IZ) - TO(IX,IY,IZ-1) ) / DZ
          ELSE IF ( WWT. LT. 0.0D0 ) THEN
            CNVTZ = WWT*( TO(IX,IY,IZ+1) - TO(IX,IY,IZ) ) / DZ
          END IF
*         拡散項(DIFT)の計算
          DIFT = ALP*(
     $    +( TO(IX-1,IY,IZ)-2.0D0*TO(IX,IY,IZ)+TO(IX+1,IY,IZ) )/DX**2
     $    +( TO(IX,IY-1,IZ)-2.0D0*TO(IX,IY,IZ)+TO(IX,IY+1,IZ) )/DY**2
     $    +( TO(IX,IY,IZ-1)-2.0D0*TO(IX,IY,IZ)+TO(IX,IY,IZ+1) )/DZ**2
     $                )
*         次の時間のTの計算
          TN(IX,IY,IZ) = TO(IX,IY,IZ) + DT*( -CNVTX-CNVTY-CNVTZ+DIFT )
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
*境界条件の処理
      CALL TBND
      RETURN
      END
*********************************************************************
*                 速度の境界条件の処理
*********************************************************************
      SUBROUTINE VELBND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
*U（右面）
      DO 10 IY = 1,NY
        DO 20 IZ = 1,NZ
          UN(NX,IY,IZ) = 0.0D0
   20   CONTINUE
   10 CONTINUE
*U（左面）
      DO 30 IY = 1,NY
        DO 40 IZ = 1,NZ
          UN(0,IY,IZ) = 0.0D0
   40   CONTINUE
   30 CONTINUE
*U（後面）
      DO 50 IX = 0,NX
        DO 60 IY = 1,NY
          UN(IX,IY,NZ+1) = -UN(IX,IY,NZ)
   60   CONTINUE
   50 CONTINUE
*U（前面）
      DO 70 IX = 0,NX
        DO 80 IY = 1,NY
          UN(IX,IY,0) = -UN(IX,IY,1)
   80   CONTINUE
   70 CONTINUE
*U（上面）
      DO 90 IX = 0,NX
        DO 100 IZ = 0,NZ+1
          UN(IX,NY+1,IZ) = -UN(IX,NY,IZ)
  100   CONTINUE
   90 CONTINUE
*U（下面）
      DO 110 IX = 0,NX
        DO 120 IZ = 0,NZ+1
          UN(IX,0,IZ) = -UN(IX,1,IZ)
  120   CONTINUE
  110 CONTINUE
*V（右面）
      DO 200 IY = 1,NY-1
        DO 210 IZ = 1,NZ
          VN(NX+1,IY,IZ) = -VN(NX,IY,IZ)
  210   CONTINUE
  200 CONTINUE
*V（左面）
      DO 220 IY = 1,NY-1
        DO 230 IZ = 1,NZ
          VN(0,IY,IZ) = -VN(1,IY,IZ)
  230   CONTINUE
  220 CONTINUE
*V（上面）
      DO 240 IX = 0,NX+1
        DO 250 IZ = 0,NZ+1
          VN(IX,NY,IZ) = 0.0D0
  250   CONTINUE
  240 CONTINUE
*V（下面）
      DO 260 IX = 0,NX+1
        DO 270 IZ = 0,NZ+1
          VN(IX,0,IZ) = 0.0D0
  270   CONTINUE
  260 CONTINUE
*V（後面）
      DO 280 IX = 0,NX+1
        DO 290 IY = 1,NY-1
          VN(IX,IY,NZ+1) = -VN(IX,IY,NZ)
  290   CONTINUE
  280 CONTINUE
*V（前面）
      DO 300 IX = 0,NX+1
        DO 310 IY = 1,NY-1
          VN(IX,IY,0) = -VN(IX,IY,1)
  310   CONTINUE
  300 CONTINUE
*W（右面）
      DO 410 IY = 1,NY
        DO 420 IZ = 1,NZ-1
          WN(NX+1,IY,IZ) = -WN(NX,IY,IZ)
  420   CONTINUE
  410 CONTINUE
*W（左面）
      DO 430 IY = 1,NY
        DO 440 IZ = 1,NZ-1
          WN(0,IY,IZ) = -WN(1,IY,IZ)
  440   CONTINUE
  430 CONTINUE
*W（後面）
      DO 450 IX = 0,NX
        DO 460 IY = 0,NY+1
          WN(IX,IY,NZ) = 0.0D0
  460   CONTINUE
  450 CONTINUE
*W（前面）
      DO 470 IX = 0,NX
        DO 480 IY = 0,NY+1
          WN(IX,IY,0) = 0.0D0
  480   CONTINUE
  470 CONTINUE
*W（上面）
      DO 490 IX = 0,NX+1
        DO 500 IZ = 1,NZ-1
          WN(IX,NY+1,IZ) = -WN(IX,NY,IZ)
  500   CONTINUE
  490 CONTINUE
*W（下面）
      DO 510 IX = 0,NX+1
        DO 520 IZ = 1,NZ-1
          WN(IX,0,IZ) = -WN(IX,1,IZ)
  520   CONTINUE
  510 CONTINUE
      RETURN
      END
*********************************************************************
*                 温度の境界条件の処理
*********************************************************************
      SUBROUTINE TBND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
*右面
      DO 10 IY = 0,NY+1
        DO 20 IZ = 0,NZ+1
          TN(NX+1,IY,IZ) = 2.0D0 * ( -0.5D0 ) - TN(NX,IY,IZ)
   20   CONTINUE
   10 CONTINUE
*左面
      DO 30 IY = 0,NY+1
        DO 40 IZ = 0,NZ+1
          TN(0,IY,IZ) = 2.0D0 * ( +0.5D0 ) - TN(1,IY,IZ)
   40   CONTINUE
   30 CONTINUE
*上面
      DO 50 IX = 1,NX
        DO 60 IZ = 1,NZ
          TN(IX,NY+1,IZ) = TN(IX,NY,IZ)
   60   CONTINUE
   50 CONTINUE
*下面
      DO 70 IX = 1,NX
        DO 80 IZ = 1,NZ
          TN(IX,0,IZ) = TN(IX,1,IZ)
   80   CONTINUE
   70 CONTINUE
*後面
      DO 90 IX = 1,NX
        DO 100 IY = 0,NY+1
          TN(IX,IY,NZ+1) = TN(IX,IY,NZ)
  100   CONTINUE
   90 CONTINUE
*前面
      DO 110 IX = 1,NX
        DO 120 IY = 0,NY+1
          TN(IX,IY,0) = TN(IX,IY,1)
  120   CONTINUE
  110 CONTINUE
      RETURN
      END
*********************************************************************
*                        データ出力
*********************************************************************
      SUBROUTINE PROUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
       WRITE (11) UN
       WRITE (12) VN
       WRITE (13) WN
       WRITE (14) PO
       WRITE (15) TN
      RETURN
      END
*********************************************************************
*                  Tecplot用データ出力
*********************************************************************
      SUBROUTINE TECPLT (FNAME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
       CHARACTER*20 FNAME
*
      OPEN (21,FILE=FNAME,STATUS='NEW')
      WRITE (21,*) 'VARIABLES = "X", "Y", "Z", "U", "V", "W", "T"'
      NX1 = NX+1
      NY1 = NY+1
      NZ1 = NZ+1
      WRITE (21,4000) NX1,NY1,NZ1
 4000 FORMAT (1H ,'ZONE I=',I3,',J=',I3,',K=',I3,', F=POINT')
      DO 10 IZ = 0,NZ
         DO 20 IY = 0,NY
           DO 30 IX = 0,NX
             X = DX * FLOAT(IX)
             Y = DY * FLOAT(IY)
             Z = DZ * FLOAT(IZ)
             U = ( UN(IX  ,IY,  IZ  )+UN(IX  ,IY+1,IZ  )
     $            +UN(IX  ,IY,  IZ+1)+UN(IX  ,IY+1,IZ+1) )/4.0D0
             V = ( VN(IX  ,IY,  IZ  )+VN(IX+1,IY  ,IZ  )
     $            +VN(IX  ,IY,  IZ+1)+VN(IX+1,IY,  IZ+1) )/4.0D0
             W = ( WN(IX  ,IY  ,IZ  )+WN(IX  ,IY+1,IZ  )
     $            +WN(IX+1,IY  ,IZ  )+WN(IX+1,IY+1,IZ  ) )/4.0D0
             T = ( TN(IX  ,IY  ,IZ  )+TN(IX+1,IY  ,IZ  )
     $            +TN(IX  ,IY+1,IZ  )+TN(IX+1,IY+1,IZ  )
     $            +TN(IX  ,IY  ,IZ+1)+TN(IX+1,IY  ,IZ+1)
     $            +TN(IX  ,IY+1,IZ+1)+TN(IX+1,IY+1,IZ+1) )/8.0D0
            WRITE (21,4010) X,Y,Z,U,V,W,T
 4010       FORMAT (1H ,7(1PE11.3))
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
      CLOSE(21)
      RETURN
      END
