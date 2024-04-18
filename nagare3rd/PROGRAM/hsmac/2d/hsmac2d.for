*********************************************************************
* ファイル名  ：HSMAC2D.FOR                                         *
* タイトル    ：HSMAC法による2次元熱流動解析プログラム              *
* 製作者      ：平野　博之                                          *
* 所属        ：岡山理科大学 工学部 応用化学科                      *
* 製作日      ：2003.12.25                                          *
* 言語        ：FORTRAN (FORTRAN77でも実行可能)                     *
*********************************************************************
*                                                                   *
*    本プログラムでは，対流項は非保存形で1次精度風上差分を，拡散項  *
*  は2次精度中心差分を用いて離散化している．さらに高精度の近似を行う*
*  場合は，適宜変更のこと．                                         *
*  格子分割数を変更するときは，PARAMETER文中のNX0,NY0をすべて変更． *
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
*    －－＝－▽Ｐ＋ＶＩＳ（▽　Ｖ）＋ＢＵＯ×Ｔ(0,1)                *
*    Ｄｔ          ~~~~~~            ~~~~~~                         *
*    ＤＴ            ２                                             *
*    －－＝ＡＬＰ（▽  Ｔ）                                         *
*    Ｄｔ  ~~~~~~                                                   *
*    方程式に応じて，"BIS","BUO","ALP"を定義して与える．            *
*                                                                   *
*  ●格子分割について (NX=3, NY=3の例)                              *
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
*     ｙ    +--------＋-------+--------+--------＋-------+          *
*     ↑       ↑    ↑                         ↑    ↑            *
*  ↓ ｜  仮想セル  左側境界                  右側境界 仮想セル     *
*  ｇ ｜            (IX=0)                    (IX=NX)               *
*     ＋－→ｘ                                                      *
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
*    |         U(i-1,j)       U(i,j)          | ・ P,T 定義点       *
*    |            |    V(i,j-1) |             | (注）               *
*    +------------+------↑-----+-------------+ プログラム中のTは， *
*    |            |             |             | 本文中ではΘとなって*
*    |            |    P(i,j-1) |             | いる．              *
*    |            |      ・     |             |                     *
*    |            |             |             |                     *
*    |            |             |             |                     *
*    +------------+-------------+-------------+                     *
*                                                                   *
*  ●パラメーターファイル（ｉｎ２ｄ．ｍａｃ）について               *
*    プログラムを実行すると，”ｉｎ２ｄ．ｍａｃ”というパラメータ   *
*    ファイルを読みに行くので，あらかじめ作成しておく．             *
*                                                                   *
*  "in2d.mac"のリスト                                               *
* 1: U.NEW ....Ｕの計算結果の出力ファイル名                         *
* 2: V.NEW ....Ｖの計算結果の出力ファイル名                         *
* 3: P.NEW ....Ｐの計算結果の出力ファイル名                         *
* 4: T.NEW ....Ｔの計算結果の出力ファイル名                         *
* 5: U.OLD ....継続計算のＵの入力データ                             *
* 6: V.OLD ....継続計算のＶの入力データ                             *
* 7: P.OLD ....継続計算のＰの入力データ                             *
* 8: T.OLD ....継続計算のＴの入力データ                             *
* 9: UVT.NEW ....Tecplot用データ                                    *
*10: +============================+                                 *
*11: | ITYPE   1===>   isothermal |                                 *
*12: |         2===>nonisothermal |                                 *
*13: +============================+                                 *
*14: -----------------------------------------------------------    *
*15: ITYPE     ICYCLE   NITR    NCYCLE                              *
*16: 2         0        10000   1000                                *
*17: -----------------------------------------------------------    *
*18: EPSP      OMG                                                  *
*19: 1.0e-3    1.7e+0                                               *
*20: -----------------------------------------------------------    *
*21: DT        RE       PR       GR                                 *
*22: 1.0e-4    0.0e+0   7.1e-1   1.0e+5                             *
*23: -----------------------------------------------------------    *
*24: DLX       DLY      IRELP   METHOD                              *
*25: 1.0e+0    1.0e+0   0       5                                   *
*                                                                   *
*    ITYPE.....１１，１２行を参照                                   *
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
*            格子幅（等間隔）DX,DYはプログラムの中で求める．        *
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
*      本プログラムではPO(1,1)=0となるように設定してある．          *
*      サブルーチン PRESS を参照．                                  *
*                                                                   *
*  ●変数・配列の説明                                               *
*    ICYCLE -----> 時間進行のためのカウンタ                         *
*    ITR    -----> 圧力計算のための反復回数のカウンタ               *
*    IX,IY  -----> 上の図を参照                                     *
*    UO,VO,TO----->圧力の反復計算を行う前の値                       *
*    UN,VN,TN----->収束した新たな圧力を用いて計算された値           *
*                                                                   *
*********************************************************************
      PROGRAM HSMAC2D
*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
      CHARACTER FNAME(10)*20
*x方向の格子分割数
      NX  = NX0
*y方向の格子分割数
      NY  = NY0
*パラメータファイルのオープン
      OPEN (10,FILE='IN2D.MAC',STATUS='OLD')
*出力ファイル名の読み込み
      DO 10 I = 1,9
        READ (10,'(A20)') FNAME(I)
   10 CONTINUE
*Uの計算結果出力用ファイルオープン(書式なし形式)
*書式なし形式はコンパイラーに依存するので注意
      OPEN (11, FILE=FNAME(1), STATUS='NEW', FORM='UNFORMATTED')
*Vの計算結果出力用ファイルオープン(書式なし形式)
      OPEN (12, FILE=FNAME(2), STATUS='NEW', FORM='UNFORMATTED')
*Pの計算結果出力用ファイルオープン(書式なし形式)
      OPEN (13, FILE=FNAME(3), STATUS='NEW', FORM='UNFORMATTED')
*Tの計算結果出力用ファイルオープン(書式なし形式)
      OPEN (14, FILE=FNAME(4), STATUS='NEW', FORM='UNFORMATTED')
*in2d.mac中のコメント行(10-15行目)のスキップ
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) ITYPE, ICYCLE, NITR, NCYCLE
*in2d.mac中のコメント行(17-18行目)のスキップ
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) EPSP, OMG
*in2d.mac中のコメント行(20-21行目)のスキップ
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) DT,RE, PR, GR
*in2d.mac中のコメント行(23-24行目)のスキップ
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) DLX, DLY, IRELP, METHOD
*継続の計算の場合
      IF (ICYCLE.NE.0) THEN
*Uデータファイルのオープン(書式なし形式)
        OPEN (15, FILE=FNAME(5), STATUS='OLD', FORM='UNFORMATTED')
*Vデータファイルのオープン(書式なし形式)
        OPEN (16, FILE=FNAME(6), STATUS='OLD', FORM='UNFORMATTED')
*Pデータファイルのオープン(書式なし形式)
        OPEN (17, FILE=FNAME(7), STATUS='OLD', FORM='UNFORMATTED')
*Tデータファイルのオープン(等温場でもT=0.0のデータを読み込む)(書式なし形式)
        OPEN (18, FILE=FNAME(8), STATUS='OLD', FORM='UNFORMATTED')
      END IF
*x方向の格子幅
      DX  = DLX / FLOAT(NX)
*y方向の格子幅
      DY  = DLY / FLOAT(NY)
*運動方程式中の拡散項の係数(ここではPr)
      VIS = PR
*エネルギー方程式中の拡散項の係数(ここでは1)
      ALP = 1.0D+0
*浮力項の係数(ここでは Gr * Pr**2)
      BUO = GR * PR**2
*等温場なら浮力項の係数はゼロに設定
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
*       圧力計算の反復回数がNITRとなったら発散とみなして計算終了
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
      CALL TECPLT(FNAME(9))
      CLOSE (10)
      CLOSE (11)
      CLOSE (12)
      CLOSE (13)
      CLOSE (14)
      STOP
      END
*********************************************************************
*                        初期設定
*********************************************************************
      SUBROUTINE CINITI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
*新規計算の場合
      IF ( ICYCLE. EQ. 0 ) THEN
*       Uの初期値設定
        DO 10 IX = 0,NX
          DO 20 IY = 0,NY+1
            UN(IX,IY) = 0.0D0
   20     CONTINUE
   10   CONTINUE
*       Vの初期値設定
        DO 30 IX = 0,NX+1
          DO 40 IY = 0,NY
            VN(IX,IY) = 0.0D0
   40     CONTINUE
   30   CONTINUE
*       Pの初期値設定
        DO 50 IX = 0,NX+1
          DO 60 IY = 0,NY+1
            PO(IX,IY) = 0.0D0
   60     CONTINUE
   50   CONTINUE
*----------------------------------------------------------------------*
*（注意）浮力項の計算で温度の配列を使用しているので等温場でもT=0として *
* 初期条件だけは設定する必要がある．ゼロ以外の値を入れると浮力項が計算 *
* される可能性があるので注意．                                         *
*----------------------------------------------------------------------*
*       Tの初期値設定(領域内は高温(+0.5)と低温(-0.5)の中間温度)
        DO 61 IX = 0,NX+1
          DO 62 IY = 0,NY+1
            TN(IX,IY) = 0.0D0
   62     CONTINUE
   61   CONTINUE
*       Tの境界：右面（冷却）T=-0.5
        DO 70 IY = 0,NY+1
          TN(NX+1,IY) = 2.0D0 * ( -0.5D0 ) - TN(NX,IY)
   70   CONTINUE
*       Tの境界：左面（加熱）T=+0.5
        DO 80 IY = 0,NY+1
          TN(0,IY) = 2.0D0 * ( +0.5D0 ) - TN(1,IY)
   80   CONTINUE
*       Tの境界：上面（断熱）
        DO 90 IX = 1,NX
          TN(IX,NY+1) = TN(IX,NY)
   90   CONTINUE
*       Tの境界：下面（断熱）
        DO 95 IX = 1,NX
          TN(IX,0) = TN(IX,1)
   95   CONTINUE
*継続計算（すでにある計算結果からスタート）の場合
      ELSE
*       Uデータファイルからの読み込み[Unit No.=15](書式なし形式)
        READ (15) UN
*       Vデータファイルからの読み込み[Unit No.=16](書式なし形式)
        READ (16) VN
*       Pデータファイルからの読み込み[Unit No.=17](書式なし形式)
        READ (17) PO
*---------------------------------------------------------------------*
*    (注意) 等温場の計算でもT(=0)のファイルを読み込む必要がある       *
*---------------------------------------------------------------------*
*       Tデータファイルからの読み込み[Unit No.=18](書式なし形式)
        READ (18) TN
        CLOSE (15)
        CLOSE (16)
        CLOSE (17)
        CLOSE (18)
      END IF
      RETURN
      END
*********************************************************************
*                         時間進行
*********************************************************************
      SUBROUTINE ADV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
      TIME   = DT*FLOAT(ICYCLE)
      ICYCLE = ICYCLE + 1
*時間進行カウンタ(ICYCLE)を100回毎に表示
      IF (MOD(ICYCLE,100).EQ.0) THEN
        WRITE (6,2000) ICYCLE
 2000   FORMAT ('  CYC = ',I8)
      END IF
*--------------------------------------------------------------------
* UN -> UO : 必要なら入れ替える前にUNとUOから変動量を求める
* UN : 前の時間ステップにおいて最終的に得られた値，圧力補正の度に更新される
* UO : 新しい時間ステップでの初期値．UNを保存．
*--------------------------------------------------------------------
      DO 70 IX = 0,NX
        DO 80 IY = 0,NY+1
          UO(IX,IY) = UN(IX,IY)
   80   CONTINUE
   70 CONTINUE
*--------------------------------------------------------------------
* VN -> VO : 必要なら入れ替える前にVNとVOから変動量を求める
* VN : 前の時間ステップにおいて最終的に得られた値，圧力補正の度に更新される
* VO : 新しい時間ステップでの初期値．VNを保存．
*--------------------------------------------------------------------
      DO 90 IX = 0,NX+1
        DO 100 IY = 0,NY
          VO(IX,IY) = VN(IX,IY)
  100   CONTINUE
   90 CONTINUE
*--------------------------------------------------------------------
* TN -> TO : 必要なら入れ替える前にTNとTOから変動量を求める
* TN : 前の時間ステップでの値
* TO : 新しい時間ステップでの初期値．TNを保存．
*--------------------------------------------------------------------
      DO 110 IX = 0,NX+1
        DO 120 IY = 0,NY+1
          TO(IX,IY) = TN(IX,IY)
  120   CONTINUE
  110 CONTINUE
      RETURN
      END
*********************************************************************
*                    速度場の計算
*********************************************************************
      SUBROUTINE CALVEL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
*--------------------------------------------------------------------
*                     U(IX,IY)の計算
*--------------------------------------------------------------------
      DO 10 IX = 1,NX-1
        DO 20 IY = 1,NY
*       VVはU(IX,IY)におけるVの補間値
        VV = (  VO(IX,IY  )+VO(IX+1,IY  )
     $         +VO(IX,IY-1)+VO(IX+1,IY-1) )/4.0D0
*       対流項(CNVUX,CNVUY)を１次精度風上差分にて計算
        IF ( UO(IX,IY). GE. 0.0D0 ) THEN
          CNVUX = UO(IX,IY)*( UO(IX,IY) - UO(IX-1,IY) ) / DX
        ELSE IF ( UO(IX,IY). LT. 0.0D0 ) THEN
          CNVUX = UO(IX,IY)*( UO(IX+1,IY) - UO(IX,IY) ) / DX
        END IF
        IF ( VV. GE. 0.0D0 ) THEN
          CNVUY = VV*( UO(IX,IY) - UO(IX,IY-1) ) / DY
        ELSE IF ( VV. LT. 0.0D0 ) THEN
          CNVUY = VV*( UO(IX,IY+1) - UO(IX,IY) ) / DY
        END IF
*       x方向の浮力項(BUOU)はゼロ
        TU = 0.0D0
        BUOU = BUO * TU
*       拡散項(DIFU)の計算
        DIFU = VIS*(
     $        ( UO(IX-1,IY)-2.0D0*UO(IX,IY)+UO(IX+1,IY) )/DX**2
     $       +( UO(IX,IY-1)-2.0D0*UO(IX,IY)+UO(IX,IY+1) )/DY**2
     $              )
*       仮の速度(U)の計算
        UN(IX,IY) = UO(IX,IY)
     $  + DT*( -CNVUX-CNVUY+DIFU+BUOU+( PO(IX,IY)-PO(IX+1,IY) )/DX )
   20   CONTINUE
   10 CONTINUE
*--------------------------------------------------------------------
*                       V(IX,IY)の計算
*--------------------------------------------------------------------
      DO 30 IX = 1,NX
        DO 40 IY = 1,NY-1
*       UUはV(IX,IY)におけるUの補間値
        UU = (  UO(IX-1,IY  )+UO(IX,IY  )
     $         +UO(IX-1,IY+1)+UO(IX,IY+1) )/4.0D0
*       対流項(CNVVX,CNVVY)を１次精度風上差分にて計算
        IF ( UU. GE. 0.0D0 ) THEN
          CNVVX = UU*( VO(IX,IY) - VO(IX-1,IY) ) / DX
        ELSE IF ( UU. LT. 0.0D0 ) THEN
          CNVVX = UU*( VO(IX+1,IY) - VO(IX,IY) ) / DX
        END IF
        IF ( VO(IX,IY). GE. 0.0D0 ) THEN
          CNVVY = VO(IX,IY)*( VO(IX,IY) - VO(IX,IY-1) ) / DY
        ELSE IF ( VO(IX,IY). LT. 0.0D0 ) THEN
          CNVVY = VO(IX,IY)*( VO(IX,IY+1) - VO(IX,IY) ) / DY
        END IF
*       浮力項(BUOV)の計算
        TV = ( TO(IX,IY) + TO(IX,IY+1) )/2.0D0
        BUOV = BUO*TV
*       拡散項(DIFV)の計算
        DIFV = VIS*(
     $          ( VO(IX-1,IY)-2.0D0*VO(IX,IY)+VO(IX+1,IY) )/DX**2
     $         +( VO(IX,IY-1)-2.0D0*VO(IX,IY)+VO(IX,IY+1) )/DY**2
     $             )
*       仮の速度(V)の計算
        VN(IX,IY) = VO(IX,IY)
     $  + DT*(-CNVVX-CNVVY+DIFV+BUOV+(PO(IX,IY)-PO(IX,IY+1))/DY )
   40   CONTINUE
   30 CONTINUE
*速度の境界条件の処理
      CALL VELBND
      RETURN
      END
*********************************************************************
*                        圧力場の計算
*********************************************************************
      SUBROUTINE PRESS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
      IXMAX = 0
      IYMAX = 0
      DMAX = 0.0D0
*P(IX,IY)の計算
      DO 10 IX = 1,NX
        DO 20 IY = 1,NY
          DEL = DT*( 2.0D0/DX**2 + 2.0D0/DY**2 )
          DIV = ( UN(IX,IY) - UN(IX-1,IY  ) )/DX
     $        + ( VN(IX,IY) - VN(IX  ,IY-1) )/DY
          IF ( DABS(DIV). GE. DABS(DMAX) ) THEN
            IXMAX = IX
            IYMAX = IY
            DMAX = DIV
          END IF
          DELP = - OMG * DIV / DEL
          PO(IX  ,IY  ) = PO(IX  ,IY  ) + DELP
          UN(IX  ,IY  ) = UN(IX  ,IY  ) + DT/DX*DELP
          UN(IX-1,IY  ) = UN(IX-1,IY  ) - DT/DX*DELP
          VN(IX,  IY  ) = VN(IX,  IY  ) + DT/DY*DELP
          VN(IX,  IY-1) = VN(IX,  IY-1) - DT/DY*DELP
   20   CONTINUE
   10 CONTINUE
*
* 圧力の相対性に関する処理(IRELP=1なら以下の処理を行う)
      IF (IRELP.EQ.1) THEN
        POSTN = PO(1,1)
        DO 30 IX = 1,NX
          DO 40 IY = 1,NY
            PO(IX,IY) = PO(IX,IY)-POSTN
   40     CONTINUE
   30   CONTINUE
      END IF
*
* IFLG=1なら，連続の式を満たしていないと判定し再び圧力計算を行う．
      IF ( DABS(DMAX). GE. EPSP ) IFLG = 1
* 圧力計算の回数を100回ごとに表示
      IF (MOD(ITR,100).EQ.0) THEN
        WRITE (6,2000) ITR, IXMAX, IYMAX, DMAX
 2000   FORMAT (' Iteration=',I8,'   Div(max)(',2I6,')=',1PE13.5)
      END IF
*新たに得られた速度を用いて境界条件を処理する
      CALL VELBND
      RETURN
      END
*********************************************************************
*                      温度場の計算
*********************************************************************
      SUBROUTINE CALTEM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
*T(IX,IY)の計算
      DO 10 IX = 1,NX
        DO 20 IY = 1,NY
*         UUT,VVTはそれぞれT(IX,IY)におけるU,Vの補間値
          UUT = ( UO(IX,IY) + UO(IX-1,IY  ) ) / 2.0D0
          VVT = ( VO(IX,IY) + VO(IX  ,IY-1) ) / 2.0D0
*         対流項(CNVTX,CNVTY)を１次精度風上差分にて計算
          IF ( UUT. GE. 0.0D0 ) THEN
            CNVTX = UUT*( TO(IX,IY) - TO(IX-1,IY) ) / DX
          ELSE IF ( UUT. LT. 0.0D0 ) THEN
            CNVTX = UUT*( TO(IX+1,IY) - TO(IX,IY) ) / DX
          END IF
          IF ( VVT. GE. 0.0D0 ) THEN
            CNVTY = VVT*( TO(IX,IY) - TO(IX,IY-1) ) / DY
          ELSE IF ( VVT. LT. 0.0D0 ) THEN
            CNVTY = VVT*( TO(IX,IY+1) - TO(IX,IY) ) / DY
          END IF
*         拡散項(DIFT)の計算
          DIFT = ALP*(
     $       +( TO(IX-1,IY)-2.0D0*TO(IX,IY)+TO(IX+1,IY) )/DX**2
     $       +( TO(IX,IY-1)-2.0D0*TO(IX,IY)+TO(IX,IY+1) )/DY**2
     $                )
*         次の時間のTの計算
          TN(IX,IY) = TO(IX,IY) + DT*( -CNVTX-CNVTY+DIFT )
   20   CONTINUE
   10 CONTINUE
*境界条件の処理
      CALL TBND
      RETURN
      END
*********************************************************************
*                    速度の境界条件の処理
*********************************************************************
      SUBROUTINE VELBND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
*U（右面）
      DO 10 IY = 1,NY
        UN(NX,IY) = 0.0D0
   10 CONTINUE
*U（左面）
      DO 20 IY = 1,NY
        UN(0,IY) = 0.0D0
   20 CONTINUE
*U（上面）
      DO 30 IX = 0,NX
        UN(IX,NY+1) = -UN(IX,NY)
   30 CONTINUE
*U（下面）
      DO 40 IX = 0,NX
        UN(IX,0) = -UN(IX,1)
   40 CONTINUE
*V（右面）
      DO 50 IY = 1,NY-1
        VN(NX+1,IY) = -VN(NX,IY)
   50 CONTINUE
*V（左面）
      DO 60 IY = 1,NY-1
        VN(0,IY) = -VN(1,IY)
   60 CONTINUE
*V（上面）
      DO 70 IX = 0,NX+1
        VN(IX,NY) = 0.0D0
   70 CONTINUE
*V（下面）
      DO 80 IX = 0,NX+1
        VN(IX,0) = 0.0D0
   80 CONTINUE
      RETURN
      END
*********************************************************************
*                  温度の境界条件の処理
*********************************************************************
      SUBROUTINE TBND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
*右面
      DO 10 IY = 0,NY+1
        TN(NX+1,IY) = 2.0D0 * ( -0.5D0 ) - TN(NX,IY)
   10 CONTINUE
*左面
      DO 20 IY = 0,NY+1
        TN(0,IY) = 2.0D0 * ( +0.5D0 ) - TN(1,IY)
   20 CONTINUE
*上面
      DO 30 IX = 1,NX
        TN(IX,NY+1) = TN(IX,NY)
   30 CONTINUE
*下面
      DO 40 IX = 1,NX
        TN(IX,0) = TN(IX,1)
   40 CONTINUE
      RETURN
      END
*********************************************************************
*                       データ出力
*********************************************************************
      SUBROUTINE PROUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
      WRITE (11) UN
      WRITE (12) VN
      WRITE (13) PO
      WRITE (14) TN
      RETURN
      END
*********************************************************************
*                       Tecplot用データ出力
*********************************************************************
      SUBROUTINE TECPLT (FNAME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
      CHARACTER FNAME*20
*
      OPEN (21,FILE=FNAME,STATUS='NEW')
      WRITE (21,*) 'VARIABLES = "X", "Y", "U", "V", "T"'
      NX1 = NX+1
      NY1 = NY+1
      WRITE (21,4000) NX1,NY1
 4000 FORMAT (1H ,'ZONE I=',I3,',J=',I3,', F=POINT')
      DO 10 IY = 0,NY
        DO 20 IX = 0,NX
           X = DX * FLOAT(IX)
           Y = DY * FLOAT(IY)
           U = ( UN(IX  ,IY)+UN(IX  ,IY+1) )/2.0D0
           V = ( VN(IX  ,IY)+VN(IX+1,IY  ) )/2.0D0
           T = ( TN(IX  ,IY)  +TN(IX+1,IY  )
     $          +TN(IX  ,IY+1)+TN(IX+1,IY+1) )/4.0D0
          WRITE (21,4010) X,Y,U,V,T
 4010     FORMAT (1H ,5(1PE11.3))
   20   CONTINUE
   10 CONTINUE
      CLOSE(21)
      RETURN
      END
