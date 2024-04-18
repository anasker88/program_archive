*********************************************************************
* ファイル名  ：SMAC3D.FOR                                          *
* タイトル    ：SMAC法による3次元熱流動解析プログラム               *
* 製作者      ：平野　博之                                          *
* 所属        ：岡山理科大学 工学部 バイオ・応用化学科              *
* 製作日      ：2011.11.01                                          *
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
*    ITYPE......1:等温計算 2:非等温計算(１３，１４行を参照)         *
*    ICYCLE.....計算開始のサイクル数（時間T=ICYCLE*DT）             *
*            0なら新規でプログラムにある初期条件にしたがって計算開始*
*            0以外の値なら継続計算                                  *
*            時刻 TIME = ICYCLE * DT                                *
*    NITR....圧力補正の線形システム解法(反復法・クリロフ部分空間法) *
*            のための最大反復回数                                   *
*    NCYCLE..計算終了サイクル数                                     *
*    EPSP....収束判定値                                             *
*            圧力補正の線形システム解法(反復法,クリロフ部分空間法)  *
*            の収束評価に使用                                       *
*    OMG....圧力補正のための緩和係数: (point, line-)SOR法の加速係数 *
*    DT.....時間刻み                                                *
*    RE.....レイノルズ数, PR.....プラントル数, GR.....グラスホフ数  *
*    DLX.....解析領域の横幅(DLX=NX*DX):DXは一定で等間隔格子         *
*    DLY.....解析領域の高さ(DLY=NY*DY):DYは一定で等間隔格子         *
*            格子幅（等間隔）DX,DYはプログラムの中で求める．        *
*   IRELP...圧力の基準値と解の1次独立性の設定                       *
*           0 : 行わない                                            *
*             <SMACにおける圧力補正の線形システム解法に関して>      *
*             -> 1次従属な解の1つを求めるのみ                       *
*             -> 直接法では特異行列の問題に遭遇. 学的には正しくない.*
*                丸め誤差がなければ解は得られない．                 *
*           1 : 圧力基準を反映した1次独立な解を求める．             *
*           2 : 圧力基準を設定 P(1,1)=0                             *
*               1次従属の解の1つを求め,圧力基準値P(1,1)=0を設定する.*
*         (IRELP=0の計算でPDを求め，PD(1,1)を差し引きP(1,1)=0とする)*
*    METHOD..圧力補正の線形システム解法のアルゴリズム               *
*            すべてバンドマトリックス用に最適化してある             *
*            1: 直接法 -> ガウスの消去法                            *
*               IRELP=1とする必要がある．(丸め誤差により，IREP=0,2  *
*               としても解を得られる場合もあるが数学的に正しくない.)*
*            2: 反復法1 -> point-SOR 法                             *
*            3: 反復法2 -> line-SOR 法                              *
*               OMG=1で十分．大きくしすぎると発散する               *
*            4: クリロフ部分空間法1 -> 共役残差法                   *
*            5: クリロフ部分空間法2 -> Bi-CGSTAB法                  *
*    パラメータファイルの数値について，                             *
*    FORTRAN プログラム.....できれば倍精度実数で与える"1.0D0,1.0d0" *
*            Cプログラムと共用させて"1.0E0,1.0e0"としても問題はない *
*    C プログラム....."1.0E0,1.0e0"として与える                     *
*    TECPLOT用データを除いた入出力ファイルは書式なし形式で，使用    *
*    するコンパイラーに依存する．コンパイラーに依存しない形式にする *
*    には，容量は増えるが書式付き形式にすればよい．                 *
*                                                                   *
*  ●圧力の相対性について（本文：14.5.3項を参照）                   *
*    圧力補正の線形システムを解くにあたり，本問題のような，境界にお *
*    ける速度が既知の場合，係数行列が特異(matrix singular)となり，  *
*    得られる解は一次従属となり，無数の解が存在することとなる．     *
*    これに対する方策として，                                       *
*    (1) 直接法を用いるときは，方程式の数を減らす必要がある．       *
*    サブルーチンPRESSにおいて，基準として δP(1,1)=P(1,1)=0 と固定 *
*    できるようになっている．                                       *
*    (2) クリロフ部分空間法を含む反復法においては，意識しなくてよい *
*    場合が多い．                                                   *
*    数学的に厳密にいえば，一次従属な解のうちの1つを見つけるという  *
*    ことは，正しいとはいえないが，圧力の値は，絶対値ではなく相対値 *
*    のみが問題となる．                                             *
*    ただし，計算を進行させてゆくにつれ，圧力の絶対値が大きくなって *
*    ゆくような場合は，一次従属な解のうちの1つを求めてから，基準値を*
*    設定したほうがよい．はじめから，基準値を設定して一次独立な解を *
*    得るには，通常，計算時間が長くなる．                           *
*                                                                   *
*  ●変数・配列の説明                                               *
*    ICYCLE  -----> 時間進行のためのカウンタ                        *
*    ITR     -----> 圧力補正計算のための反復回数のカウンタ          *
*                 （クリロフ部分空間法を含む反復法において使用）    *
*    IX,IY,IZ   -----> 上の図を参照                                 *
*    UO,VO,WO,TO-----> 圧力補正計算を行う前の値                     *
*    UN,VN,WN,TN-----> 新たな圧力補正を用いて計算された値           *
*    PD         -----> 圧力補正                                     *
*                                                                   *
*  ●SMAC法における線形システム解法について                         *
*    静止状態直後など，速度場がゼロであったりすると，METHOD=1,4,5に *
*    おいて解を得られない場合がある．                               *
*    METHOD=1: 一次独立な解を求めるのが数学的に正しいのでIRELP=1と  *
*              とすること．コンパイラーによっては，IRELP=0,2として  *
*              も丸め誤差などにより，解が得られる場合もある．       *
*    METHOD=4: 初期値ベクトル\vec{X}_{0}がゼロとなったりすると，探索*
*              ベクトルの計算がゼロ除算となるので適当に初期値を設定 *
*              する必要がある．本プログラムでは，初期値として最新の *
*              値を用いるようにして，可能な限り繰り返し回数が少なく *
*              なることを優先としており，ゼロ除算のときは POINT-SOR *
*              法に切り替えるようにしてあるので，このときは，OMGの  *
*              値を用いることとなる．                               *
*    METHOD=5: METHOD=4に同じ．                                     *
*                                                                   *
*    本プログラムは，HAMAC法のプログラムhsmac3d.forをもとにして作成 *
*    してある．HSMACとSMACが本質的に同等であるので，プログラムの変更*
*    は，圧力補正をニュートン法で行う(HSMAC法)代わりに，線形システム*
*    解法で行うようにすればよい．なお，できる限り，hsmac3d.forとの違*
*    いをわかりやすくするため，SMAC法において新たに付け加えた個所に *
*    は，"--- SMAC ---"としてコメントを挿入してある．               *
*                                                                   *
*********************************************************************
      PROGRAM SMAC3D
*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*     圧力補正の線形システム用
      DIMENSION A1(NE0),A2(NE0),A3(NE0),A4(NE0),A5(NE0),A6(NE0),
     $          A7(NE0),B(NE0),X(NE0)
*     作業用配列
      DIMENSION W1(NE0),W2(NE0),W3(NE0),W4(NE0),W5(NE0),W6(NE0),
     $          W7(NE0),W8(NE0),W9(NE0)
*     バンドガウス消去法のための係数行列配列
      DIMENSION AGB(-NXY0:NXY0,NE0)
*
      CHARACTER FNAME(12)*20
*
*x方向の格子分割数
      NX  = NX0
*y方向の格子分割数
      NY  = NY0
*Z方向の格子分割数
      NZ  = NZ0
*方程式の数(圧力補正に関する未知数の数) --- SMAC ---
      NE  = NE0
* NX * NY
      NXY=NXY0
*パラメータファイルのオープン
      OPEN (10,FILE='in3d.mac',STATUS='OLD')
*出力ファイル名の読み込み
      DO 10 I = 1,11
        READ (10,'(A20)') FNAME(I)
   10 CONTINUE
*Uの計算結果出力用ファイルオープン(書式なし形式)
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
      READ (10,*) ITYPE, ICYCLE0, NITR, NCYCLE
      ICYCLE = ICYCLE0
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
*浮力項の係数( ここでは Pr * Ra = Pr * ( Pr * Gr) )
      BUO = ( GR * PR**2 )
*等温場なら浮力項の係数はゼロに設定
      IF (ITYPE.EQ.1) BUO = 0.0D0
*初期値の設定
      CALL CINITI
      DO 20 I=1,NE
        A1(I)=0.0D0
        A2(I)=0.0D0
        A3(I)=0.0D0
        A4(I)=0.0D0
        A5(I)=0.0D0
        A6(I)=0.0D0
        A7(I)=0.0D0
        B(I)=0.0D0
        X(I)=0.0D0
        W1(I)=0.0D0
        W2(I)=0.0D0
        W3(I)=0.0D0
        W4(I)=0.0D0
        W5(I)=0.0D0
        W6(I)=0.0D0
        W7(I)=0.0D0
        W8(I)=0.0D0
        W9(I)=0.0D0
        DO 25 II = -NXY,NXY
          AGB(II,I)=0.0D0
   25   CONTINUE
   20 CONTINUE
*時間進行のための戻り点
  700 CONTINUE
*時間進行
      CALL ADV
*圧力補正の計算が収束したかどうかのパラメータ(IFLG)を初期化(反復法でのみ有効)
*IFLG -> 0:収束 1:発散(設定された許容回数NITR以下で解が得られない)
      IFLG = 0
C     --- SMAC ---
C     本計算においては，ICYCLE=0のときの初期条件は速度場と圧力場をゼロと
C     しているので，圧力補正の線形システムを解くのが困難となるため，
C     ICYCLE=1のときだけは温度の計算へジャンプするようにする．
C     計算条件や問題に応じて適宜削除あるいは変更
      IF (ICYCLE.EQ.1) GOTO 720
*速度場の計算
      CALL CALVEL
*圧力場の計算
C     CALL PRESS
      CALL PRESS (A1,A2,A3,A4,A5,A6,A7,B,X,
     $            AGB,W1,W2,W3,W4,W5,W6,W7,W8,W9,NXY)
C     速度場と圧力場がゼロとなるときのスキップ先--- SMAC ---
  720 CONTINUE
*     圧力場の計算が収束したとき
      IF ( IFLG. EQ. 0 ) THEN
*       非等温場計算の場合
        IF ( ITYPE.EQ.2 ) THEN
*         温度場を計算
          CALL CALTEM
        END IF
*     圧力場の計算が収束していないとき
      ELSE IF ( IFLG. EQ. 1 ) THEN
        WRITE (6,*) ' NOT CONVERGE ! '
C       データを出力して強制終了
        CALL PROUT
        GO TO 900
      END IF
*     時間進行カウンタ(ICYCLE)がNCYCLEより小さい時
      IF ( ICYCLE. LT. NCYCLE ) THEN
        GO TO 700
*     時間進行サイクルがNCYCLEになったら計算終了
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
*
*********************************************************************
*                          初期設定
*********************************************************************
      SUBROUTINE CINITI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
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
C           --- SMAC ---
            PD(IX,IY,IZ) = 0.0D0
            PN(IX,IY,IZ) = 0.0D0
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
C     --- SMAC ---
        READ (19) PN,PD
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
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*
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
* VO : 新しい時間ステップでの初期値，VNを保存．
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
* WO : 新しい時間ステップでの初期値，VNを保存．
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
* TN : 前の時間ステップでの値
* TO : 新しい時間ステップでの初期値．TNを保存．
* PN -> PO : 必要なら入れ替える前にPNとPOから変動量を求める
* PN : 前の時間ステップでの値
* PO : 新しい時間ステップでの初期値．PNを保存．
*--------------------------------------------------------------------
      DO 220 IX = 0,NX+1
        DO 230 IY = 0,NY+1
          DO 240 IZ = 0,NZ+1
            TO(IX,IY,IZ) = TN(IX,IY,IZ)
C           --- SMAC ---
            PO(IX,IY,IZ) = PN(IX,IY,IZ)
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
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*
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
*                                                                   *
*［圧力補正の線形システムに関する配列の説明］                       *
*                                                                   *
* 圧力補正(PD)の線形システムを A_{i,j} PD_{i} = B_{i} とする．      *
* 有限差分近似を用いて離散化していることから明らかなように，        *
* 1. A_{i,j}の大部分はゼロ                                          *
* 2. 非ゼロ要素は疎であり，規則的に並んでいる                       *
* したがって，A_{i,j}すべてを記憶させるのは計算容量のムダであるので *
* 本プログラムでは，以下のような規則にしたがって記憶する．          *
* 3次元問題であれば，圧力補正に関するポアソン方程式において，2階の  *
* 微分項は近隣の7点で表せる．                                       *
*                                                                   *
* ●用いる配列(1次元配列に格納)：クリロフ部分空間法を含む反復法     *
*   係数行列     ---> A1, A2, A3, A4, A5, A6, A7                    *
*   既知ベクトル ---> B                                             *
*   未知ベクトル ---> X                                             *
*   未知数の数(未知のPDの数)   ---> NX * NY * NZ                    *
*                                                                   *
*    A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE)               *
*                                                                   *
*    i=3, j=3, k=3 (m=59): PD(3,3,3) における A1,A2,A3,A4,A5,A6,A7  *
*                                                                   *
*        +------------------/-----+                                 *
* (NY)5  |    |    |    |  / |    |                                 *
*        +---------------A7-------+    通し番号の規則               *
*     4  |    |    |A4  /    |    |    m=(k-1)*(NX*NY)+(i-1)*NY+j   *
*        +-------------/----------+                                 *
*     3  |    | A1 |A3  | A5 |    |    A1(m)                        *
*        +---------/--------------+    A2(m)                        *
*     2  |    |   /|A2  |    |    |    A3(m)                        *
*        +------A6----------------+    A4(m)                        *
*     1  |      /  |    |    |    |    A5(m)                        *
*        +-----/------------------+    A6(m)                        *
*     j   i=1 /   2    3    4    5(NX) A7(m)                        *
*            /                         B(m)                         *
*         k=1,5(NZ)                    NXY = NX * NY                *
*                                                                   *
*    B(NE) : 既知ベクトル                                           *
*                                                                   *
*    X(NE) : 未知ベクトル(各サブルーチンでこれが求まる)             *
*                                                                   *
* ●用いる配列：ガウスの消去法による直接法(バンド行列として格納)    *
*    NXY = NX * NY                                                  *
*    A(-NXY:NXY,NE) : 3次元差分近似による規則的非対称行列           *
*          -NXY:NXY->上の図で考えるとM=59のときA1からA7までのmは    *
*                   A1のmは M"-NY"                                  *
*                   A2のmは M"-1 "                                  *
*                   A3のmは M"+0 "                                  *
*                   A4のmは M"+1 "                                  *
*                   A5のmは M"+NY"                                  *
*                   A6のmは M"-(NX*NY)"                             *
*                   A7のmは M"+(NX*NY)"                             *
*              このように，mを固定したとき，" "で囲まれた7個の値    *
*              (-NY,-1,0,1,NY,-NXY,NXY)                             *
*     NE -> 上述のような m は全部で NE 個定義される                 *
*                                                                   *
* ●線形システムを構成する際の注意                                  *
* 線形システムを構成する際，境界条件をどこで反映させるかによって，  *
* プログラミングが異なる．本プログラムでは，境界条件は，係数行列を  *
* 作成する際に反映させ，線形システムの解法においては専ら AX=B のみ  *
* に着目する立場をとる．                                            *
* 上述のような立場とは異なり，線形システムの解法において境界条件を  *
* 反映させることもできるが，ここでは，線形システムの解法のサブルー  *
* チンに汎用性をもたせることを優先させた．                          *
*                                                                   *
*********************************************************************
      SUBROUTINE PRESS (A1,A2,A3,A4,A5,A6,A7,B,X,
     $                  AGB,W1,W2,W3,W4,W5,W6,W7,W8,W9,NXY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*     圧力補正の線形システム用
      DIMENSION A1(NE0),A2(NE0),A3(NE0),A4(NE0),A5(NE0),
     $          A6(NE0),A7(NE0),B(NE0),X(NE0)
*     作業用配列
      DIMENSION W1(NE0),W2(NE0),W3(NE0),W4(NE0),W5(NE0),
     $          W6(NE0),W7(NE0),W8(NE0),W9(NE0)
*     バンドガウス消去法のための係数行列配列
      DIMENSION AGB(-NXY:NXY,NE0)
*
*線形システム関係の配列の初期化
      DO 1 IX = 1,NX
        DO 2 IY = 1,NY
          DO 3 IZ = 1,NZ
          K = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          B(K) = 0.0D0
          A1(K) = 0.0D0
          A2(K) = 0.0D0
          A3(K) = 0.0D0
          A4(K) = 0.0D0
          A5(K) = 0.0D0
          A6(K) = 0.0D0
          A7(K) = 0.0D0
          X(K)  = 0.0D0
    3     CONTINUE
    2   CONTINUE
    1 CONTINUE
*P(IX,IY,IZ)の計算
      DO 10 IX = 1,NX
        DO 20 IY = 1,NY
          DO 30 IZ = 1,NZ
          DIV(IX,IY,IZ) = ( UN(IX,IY,IZ) - UN(IX-1,IY  ,IZ  ) )/DX
     $                  + ( VN(IX,IY,IZ) - VN(IX  ,IY-1,IZ  ) )/DY
     $                  + ( WN(IX,IY,IZ) - WN(IX  ,IY  ,IZ-1) )/DZ
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
*     係数行列の作成
*         計算領域の1点さらに内側の点（仮想セルを含めて考えると2点を除く内側の点）
*         に関して係数行列を作成．残りは境界条件を反映させて PDBNDC で設定する
      DO 40 IX = 2,NX-1
        DO 50 IY = 2,NY-1
          DO 60 IZ = 2,NZ-1
            I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
            A3(I) = DT*( 2.0D0/DX**2 + 2.0D0/DY**2 + 2.0D0/DZ**2 )
            A4(I) = - 1.0D0 / DY**2 * DT
            A2(I) = - 1.0D0 / DY**2 * DT
            A1(I) = - 1.0D0 / DX**2 * DT
            A5(I) = - 1.0D0 / DX**2 * DT
            A6(I) = - 1.0D0 / DZ**2 * DT
            A7(I) = - 1.0D0 / DZ**2 * DT
            B(I) = -DIV(IX,IY,IZ)
            X(I) = PD(IX,IY,IZ)
   60     CONTINUE
   50   CONTINUE
   40 CONTINUE
* --- SMAC --- 境界条件を係数行列に反映させる
      CALL PDBNDC (A1,A2,A3,A4,A5,A6,A7,B)
* --- SMAC --- 圧力補正 P'(PD) に関するポアソン方程式の解法
*サブルーチンに一般性を持たせるため，NX,NY,NZ,NEを引数として渡す
*COMMON文で定義されている値でもあるので，そのまま渡せない！
      NNX = NX
      NNY = NY
      NNZ = NZ
      NNE = NE
      NXY=NX*NY
C   1. 直接法 : バンドマトリックスによるガウスの消去法
      IF (METHOD.EQ.1)
     $  CALL GB      (A1,A2,A3,A4,A5,A6,A7,B,X,NNX,NNY,NNZ,NNE,NXY,
     $                AGB)
C   2. 反復法1 : point-SOR 法
      IF (METHOD.EQ.2)
     $  CALL PSORB   (A1,A2,A3,A4,A5,A6,A7,B,X,NNX,NNY,NNZ,NNE,NXY)
C   3. 反復法2 : line-SOR 法 : OMG=1で十分．大きくしすぎると発散する
      IF (METHOD.EQ.3)
     $  CALL LSORB   (A1,A2,A3,A4,A5,A6,A7,B,X,NNX,NNY,NNZ,NNE,NXY,
     $                W1,W2,W3,W4,W5,W6,W7,W8,W9)
C   4. クリロフ部分空間法1 : 共役残差法
      IF (METHOD.EQ.4)
     $  CALL CRB     (A1,A2,A3,A4,A5,A6,A7,B,X,NNX,NNY,NNZ,NNE,NXY,
     $                W1,W2,W3,W4,W5,W6)
C   5. クリロフ部分空間法2 : Bi-CGSTAB
      IF (METHOD.EQ.5)
     $  CALL BICGB   (A1,A2,A3,A4,A5,A6,A7,B,X,NNX,NNY,NNZ,NNE,NXY,
     $                W1,W2,W3,W4,W5,W6,W7)
C   METHOD=4,5のときの探索ベクトル計算時のゼロ除算対策
      IF (IFLG.EQ.2)
     $  CALL PSORB   (A1,A2,A3,A4,A5,A6,A7,B,X,NNX,NNY,NNZ,NNE,NXY)
      DO 70 IX = 1,NX
        DO 80 IY = 1,NY
          DO 90 IZ = 1,NZ
          K = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
*         圧力の相対性の処理 : 基準値の設定
*         1次独立な解を求める場合は圧力の基準点が強制的に設定されるの
*         圧力の基準点を設けない場合 (IRELP=0,1)
*         下記(1)もしくは(2)のいずれかを選択
*         圧力の基準点を設けない場合 (IRELP=0,1)
          IF (IRELP.EQ.0. OR. IRELP.EQ.1) PD(IX,IY,IZ) = X(K)
*         1次従属な解のうちの1つを求めた後，圧力基準を設ける場合 (IRELP=2)
*         P'(1,1,1)=0 ---> P(1,1,1)=0 にセットする場合
          IF (IRELP.EQ.2) PD(IX,IY,IZ) = X(K)-X(1)
   90     CONTINUE
   80   CONTINUE
   70 CONTINUE
*圧力補正の境界条件の処理
      CALL PDBND
*速度の修正
      DO 100 IX = 1,NX-1
        DO 110 IY = 1,NY
          DO 120 IZ = 1,NZ
          UN(IX,IY,IZ)=UN(IX,IY,IZ)+(PD(IX,IY,IZ)-PD(IX+1,IY,IZ))/DX*DT
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE
      DO 130 IX = 1,NX
        DO 140 IY = 1,NY-1
          DO 150 IZ = 1,NZ
          VN(IX,IY,IZ)=VN(IX,IY,IZ)+(PD(IX,IY,IZ)-PD(IX,IY+1,IZ))/DY*DT
  150     CONTINUE
  140   CONTINUE
  130 CONTINUE
      DO 160 IX = 1,NX
        DO 170 IY = 1,NY-1
          DO 180 IZ = 1,NZ
          WN(IX,IY,IZ)=WN(IX,IY,IZ)+(PD(IX,IY,IZ)-PD(IX,IY,IZ+1))/DZ*DT
  180     CONTINUE
  170   CONTINUE
  160 CONTINUE
*新たに得られた速度を用いて境界条件を処理する
      CALL VELBND
*圧力の修正
      DO 190 IX = 1,NX
        DO 200 IY = 1,NY
          DO 210 IZ = 1,NZ
            PN(IX,IY,IZ) = PO(IX,IY,IZ) + PD(IX,IY,IZ)
  210     CONTINUE
  200   CONTINUE
  190 CONTINUE
*新たに得られた圧力を用いて境界条件を処理する
      CALL PBND
      RETURN
      END
*********************************************************************
*                     温度場の計算
*********************************************************************
      SUBROUTINE CALTEM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*     圧力補正の線形システム用
      DIMENSION A1(NE0),A2(NE0),A3(NE0),A4(NE0),A5(NE0),
     $          A6(NE0),A7(NE0),B(NE0),X(NE0)
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
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*     圧力補正の線形システム用
      DIMENSION A1(NE0),A2(NE0),A3(NE0),A4(NE0),A5(NE0),
     $          A6(NE0),A7(NE0),B(NE0),X(NE0)
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
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*     圧力補正の線形システム用
      DIMENSION A1(NE0),A2(NE0),A3(NE0),A4(NE0),A5(NE0),
     $          A6(NE0),A7(NE0),B(NE0),X(NE0)
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
*                                                                   *
*                圧力補正の係数行列と境界条件                       *
*                                                                   *
* 本計算プログラムでは，圧力補正の境界条件を，係数行列作成時に考慮  *
* する．また，解の1次独立性の確保もここで処理する．                 *
*                                                                   *
*********************************************************************
      SUBROUTINE PDBNDC (A1,A2,A3,A4,A5,A6,A7,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      DIMENSION A1(NE0),A2(NE0),A3(NE0),A4(NE0),A5(NE0),
     $          A6(NE0),A7(NE0),B(NE0)
*
C     --- SMAC ---
*     係数行列に境界条件を反映させる
*     計算領域の左側の係数行列（境界条件を反映）
      IX = 1
      DO 10 IY = 2,NY-1
        IZ=1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 1.0D0/DX**2 + 2.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
        DO 20 IZ = 2,NZ-1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 1.0D0/DX**2 + 2.0D0/DY**2 + 2.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   20   CONTINUE
        IZ=NZ
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 1.0D0/DX**2 + 2.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   10 CONTINUE
*     計算領域の右側の係数行列（境界条件を反映）
      IX = NX
      DO 30 IY = 2,NY-1
        IZ=1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 1.0D0/DX**2 + 2.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
        DO 40 IZ = 2,NZ-1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 1.0D0/DX**2 + 2.0D0/DY**2 + 2.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   40   CONTINUE
        IZ=NZ
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 1.0D0/DX**2 + 2.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   30 CONTINUE
*     計算領域の下側の係数行列（境界条件を反映）
      IY = 1
      DO 50 IX = 2,NX-1
        IZ=1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
        DO 60 IZ = 2,NZ-1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 2.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   60   CONTINUE
        IZ=NZ
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   50 CONTINUE
*     計算領域の上側の係数行列（境界条件を反映）
      IY = NY
      DO 70 IX = 2,NX-1
        IZ=1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
        DO 80 IZ = 2,NZ-1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 2.0D0/DZ**2 )
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   80   CONTINUE
        IZ=NZ
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   70 CONTINUE
*     計算領域の後面の係数行列（境界条件を反映）
      IZ = 1
      DO 90 IX = 2,NX-1
        IY=1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
        DO 100 IY = 2,NY-1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 2.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
  100   CONTINUE
        IY=NY
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   90 CONTINUE
*     計算領域の前面の係数行列（境界条件を反映）
      IZ = NZ
      DO 110 IX = 2,NX-1
        IY=1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
        DO 120 IY = 2,NY-1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 2.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
  120   CONTINUE
        IY=NY
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
  110 CONTINUE
*     計算領域の左上点列の係数行列（境界条件を反映）
      IX=1
      IY=NY
      IZ=1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A2(I) = - 1.0D0 / DY**2 * DT
        A5(I) = - 1.0D0 / DX**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
      DO 130 IZ = 2,NZ-1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 2.0D0/DZ**2 )
        A2(I) = - 1.0D0 / DY**2 * DT
        A5(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
  130 CONTINUE
      IZ=NZ
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A2(I) = - 1.0D0 / DY**2 * DT
        A5(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
*     計算領域の左下点列の係数行列（境界条件を反映）
      IX=1
      IY=1
      IZ=1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A4(I) = - 1.0D0 / DY**2 * DT
        A5(I) = - 1.0D0 / DX**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
      DO 140 IZ = 2,NZ-1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 2.0D0/DZ**2 )
        A4(I) = - 1.0D0 / DY**2 * DT
        A5(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
  140 CONTINUE
      IZ=NZ
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A4(I) = - 1.0D0 / DY**2 * DT
        A5(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
*     計算領域の右上点列の係数行列（境界条件を反映）
      IX=NX
      IY=NY
      IZ=1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A2(I) = - 1.0D0 / DY**2 * DT
        A1(I) = - 1.0D0 / DX**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
      DO 150 IZ = 2,NZ-1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 2.0D0/DZ**2 )
        A2(I) = - 1.0D0 / DY**2 * DT
        A1(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
  150 CONTINUE
      IZ=NZ
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A2(I) = - 1.0D0 / DY**2 * DT
        A1(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
*     計算領域の右下点列の係数行列（境界条件を反映）
      IX=NX
      IY=1
      IZ=1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A4(I) = - 1.0D0 / DY**2 * DT
        A1(I) = - 1.0D0 / DX**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
      DO 160 IZ = 2,NZ-1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 2.0D0/DZ**2 )
        A4(I) = - 1.0D0 / DY**2 * DT
        A1(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
  160 CONTINUE
      IZ=NZ
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A4(I) = - 1.0D0 / DY**2 * DT
        A1(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
*１次独立な解を得るための処理 : 直接法においては必須
*(IX=1,IY=1,IZ=1 ---> M=1を基準点とし，常にPN(1,1,1)=PD(1,1,1)=0とする)
      IF (IRELP.EQ.1) THEN
        A3(1) = 1.0D0
        A4(1) = 0.0D0
        A5(1) = 0.0D0
        A7(1) = 0.0D0
        B(1)  = 0.0D0
*       M=2 の点の処理(M=1とのリンクを断つ)
        A3(2) = DT*( 1.0D0/DX**2 + 2.0D0/DY**2 + 1.0D0/DZ**2 )
        A4(2) = - 1.0D0 / DY**2 * DT
        A2(2) = 0.0D0
        A5(2) = - 1.0D0 / DX**2 * DT
        A7(2) = - 1.0D0 / DZ**2 * DT
*       M=1+NY の点の処理(M=1とのリンクを断つ)
        A3(1+NY) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A4(1+NY) = - 1.0D0 / DY**2 * DT
        A1(1+NY) = 0.0D0
        A5(1+NY) = - 1.0D0 / DX**2 * DT
        A7(1+NY) = - 1.0D0 / DZ**2 * DT
*       M=1+NX*NY の点の処理(M=1とのリンクを断つ)
        A3(1+NX*NY) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A4(1+NX*NY) = - 1.0D0 / DY**2 * DT
        A5(1+NX*NY) = - 1.0D0 / DX**2 * DT
        A6(1+NX*NY) = 0.0D0
        A7(1+NX*NY) = - 1.0D0 / DZ**2 * DT
      END IF
*
      RETURN
      END
*********************************************************************
*                                                                   *
*                  圧力補正の境界条件の処理                         *
*                                                                   *
* 本計算プログラムでは，係数行列作成時に境界条件を考慮しているので，*
* 実質的には，最終的に仮想セルの値を配列に格納しているにすぎない．  *
* このルーチンはなくともよい．                                      *
* 反復過程で境界条件を考慮する場合は必須で，反復のたびにこの処理を  *
* 行う必要が生じる．                                                *
*                                                                   *
*********************************************************************
      SUBROUTINE PDBND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*右面
      DO 10 IY = 0,NY+1
        DO 20 IZ = 0,NZ+1
          PD(NX+1,IY,IZ) = PD(NX,IY,IZ)
   20   CONTINUE
   10 CONTINUE
*左面
      DO 30 IY = 0,NY+1
        DO 40 IZ = 0,NZ+1
          PD(0,IY,IZ) = PD(1,IY,IZ)
   40   CONTINUE
   30 CONTINUE
*上面
      DO 50 IX = 1,NX
        DO 60 IZ = 1,NZ
          PD(IX,NY+1,IZ) = PD(IX,NY,IZ)
   60   CONTINUE
   50 CONTINUE
*下面
      DO 70 IX = 1,NX
        DO 80 IZ = 1,NZ
          PD(IX,0,IZ) = PD(IX,1,IZ)
   80   CONTINUE
   70 CONTINUE
*後面
      DO 90 IX = 1,NX
        DO 100 IY = 0,NY+1
          PD(IX,IY,NZ+1) = PD(IX,IY,NZ)
  100   CONTINUE
   90 CONTINUE
*前面
      DO 110 IX = 1,NX
        DO 120 IY = 0,NY+1
          PD(IX,IY,0) = PD(IX,IY,1)
  120   CONTINUE
  110 CONTINUE
      RETURN
      END
*********************************************************************
*                     圧力の境界条件の処理
*********************************************************************
      SUBROUTINE PBND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*右面
      DO 10 IY = 0,NY+1
        DO 20 IZ = 0,NZ+1
          PN(NX+1,IY,IZ) = PN(NX,IY,IZ)
   20   CONTINUE
   10 CONTINUE
*左面
      DO 30 IY = 0,NY+1
        DO 40 IZ = 0,NZ+1
          PN(0,IY,IZ) = PN(1,IY,IZ)
   40   CONTINUE
   30 CONTINUE
*上面
      DO 50 IX = 1,NX
        DO 60 IZ = 1,NZ
          PN(IX,NY+1,IZ) = PN(IX,NY,IZ)
   60   CONTINUE
   50 CONTINUE
*下面
      DO 70 IX = 1,NX
        DO 80 IZ = 1,NZ
          PN(IX,0,IZ) = PN(IX,1,IZ)
   80   CONTINUE
   70 CONTINUE
*後面
      DO 90 IX = 1,NX
        DO 100 IY = 0,NY+1
          PN(IX,IY,NZ+1) = PN(IX,IY,NZ)
  100   CONTINUE
   90 CONTINUE
*前面
      DO 110 IX = 1,NX
        DO 120 IY = 0,NY+1
          PN(IX,IY,0) = PN(IX,IY,1)
  120   CONTINUE
  110 CONTINUE
      RETURN
      END
*********************************************************************
*                        データ出力
*********************************************************************
      SUBROUTINE PROUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*
       WRITE (11) UN
       WRITE (12) VN
       WRITE (13) WN
       WRITE (14) PN,PD
       WRITE (15) TN
      RETURN
      END
*********************************************************************
*                  Tecplot用データ出力
*********************************************************************
      SUBROUTINE TECPLT (FNAME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
      CHARACTER FNAME*20
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
*********************************************************************
*                    各種の線形システム解法                         *
*                                                                   *
*  1. ガウスの消去法                                                *
*  2. point-SOR 法                                                  *
*  3. line-SOR 法                                                   *
*  4. 共役残差法                                                    *
*  5. BiCGSTAB法                                                    *
*                                                                   *
*  いずれも，3次元のポアソン方程式を7点差分近似にて離散化した       *
*  線形システムを解くためのもので，最適化してある                   *
*  いずれのサブルーチンも同一引数としてある                         *
*                                                                   *
* [注意]                                                            *
*  1. 収束判定条件は適宜変更のこと．                                *
*  2. 引数の NX,NY,NE と，配列宣言文の A1,A2,A3,A4,A5,A6,A7,B,X は，*
*     いずれもこのサブルーチンがコールされている PRESS において対応 *
*     するものと同じ名前としてあるが，COMMON文で定義していないので，*
*     計算機の中では異なる変数として定義される．本計算プログラムに  *
*     おいては，できるだけ線形システム解法のサブルーチンに汎用性を  *
*     もたせるため，あえて，COMMON文は使用していない．また，分かり  *
*     やすくするため，サブルーチンがコールされている個所と同じ名前で*
*     それぞれの引数を定義してある．以降同様．                      *
*                                                                   *
*********************************************************************
*
*********************************************************************
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
*********************************************************************
      SUBROUTINE GB (A1,A2,A3,A4,A5,A6,A7,B,X,NX,NY,NZ,NE,NXY,
     $               A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
C     --- SMAC ---
      COMMON / D9 / NITR
*
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE),
     $          B(NE),X(NE)
      DIMENSION A(-NXY:NXY,NE)
*
*直接法のときは(係数行列が特異でなければ)反復なしで，必ず解を得る
      IFLG = 0
*マトリックスAのゼロクリア
      DO 10 INE = 1,NE
        DO 20 I = -NXY,NXY
          A(I,INE) = 0.0D0
   20   CONTINUE
   10 CONTINUE
*
*必要なところにA1からA7までを格納する
      DO 30 INE = 1,NE
        A(- NY,INE) = A1(INE)
        A(  -1,INE) = A2(INE)
        A(   0,INE) = A3(INE)
        A(   1,INE) = A4(INE)
        A(  NY,INE) = A5(INE)
        A(-NXY,INE) = A6(INE)
        A( NXY,INE) = A7(INE)
   30 CONTINUE
*
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
*
*後退代入
*     係数行列の特異性を判定
      IF ( DABS(A(0,NE)).LE.1.0D-50 ) THEN
        WRITE (6,*) ' Matrix is singular : |A(0,NE)| < 1E-50 '
        IFLG = 1
      END IF
      X(NE) = B(NE) / A(0,NE)
*
      DO 90 I = NE-1,1,-1
        S = 0.0D0
        IF ( I.GT.NE-NXY ) THEN
          DO 100 N = 1,NE-I
            S = S + A(N,I)*X(I+N)
  100     CONTINUE
          X(I) = ( B(I)-S ) / A(0,I)
        ELSE
          DO 110 N = 1,NXY
            S = S + A(N,I)*X(I+N)
  110     CONTINUE
          X(I) = ( B(I)-S ) / A(0,I)
        END IF
   90 CONTINUE
*
*サブルーチン終了
      RETURN
      END
*********************************************************************
*                                                                   *
*  point-SOR 法による非対称行列 A を含む線形システム解法サブルーチン*
*  3次元ラプラシアン離散化による7点差分近似用                       *
*                                                                   *
*    線形システム --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    係数行列の計算容量節約：詳細はサブルーチン PRESS を参照        *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE)  *
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
*********************************************************************
      SUBROUTINE PSORB (A1,A2,A3,A4,A5,A6,A7,B,X,NX,NY,NZ,NE,NXY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
C     --- SMAC ---
      COMMON / D9 / NITR
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE),
     $          B(NE),X(NE)
*
*     IFLGの初期値は"収束せず"
      IFLG=1
*
      DO 10 J = 1,NITR
*
        RNORM1 = 0.0D0
*
        I=1
          XOLD = X(I)
          SUM =                           +A4(I)*X(I+1)+A5(I)*X(I+NY)
     $                       +A7(I)*X(I+NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
          RNORM1 = RNORM1 + ( XNEW - XOLD )**2
*
        DO 20 I=2,NY
          XOLD = X(I)
          SUM =               A2(I)*X(I-1)+A4(I)*X(I+1)+A5(I)*X(I+NY)
     $                       +A7(I)*X(I+NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
          RNORM1 = RNORM1 + ( XNEW - XOLD )**2
   20   CONTINUE
*
        DO 30 I=NY+1,NXY
          XOLD = X(I)
          SUM = A1(I)*X(I-NY)+A2(I)*X(I-1)+A4(I)*X(I+1)+A5(I)*X(I+NY)
     $                       +A7(I)*X(I+NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
          RNORM1 = RNORM1 + ( XNEW - XOLD )**2
   30   CONTINUE
*
        DO 40 I=NXY+1,NE-NXY
          XOLD = X(I)
          SUM = A1(I)*X(I-NY   )+A2(I)*X(I-1)+A4(I)*X(I+1)+A5(I)*X(I+NY)
     $         +A6(I)*X(I-NXY)+A7(I)*X(I+NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
          RNORM1 = RNORM1 + ( XNEW - XOLD )**2
   40   CONTINUE
*
        DO 50 I=NE-NXY+1,NE-NY
          XOLD = X(I)
          SUM = A1(I)*X(I-NY   )+A2(I)*X(I-1)+A4(I)*X(I+1)+A5(I)*X(I+NY)
     $         +A6(I)*X(I-NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
          RNORM1 = RNORM1 + ( XNEW - XOLD )**2
   50   CONTINUE
*
        DO 60 I=NE-NY+1,NE-1
          XOLD = X(I)
          SUM = A1(I)*X(I-NY   )+A2(I)*X(I-1)+A4(I)*X(I+1)
     $         +A6(I)*X(I-NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
          RNORM1 = RNORM1 + ( XNEW - XOLD )**2
   60   CONTINUE
*
        I=NE
          XOLD = X(I)
          SUM = A1(I)*X(I-NY)+A2(I)*X(I-1)
     $         +A6(I)*X(I-NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
          RNORM1 = RNORM1 + ( XNEW - XOLD )**2
*
*       収束判定; 収束なら IFLG=0 に設定
        IF (RNORM1.LE.EPSP) THEN
          IFLG=0
          ITR = J
          GO TO 700
        END IF
   10 CONTINUE
*
* 収束と判定されたときの分岐点
  700 CONTINUE
*
*サブルーチン終了
      RETURN
      END
*********************************************************************
*                                                                   *
*  line-SOR 法による非対称行列 A を含む線形システム解法サブルーチン *
*  3次元ラプラシアン離散化による7点差分近似用                       *
*                                                                   *
*    線形システム --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    係数行列の計算容量節約：詳細はサブルーチン PRESS を参照        *
*    A_{i,j} ---> AT1(NE),AT2(NE),AT3(NE),AT4(NE),AT5(NE),          *
*                 AT6(NE),AT7(NE)                                   *
*                                                                   *
*    B -> BX(NE) : 既知ベクトル                                     *
*    X -> XN(NE) : 未知ベクトル ---> これを求める                   *
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
* [配列の説明]                                                      *
* XN...各方向への掃引後のX(番号付けは不変)                          *
*      はじめにこのサブルーチンへ渡されるXでもある                  *
* X1...各方向への掃引後のX(番号付けは軸方向に異なる)                *
* XO...各方向への掃引前のX(番号付けは不変)                          *
* XOLD...このサブルーチンに入る前のX(番号付けは不変)                *
*                                                                   *
*********************************************************************
      SUBROUTINE LSORB (AT1,AT2,AT3,AT4,AT5,AT6,AT7,
     $                  BX,XN,NX,NY,NZ,NE,NXY,
     $                  X1,XO,XOLD,A,B,C,D,U,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
C     --- SMAC ---
      COMMON / D9 / NITR
      DIMENSION AT1(NE),AT2(NE),AT3(NE),AT4(NE),AT5(NE),
     $          AT6(NE),AT7(NE)
      DIMENSION XN(NE),X1(NE),BX(NE),XO(NE),XOLD(NE)
      DIMENSION A(NE),B(NE),C(NE),D(NE),U(NE),Y(NE)
*
*     IFLGの初期値は"収束せず"
      IFLG=1
*
      DO 10 K=1,NITR
*     x 軸方向への掃引 : トーマス法による
      INX = 1
      DO 100 IZ = 1,NZ
        DO 110 IY = 1,NY
          DO 120 IX = 1,NX
          INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*         トーマス法のための係数A,B,C,Dの設定 : XNは最新のX
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
*         トーマス法で答えを求める前のXNをXOに保存
          XO(INZ) = XN(INZ)
*         x方向への掃引の前の値をXOLDに保存
          XOLD(INZ) = XN(INZ)
          INX = INX + 1
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE
*     Ly=b を解く
      U(1) = C(1) / B(1)
      DO 130 J = 2,NE-1
        U(J) = C(J) / ( B(J)-A(J)*U(J-1) )
  130 CONTINUE
      Y(1) = D(1) / B(1)
      DO 140 J = 2,NE
        Y(J) = ( D(J)-A(J)*Y(J-1) ) / ( B(J)-A(J)*U(J-1) )
  140 CONTINUE
*     Ux=y を解く
      X1(NE) = Y(NE)
      DO 150 J = NE-1,1,-1
        X1(J) = Y(J) - U(J)*X1(J+1)
  150 CONTINUE
      INX = 1
      DO 160 IZ = 1,NZ
        DO 170 IY = 1,NY
          DO 180 IX = 1,NX
            INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*           得られたX1と反復前のXOにより最新のXNを緩和
            XN(INZ)=(1.0D0-OMG)*XO(INZ)+OMG*X1(INX)
            INX = INX + 1
  180     CONTINUE
  170   CONTINUE
  160 CONTINUE
*     y 軸方向への掃引 : トーマス法による
      INY = 1
      DO 200 IZ = 1,NZ
        DO 210 IX = 1,NX
          DO 220 IY = 1,NY
          INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*         トーマス法のための係数A,B,C,Dの設定 : XNは最新のX
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
*         トーマス法で答えを求める前のXNをXOに保存
          XO(INZ) = XN(INZ)
          INY = INY + 1
  220     CONTINUE
  210   CONTINUE
  200 CONTINUE
*     Ly=b を解く
      U(1) = C(1) / B(1)
      DO 230 J = 2,NE-1
        U(J) = C(J)/( B(J)-A(J)*U(J-1) )
  230 CONTINUE
      Y(1) = D(1) / B(1)
      DO 240 J = 2,NE
        Y(J) = ( D(J)-A(J)*Y(J-1) ) / ( B(J)-A(J)*U(J-1) )
  240 CONTINUE
*     Ux=y を解く
      X1(NE) = Y(NE)
      DO 250 J = NE-1,1,-1
        X1(J) = Y(J) - U(J)*X1(J+1)
  250 CONTINUE
      INY = 1
      DO 260 IZ = 1,NZ
        DO 270 IX = 1,NX
          DO 280 IY = 1,NY
            INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*           得られたX1と反復前のXOにより最新のXNを緩和
            XN(INZ)=(1.0D0-OMG)*XO(INZ)+OMG*X1(INY)
            INY = INY + 1
  280     CONTINUE
  270   CONTINUE
  260 CONTINUE
*
*     z 軸方向への掃引 : トーマス法による
      INZZ = 1
      DO 300 IY = 1,NY
        DO 310 IX = 1,NX
          DO 320 IZ = 1,NZ
          INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*         トーマス法のための係数A,B,C,Dの設定 : XNは最新のX
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
*         トーマス法で答えを求める前のXNをXOに保存
          XO(INZ) = XN(INZ)
          INZZ = INZZ + 1
  320     CONTINUE
  310   CONTINUE
  300 CONTINUE
*     Ly=b を解く
      U(1) = C(1) / B(1)
      DO 330 J = 2,NE-1
        U(J) = C(J)/( B(J)-A(J)*U(J-1) )
  330 CONTINUE
      Y(1) = D(1) / B(1)
      DO 340 J = 2,NE
        Y(J) = ( D(J)-A(J)*Y(J-1) ) / ( B(J)-A(J)*U(J-1) )
  340 CONTINUE
*     Ux=y を解く
      X1(NE) = Y(NE)
      DO 350 J = NE-1,1,-1
        X1(J) = Y(J) - U(J)*X1(J+1)
  350 CONTINUE
      INZZ = 1
      DO 360 IY = 1,NY
        DO 370 IX = 1,NX
          DO 380 IZ = 1,NZ
            INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*           得られたX1と反復前のXOにより最新のXNを緩和
            XN(INZ)=(1.0D0-OMG)*XO(INZ)+OMG*X1(INZZ)
            INZZ = INZZ + 1
  380     CONTINUE
  370   CONTINUE
  360 CONTINUE
*
      RNORM= 0.0D0
      DO 400 I = 1,NE
        RNORM = RNORM + (XN(I)-XOLD(I))**2
  400 CONTINUE
*
*     収束判定
      IF (RNORM.LE.EPSP) THEN
        IFLG=0
        ITR=K
        GO TO 900
      END IF
*
   10 CONTINUE
*
* 収束と判定されたときの分岐点
  900 CONTINUE
*
      RETURN
      END
*********************************************************************
*                                                                   *
*    共役残差(Conjugate Residual)法による非対称行列 A を含む        *
*  線形システム解法サブルーチン                                     *
*                        AX=B                                       *
*                                                                   *
*    係数行列の計算容量節約：詳細はサブルーチン PRESS を参照        *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE)  *
*                                                                   *
*  B : 既知ベクトル                                                 *
*  X : 未知ベクトル ---> これを求める ---> ここでは便宜上配列 XP(NE)*
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
*    R(NE) : r_{k} = B - A x_{k}                                    *
*    P(NE) :  p_{k+1} = r_{k+1} + β_{k} p_{k}, p_{0} = r_{0}       *
*    AP(NE) : A * P                                                 *
*    AR(NE) : A * R                                                 *
*                                                                   *
* [収束判定条件]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : 前の反復による値                                 *
*  \vec{x}^{new} : 新しい反復による値                               *
*                                                                   *
*********************************************************************
      SUBROUTINE CRB (A1,A2,A3,A4,A5,A6,A7,B,XP,NX,NY,NZ,NE,NXY,
     $                R,P,AP,AR,X,XOLD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
C     --- SMAC ---
      COMMON / D9 / NITR
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE),
     $          B(NE),XP(NE)
*作業用配列
      DIMENSION R(NE),P(NE),AP(NE),AR(NE),X(NE),XOLD(NE)
*
* R に AX を代入
      CALL PROMV (A1,A2,A3,A4,A5,A6,A7,X,R,NX,NY,NZ,NE,NXY)
*
      DO 40 I = 1, NE
*       r_{0}とp_{0}(初期値)の設定
        R(I) = B(I) - R(I)
        P(I) = R(I)
*       前の時刻のXをXOLDに代入
        XOLD(I) = XP(I)
   40 CONTINUE
*
* APに A p_{0} を代入
      CALL PROMV (A1,A2,A3,A4,A5,A6,A7,P,AP,NX,NY,NZ,NE,NXY)
*
* 反復計算
      DO 50 K = 1,NITR
*       ( r_{k}, A p_{k} )の計算 => RAP
        CALL PROVV (R,AP,RAP,NE)
*       ( A p_{k}, A p_{k} )の計算 => APAP
        CALL PROVV (AP,AP,APAP,NE)
*       α_{k} = ( r_{k}, A p_{k} ) / ( A p_{k}, A p_{k} )
**********************************************************
*.......探索方向のための計算が 0除算 なら計算終了(--- SMAC ---)
        IF (DABS(APAP).LT.1.0D-50) THEN
          WRITE (6,*) ' 0 division : ALPHA_{K} in Conjugate Residual '
          IFLG = 2
          RETURN
        ELSE
          ALP = RAP / APAP
        END IF
**********************************************************
*
        RNORM = 0.0D0
*
        DO 70 I = 1,NE
*         x_{k+1}=x_{k}+α_{k}p_{k}
          X(I) = X(I) + ALP*P(I)
*         r_{k+1}=r_{k}-α_{k}Ap_{k}
          R(I) = R(I) - ALP*AP(I)
*         前の反復との差のノルムの計算
          RNORM = RNORM + (X(I)-XOLD(I))**2
*         得られたXをXOLDに代入
          XOLD(I) = X(I)
   70   CONTINUE
*
* RNORM が EPSP 以下なら収束とみなして 700 へ
        IF (RNORM.LE.EPSP) THEN
          IFLG=0
          ITR=K
          GO TO 700
        END IF
*
* 収束せずの場合
*       A r_{k+1} の計算 => AR(NE)
        CALL PROMV (A1,A2,A3,A4,A5,A6,A7,R,AR,NX,NY,NZ,NE,NXY)
*       ( A r_{k+1}, A p_{k} )の計算 => ARAP
        CALL PROVV( AR,AP,ARAP,NE)
*       β_{k} = - ( A r_{k+1}, A p_{k} ) / ( A p_{k}, A p_{k} )
**********************************************************
*.......探索方向のための計算が 0除算 なら計算終了(--- SMAC ---)
        IF (DABS(APAP).LT.1.0D-50) THEN
          WRITE (6,*) ' 0 division : BETA_{K} in Conjugate Residual'
          IFLG = 2
          RETURN
        ELSE
          BETA = -ARAP / APAP
        END IF
**********************************************************
*
        DO 90 I = 1, NE
*         p_{k+1} = r_{k+1} + β_{k} p_{k}
          P(I) = R(I) + BETA*P(I)
*         A p_{k+1} = A r_{k+1} + β_{k}A p_{k}
          AP(I) = AR(I) + BETA*AP(I)
   90   CONTINUE
*
   50 CONTINUE
*     NITR まで計算しても収束せず
      IFLG=1
*
* 収束と判定されたときの分岐点
  700 CONTINUE
*
      DO 100 I = 1,NE
        XP(I) = X(I)
  100 CONTINUE
*
*サブルーチン終了
      RETURN
      END
*
*********************************************************************
*                                                                   *
*    ベクトル A とベクトル B の積の計算サブルーチン                 *
*                        AB=C                                       *
*                                                                   *
*［変数の説明］                                                     *
*    NE : 総格子点数(ベクトル A,B のサイズ)                         *
*    C  : A と B の積(スカラー)                                     *
*                                                                   *
*********************************************************************
      SUBROUTINE PROVV(A,B,C,NE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  A(NE),B(NE)
*
      C = 0.0D0
      DO 10 I=1,NE
        C = C + A(I)*B(I)
   10 CONTINUE
      RETURN
      END
*
*********************************************************************
*                                                                   *
*    マトリックス A とベクトル B の積の計算サブルーチン             *
*                        AB=C                                       *
*                                                                   *
*［変数の説明］                                                     *
*    NE : 総格子点数(正方マトリックス A,B,C のサイズ)               *
*    C  : A と B の積(ベクトル)                                     *
*                                                                   *
*********************************************************************
      SUBROUTINE PROMV (A1,A2,A3,A4,A5,A6,A7,B,C,NX,NY,NZ,NE,NXY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE),
     $          B(NE),C(NE)
*
      I=1
        C(I) =                              A3(I)*B(I)
     $        +A4(I)*B(I+1)  +A5(I)*B(I+NY)
     $                       +A7(I)*B(I+NXY)
*
      DO 10 I=2,NY
        C(I) = A2(I)*B(I-1)                +A3(I)*B(I)
     $        +A4(I)*B(I+1)  +A5(I)*B(I+NY)
     $                       +A7(I)*B(I+NXY)
   10 CONTINUE
*
      DO 20 I=NY+1,NXY
        C(I) = A1(I)*B(I-NY) +A2(I)*B(I-1 )+A3(I)*B(I)
     $        +A4(I)*B(I+1)  +A5(I)*B(I+NY)
     $                       +A7(I)*B(I+NXY)
   20 CONTINUE
*
      DO 30 I=NXY+1,NE-NXY
        C(I) = A1(I)*B(I-NY) +A2(I)*B(I-1 )+A3(I)*B(I)
     $        +A4(I)*B(I+1)  +A5(I)*B(I+NY)
     $        +A6(I)*B(I-NXY)+A7(I)*B(I+NXY)
   30 CONTINUE
*
      DO 40 I=NE-NXY+1,NE-NY
        C(I) = A1(I)*B(I-NY) +A2(I)*B(I-1 )+A3(I)*B(I)
     $        +A4(I)*B(I+1)  +A5(I)*B(I+NY)
     $        +A6(I)*B(I-NXY)
   40 CONTINUE
*
      DO 50 I=NE-NY+1,NE-1
        C(I) = A1(I)*B(I-NY) +A2(I)*B(I-1 )+A3(I)*B(I)
     $        +A4(I)*B(I+1)
     $        +A6(I)*B(I-NXY)
   50 CONTINUE
*
      I=NE
        C(I) = A1(I)*B(I-NY) +A2(I)*B(I-1 )+A3(I)*B(I)
     $        +A6(I)*B(I-NXY)
*
*サブルーチン終了
      RETURN
      END
*
*********************************************************************
*                                                                   *
*  Bi-CGSTAB 法による非対称行列 A を含む                            *
*  線形システム解法サブルーチン                                     *
*                        AX=B                                       *
*                                                                   *
*    係数行列の計算容量節約：詳細はサブルーチン PRESS を参照        *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE)  *
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
*［配列の説明］                                                     *
*                                                                   *
* [収束判定条件]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : 前の反復による値                                 *
*  \vec{x}^{new} : 新しい反復による値                               *
*                                                                   *
*********************************************************************
      SUBROUTINE BICGB (A1,A2,A3,A4,A5,A6,A7,B,X,NX,NY,NZ,NE,NXY,
     $                  R,AT,AP,P,S,T,XOLD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
C     --- SMAC ---
      COMMON / D9 / NITR
*
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE),
     $          B(NE),X(NE)
* 作業用配列
      DIMENSION R(NE),AT(NE),AP(NE),P(NE),S(NE),T(NE),XOLD(NE)
*
      DO 5 I = 1,NE
        XOLD(I) = X(I)
    5 CONTINUE
*
* R に AX を代入
      CALL PROMV (A1,A2,A3,A4,A5,A6,A7,X,R,NX,NY,NZ,NE,NXY)
* r_{0}とp_{0}(初期値)，そして s=r_{0} の設定
      DO 10 I=1,NE
        R(I) = B(I) - R(I)
        P(I) = R(I)
        S(I) = R(I)
   10 CONTINUE
*
* 繰り返し計算
      DO 20 J =1,NITR
*       ( s, r_{k} ) の計算 => SR1
        CALL PROVV (S,R,SR1,NE)
*       A p_{k} の計算 => AP(NE)
        CALL PROMV (A1,A2,A3,A4,A5,A6,A7,P,AP,NX,NY,NZ,NE,NXY)
*       ( s, A p_{k} ) の計算 => SAP
        CALL PROVV (S,AP,SAP,NE)
*       α_{k} = ( s, r_{k} ) / ( s, A p_{k} )
**********************************************************
*.......探索方向のための計算が 0除算 なら計算終了(--- SMAC ---)
        IF (DABS(SAP).LT.1.0D-50) THEN
          WRITE (6,*) ' 0 division : ALPHA_{K} in Bi-CGSTAB'
          IFLG = 2
          RETURN
        ELSE
          ALPHA = SR1/SAP
        END IF
**********************************************************
        DO 50 I=1,NE
*         t_{k} = r_{k} - α_{k} A p_{k}
          T(I) = R(I) - ALPHA*AP(I)
   50   CONTINUE
*       A t_{k} の計算 => AT(NE)
        CALL PROMV (A1,A2,A3,A4,A5,A6,A7,T,AT,NX,NY,NZ,NE,NXY)
*       ( A t_{k}, t_{k} ) の計算 => ATT
        CALL PROVV (AT,T,ATT,NE)
*       ( A t_{k}, A t_{k} ) の計算 => ATAT
        CALL PROVV (AT,AT,ATAT,NE)
*       ξ_{k} = ( A t_{k}, t_{k} ) / ( A t_{k}, A t_{k} )
**********************************************************
*.......探索方向のための計算が 0除算 なら計算終了(--- SMAC ---)
        IF (DABS(ATAT).LT.1.0D-50) THEN
          WRITE (6,*) ' 0 division : XI_{K} in Bi-CGSTAB'
          IFLG = 2
          RETURN
        ELSE
          XI = ATT/ATAT
        END IF
**********************************************************
        RNORM = 0.0D0
        DO 60 I=1,NE
*         x_{k+1} = x_{k} + α_{k} p_{k} + ξ_{k} t_{k}
          X(I) = X(I) + ALPHA*P(I) + XI*T(I)
*         r_{k+1} = t_{k} - ξ_{k} A t_{k}
          R(I) = T(I) - XI*AT(I)
*         前の反復との差のノルムの計算
          RNORM = RNORM + (X(I)-XOLD(I))**2
*         得られたXをXOLDに代入
          XOLD(I) = X(I)
   60   CONTINUE
* RNORM が EPSP 以下なら収束とみなして 900 へ
        IF (RNORM.LE.EPSP) THEN
          IFLG=0
          ITR=J
          GO TO 900
        END IF
*       収束せずの場合
*       ( s, r_{k+1} ) の計算 => SR2
        CALL PROVV (S,R,SR2,NE)
*       β_{k} = ( α_{k} / ξ_{k} ) * ( s, r_{k+1} )/( s, r_{k} )
        BETA = (ALPHA / XI) * (SR2 / SR1)
        DO 70 I=1,NE
*         p_{k+1} = r_{k+1} + β_{k}( p_{k}-ξ_{k} A p_{k} )
          P(I) = R(I) + BETA * ( P(I) - XI*AP(I) )
   70   CONTINUE
   20 CONTINUE
*     NITR まで計算しても収束せず
      IFLG=1
*
  900 CONTINUE
      RETURN
      END
*
