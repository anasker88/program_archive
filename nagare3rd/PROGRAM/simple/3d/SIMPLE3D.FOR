*********************************************************************
* ファイル名  ：SIMPLE3D.FOR                                        *
* タイトル    ：SIMPLE法による3次元熱流動解析プログラム             *
* 製作者      ：平野　博之                                          *
* 所属        ：岡山理科大学 工学部 バイオ・応用化学科              *
* 製作日      ：2011.11.01                                          *
* 言語        ：FORTRAN (FORTRAN77でも実行可能)                     *
*********************************************************************
*                                                                   *
*    本プログラムでは，SIMPLE法のアルゴリズムに基づき,対流項・拡散項*
*  風上法あるいはべき乗法を用いて離散化している．高精度の計算を行う *
*  場合などは，適宜変更のこと．                                     *
*  格子分割数を変更するときはPARAMETER文中の                        *
*  NX0,NY0,NZ0,NEU0,NEV0,NEW0,NEP0をすべて変更.                     *
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
*    上の格子分割図を参考に，たとえば，NX=20,NY=20の場合，PARAMETER *
*  文の変数は以下のように計算する．                                 *
*  NEU0: 境界を除く未知のUの数  : (NX-1) *  NY    *  NZ    = 7600   *
*  NEV0: 境界を除く未知のVの数  :  NX    * (NY-1) *  NZ    = 7600   *
*  NEW0: 境界を除く未知のWの数  :  NX    *  NY    * (NZ-1) = 7600   *
*  NEP0: 境界を除く未知のP,Tの数:  NX    *  NY    *  NZ    = 8000   *
*                                                                   *
*  ●パラメーターファイル（ｉｎ３ｄ．ｓｉｍ）について               *
*    プログラムを実行すると，”ｉｎ３ｄ．ｓｉｍ”というパラメータ   *
*    ファイルを読みに行くので，あらかじめ作成しておく．             *
*                                                                   *
*  "in3d.sim"のリスト                                               *
* 1: U.1......U の計算結果を出力するファイル名                      *
* 2: V.1......V の計算結果を出力するファイル名                      *
* 3: W.1......W の計算結果を出力するファイル名                      *
* 4: P.1......PS(推測した圧力P*)とPD(圧力補正値P') の計算結果を出力 *
*             するファイル名                                        *
* 5: T.1......Θの計算結果を出力するファイル名                      *
* 6: U.0......継続計算用のUの入力データ　　　新規の計算を行う場合は *
* 7: V.0......継続計算用のVの入力データ　　　これらのファイルは不要 *
* 8: W.0......継続計算用のWの入力データ　　　であるが，削除は不可． *
* 9: P.0......継続計算用のPとPDの入力データ                         *
*10: T.0......継続計算用のΘの入力データ                            *
*11: UVT.1....Tecplot(第3部で述べる可視化ソフト)の入力用データ      *
*12: +============================+                                 *
*13: | ITYPE   1===>   isothermal |                                 *
*14: |         2===>nonisothermal |                                 *
*15: +============================+                                 *
*16: -----------------------------------------------------------    *
*17: ITYPE   ICYCLE   NITR    NCYCLE   NCNT                         *
*18: 2       0        100000  10000    1000                         *
*19: -----------------------------------------------------------    *
*20: EPSP    EPSC    OMG     ALPHAP  ALPHAU  ALPHAV  ALPHAW  ALPHAT *
*21: 1.0e-3  1.0e-3  1.0e+0  1.0e-1  0.3e0   0.3e0   0.5e0   0.5e0  *
*22: -----------------------------------------------------------    *
*23: DT        RE       PR       GR                                 *
*24: 1.0e-4    0.0e+0   7.1e-1   1.0e+5                             *
*25: -----------------------------------------------------------    *
*26: DLX    DLY    DLZ    IRELP  METHOD  NMT(1:UPWIND,ELSE:SHISUU)  *
*27: 1.0e0  1.0e0  1.0e0  0      2       2                          *
*                                                                   *
*    ITYPE...1:等温場(浮力項と温度場の計算をしない)  2:非等温場     *
*    ICYCLE.....計算開始のサイクル数（時間T=ICYCLE*DT）             *
*            0なら新規でプログラムにある初期条件にしたがって計算開始*
*            0以外の値なら継続計算                                  *
*    NITR.....(U*,V*,P',T)の各線形システム解法(反復法とクリロフ部分 *
*             空間法)のための最大反復回数．                         *
*    NCYCLE.....計算終了サイクル数                                  *
*    NCNT.......連続の式の収束判定の許容最大反復回数                *
*    EPSP.......収束判定値（速度，圧力，温度の線形システム解法で    *
*               クリロフ部分空間法を含む反復法において使用．）      *
*    EPSC.......連続の式の収束判定値                                *
*    OMG........線形システム解法のうち，(point,line - SOR)を        *
*               用いるときの緩和係数                                *
*    ALPHAP.....圧力の緩和係数 -+  ある時間ステップ(ICYCLE)において *
*    ALPHAU.....Uの緩和係数     |->速度，圧力，温度を修正しながら   *
*    ALPHAV.....Vの緩和係数     |  連続の式を満足する値を求めるとき *
*    ALPHAW.....Wの緩和係数     |  に使用する(NCNT以下であれば収束) *
*    ALPHAT.....温度の緩和係数 -+                                   *
*    DT.....時間刻み                                                *
*    RE.....レイノルズ数, PR.....プラントル数, GR.....グラスホフ数  *
*    DLX.....解析領域の横幅(DLX=NX*DX)                              *
*    DLY.....解析領域の高さ(DLY=NY*DY)                              *
*    DLZ.....解析領域の高さ(DLZ=NZ*DZ)                              *
*            格子幅（等間隔）DX,DY,DZはプログラムの中で求める．     *
*    IRELP....圧力の基準値と解の1次独立性の設定                     *
*             0 : 行わない                                          *
*                   直接法では特異行列の問題に遭遇                  *
*                   数学的には正しくない．                          *
*                   丸め誤差がなければ解は得られない．              *
*             1 : 圧力基準を反映した1次独立な解を求める．           *
*             2 : 圧力基準を設定 P(1,1)=0                           *
*                   1次従属の解のうちの1つを求めて，                *
*                   圧力基準値P(1,1)=0を設定する．                  *
*    METHOD...圧力補正の線形システム解法に用いるアルゴリズム        *
*             すべてバンドマトリックス用に最適化してある            *
*             1: 直接法 -> ガウスの消去法                           *
*                  IRELP=1とする必要がある.(丸め誤差により,IREP=0,2 *
*                  としても解を得られる場合もあるが，数学的には誤り)*
*             2: 反復法1 -> point-SOR 法                            *
*             3: 反復法2 -> line-SOR 法                             *
*                             OMG=1で十分.大きくしすぎると発散.     *
*             4: クリロフ部分空間法1 -> 共役残差法                  *
*             5: クリロフ部分空間法2 -> Bi-CGSTAB法                 *
*    NMT......Simple法における対流・拡散項の離散化方法              *
*             1:風上法 / 2:指数法                                   *
*    パラメータファイルの数値について，                             *
*    FORTRAN プログラム.....できれば倍精度実数で与える"1.0D0,1.0d0" *
*            Cプログラムと共用させて"1.0E0,1.0e0"としても問題はない *
*    C プログラム....."1.0E0,1.0e0"として与える                     *
*    TECPLOT用データを除いた入出力ファイルは書式なし形式で，使用    *
*    するコンパイラーに依存する．コンパイラーに依存しない形式にする *
*    には，容量は増えるが書式付き形式に変更すればよい．             *
*                                                                   *
*  ●圧力の相対性について                                           *
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
*                   NCYCLEになれば計算を終了                        *
*    ITR     -----> 速度，圧力，温度の線形システム解法のための      *
*                   反復回数のカウンタ                              *
*                 （クリロフ部分空間法を含む反復法において使用）    *
*    ITR     -----> 線形システム解法のための反復回数のカウンタ      *
*    ITRCNT  -----> ある時刻における連続の式の評価回数のカウンタ    *
*                   NCNT回以下でEPSCの精度で連続の式を満たせば次の  *
*                   時間ステップへ                                  *
*    IX,IY,IZ    -----> 上の図を参照                                *
*    UO,VO,WO,TO -----> 前の時刻の値                                *
*    U ,V ,W ,T  -----> 新たな圧力補正を用いて計算された値          *
*    PS          -----> 圧力推測値，あるいは補正を反映した圧力      *
*    PD          -----> 圧力補正                                    *
*                                                                   *
*  ●本プログラムは，HSMAC法のプログラムhsmac3d.forや，SMAC法のプロ *
*    グラムsmac3d.forと可能な限り対応するように作成してあり，これら *
*    プログラムとの違いをわかりやすくするため，SIMPLE法において変更 *
*    あるいは新たに付け加えた個所には,"--- SIMPLE ---"としてコメント*
*    を挿入してある．                                               *
*                                                                   *
*  ●SIMPLE法における線形システム解法について                       *
*    SIMPLE法の線形システムの係数行列は速度が影響するため，静止状態 *
*    直後の場合など，速度場がゼロであったりすると，METHOD=1,4,5にお *
*    いて解を得られない場合がある．                                 *
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
*  ●SIMPLE法における収束解について                                 *
*    もともと線形システム解法における収束判定値は，それぞれの解法の *
*    途中で，解がその判定値”以下”となるかどうかを判定するものであ *
*    る．したがって，得られる解はその収束判定値に依存する場合がある.*
*    SIMPLE法ではとくに，線形システム解法で得られる速度場が係数行列 *
*    にも影響するので，線形システム解法の種類によって解が若干異なる *
*    場合がある．                                                   *
*                                                                   *
*  ●ほか                                                           *
*    この計算プログラム例では，計算終了時に，                       *
*      1. 圧力補正の線形システム解法を行った回数の合計              *
*      2. 連続の式を満たすのにかかった回数の合計                    *
*    を画面表示しますので，各解法による違いを確認できます．         *
*                                                                   *
* SIMPLE法については以下の文献を参照                                *
* [1] スハス.V.パタンカー 著, 水谷幸夫・香月正司 共訳               *
*     コンピュータによる熱移動と流れの数値解析, 森北出版 (1985).    *
*                                                                   *
* SIMPLE法を用いたより汎用的なプログラミングについては以下の文献を  *
* 参照                                                              *
* [2] 香月正司・中山　顕 共著                                       *
*     熱流動の数値シミュレーション, 森北出版 (1990)                 *
*                                                                   *
*********************************************************************
      PROGRAM SMPL3D
*********************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NX0=20,NY0=20,NZ0=20,NXY0=400,NEU0=7600,
     $           NEV0=7600,NEW0=7600,NEP0=8000)
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP,EPSC,
     $              ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT
      COMMON / D5 / ICYCLE,ITR,IFLG,IFLGC,IRELP,METHOD,NMT
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
C     --- SIMPLE ---
      COMMON / D8 / NEU,NEV,NEW,NEP
      COMMON / D9 / NITR, NCNT
      COMMON / D10 / RO
      COMMON / ARRAY1 / AU3(0:NX0,0:NY0+1,0:NZ0+1),
     $AV3(0:NX0+1,0:NY0,0:NZ0+1),AW3(0:NX0+1,0:NY0+1,0:NZ0),
     $ATP1(NEP0),ATP2(NEP0),ATP3(NEP0),ATP4(NEP0),ATP5(NEP0),
     $ATP6(NEP0),ATP7(NEP0),BP(NEP0),
     $ U(0:NX0,  0:NY0+1,0:NZ0+1),  V(0:NX0+1,0:NY0,  0:NZ0+1),
     $ W(0:NX0+1,0:NY0+1,0:NZ0  ),
     $US(0:NX0,  0:NY0+1,0:NZ0+1), VS(0:NX0+1,0:NY0,  0:NZ0+1),
     $WS(0:NX0+1,0:NY0+1,0:NZ0  ),
     $PS(0:NX0+1,0:NY0+1,0:NZ0+1), PD(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ T(0:NX0+1,0:NY0+1,0:NZ0+1),
     $XX(NEP0)
      COMMON / ARRAY2 /
     $UO(0:NX0  ,0:NY0+1,0:NZ0+1),VO(0:NX0+1,0:NY0,  0:NZ0+1),
     $WO(0:NX0+1,0:NY0+1,0:NZ0  ),TO(0:NX0+1,0:NY0+1,0:NZ0+1)
      CHARACTER FNAME(20)*20
*********************************************************************
* コンパイラーによっては，以下の2つのDIMENSION文は，メインプログラム*
* で必ずしも宣言されていなくてもよく，これらの配列を用いるSUBROUTINE*
* においてのみ宣言すればよい．この場合，これに対応して，CALL文や    *
* SUBROUTINE文の引数にあるこれらの配列を削除する．                  *
*********************************************************************
*     作業用配列
      DIMENSION W1(NEP0),W2(NEP0),W3(NEP0),W4(NEP0),W5(NEP0),W6(NEP0),
     $          W7(NEP0),W8(NEP0),W9(NEP0)
*     バンドガウス消去法のための係数行列配列
      DIMENSION AGB(-NXY0:NXY0,NEP0)
*********************************************************************
*
*x方向の格子分割数
      NX = NX0
*y方向の格子分割数
      NY = NY0
*z方向の格子分割数
      NZ = NZ0
* NX * NY
C     NXY=NXY0
*方程式の数( U* に関する未知数の数) --- SIMPLE ---
      NEU = NEU0
*方程式の数( V* に関する未知数の数) --- SIMPLE ---
      NEV = NEV0
*方程式の数( W* に関する未知数の数) --- SIMPLE ---
      NEW = NEW0
*方程式の数(圧力補正に関する未知数の数) --- SIMPLE ---
      NEP = NEP0
*圧力計算の線形システム解法の総反復回数のための変数
      ITRALL = 0
      ICNTAL = 0
*パラメータファイルのオープン
      OPEN (10,FILE='IN3D.SIM',STATUS='OLD')
*出力ファイル名の読み込み
      DO 20 I = 1,11
        READ (10,'(A20)') FNAME(I)
   20 CONTINUE
*書式なし形式はコンパイラーに依存するので注意
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
*パラメータファイルからパラメータを読み込む
  700 CONTINUE
      READ (10,*,ERR=700) ITYPE,ICYCLE,NITR,NCYCLE,NCNT
  710 CONTINUE
      READ (10,*,ERR=710) EPSP,EPSC,OMG,
     $                    ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT
  720 CONTINUE
      READ (10,*,ERR=720) DT,RE,PR,GR
  730 CONTINUE
      READ (10,*,ERR=730) DLX,DLY,DLZ,IRELP,METHOD,NMT
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
      DX = DLX / FLOAT(NX)
*y方向の格子幅
      DY = DLY / FLOAT(NY)
*z方向の格子幅
      DZ = DLZ / FLOAT(NZ)
*運動方程式中の拡散項の係数(ここではPr)
      VIS = PR
*浮力項の係数(ここでは Gr * Pr**2)
      BUO = GR * ( PR**2 )
*エネルギー方程式中の拡散項の係数(ここでは1)
      ALP = 1.0D0
*等温場なら浮力項の係数はゼロに設定
      IF (ITYPE.EQ.1) BUO = 0.0D0
*文献[1]にしたがった離散化を考慮してROを定義
*本書では無次元化をおこなってから計算しており，RO=1とする
      RO = 1.0D0
*初期値の設定
      CALL CINITI (AGB,W1,W2,W3,W4,W5,W6,W7,W8,W9)
*時間進行のための戻り点
  750 CONTINUE
*時間進行
      CALL ADV
*連続の式の収束判定の反復回数を表す変数
      ITRCNT = 0
*速度・圧力・温度補正のための戻り点
  760 CONTINUE
      ITRCNT = ITRCNT + 1
*速度，圧力，温度の線形システム解法が収束したかどうかのパラメータ(IFLG)を初期化
*これは，クリロフ部分空間法を含む反復法で用いる
*IFLG -> 0:収束 1:発散(設定された許容回数NITR以下で解が得られない)
      IFLG  = 0
*速度，圧力，温度を同時に補正しながら連続の式を満たすまで計算するが，
*この際，連続の式を満たすかどうかを判定する反復回数を意味する
*パラメータ(IFLGC)を初期化
*IFLGC-> 0:連続の式を満たす 1:発散(許容回数NCNT以下で解が得られない)
      IFLGC = 0
C     --- 注意 ---
C     本計算においては，ICYCLE=0のときの初期条件は速度場と圧力場をゼロと
C     しているので，これらに関わる線形システムを解くのが困難となるため，
C     ICYCLE=1のときだけは温度の計算にジャンプするようにする．
C     計算条件や問題に応じて適宜削除あるいは変更
      IF (ICYCLE.EQ.1.AND.ITRCNT.EQ.1) GOTO 770
*速度場の計算
      CALL CALVEL (AGB,W1,W2,W3,W4,W5,W6,W7,W8,W9)
*     速度計算の線形システム解法が収束していないとき
      IF ( IFLG.EQ.1 ) THEN
        WRITE (6,*) ' NOT CONVERGE ! '
C       データを出力して強制終了
        CALL PROUT
        GO TO 900
      END IF
*圧力補正の計算と，圧力場と速度場の更新
      CALL PRESS (AGB,W1,W2,W3,W4,W5,W6,W7,W8,W9)
*     圧力計算の線形システム解法の総反復回数
      ITRALL = ITRALL + ITR
C     速度場と圧力場がゼロとなるときのスキップ先
  770 CONTINUE
*     圧力計算の線形システム解法が収束していないとき
      IF ( IFLG.EQ.1 ) THEN
        WRITE (6,*) ' NOT CONVERGE ! '
C       データを出力して強制終了
        CALL PROUT
        GO TO 900
      END IF
*     非等温場計算の場合
      IF ( ITYPE.EQ.2 ) THEN
*       温度場を計算
        CALL CALTEM (AGB,W1,W2,W3,W4,W5,W6,W7,W8,W9)
*       温度計算の線形システム解法が収束していないとき
        IF ( IFLG.EQ.1 ) THEN
          WRITE (6,*) ' NOT CONVERGE ! '
C         データを出力して強制終了
          CALL PROUT
          GO TO 900
        END IF
      END IF
*     速度，圧力，温度場の計算が収束して連続の式をEPSP以下で満足していないとき(IFLGC=1)
      IF ( IFLGC.EQ.1 ) THEN
*       予め設定された許容回数NCNTに達していないときは再計算
        IF (ITRCNT.LT.NCNT) THEN
          GO TO 760
*       予め設定された許容回数NCNTに達したらデータを出力して強制終了
        ELSE
          WRITE (6,*) ' NOT CONVERGE ! '
C         データを出力して強制終了
          CALL PROUT
          GO TO 900
        END IF
      END IF
*     速度，圧力，温度場の計算が収束して連続の式をEPSP以下で満足し(IFLGC=0)
*     時間進行カウンタ(ICYCLE)がNCYCLEより小さい時
      ICNTAL = ICNTAL + ITRCNT
      IF ( ICYCLE. LT. NCYCLE ) THEN
        GO TO 750
*     時間進行カウンタがNCYCLEになったら計算終了
      ELSE
        CALL PROUT
      END IF
  900 CONTINUE
*Tecplot用データの出力
      CALL TECPLT(FNAME(11))
      CLOSE (10)
      CLOSE (11)
      CLOSE (12)
      CLOSE (13)
      CLOSE (14)
      CLOSE (15)
*
*     総反復回数を出力
      WRITE (6,*) 'Total Iteration Number'
      WRITE (6,*) 'To solve Continuity Eq.     ',ICNTAL
      WRITE (6,*) 'To solve Linear System of PD',ITRALL
*
      STOP
      END
*********************************************************************
*                        初期設定
*********************************************************************
      SUBROUTINE CINITI (AGB,W1,W2,W3,W4,W5,W6,W7,W8,W9)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NX0=20,NY0=20,NZ0=20,NXY0=400,NEU0=7600,
     $           NEV0=7600,NEW0=7600,NEP0=8000)
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP,EPSC,
     $              ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT
      COMMON / D5 / ICYCLE,ITR,IFLG,IFLGC,IRELP,METHOD,NMT
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
C     --- SIMPLE ---
      COMMON / D8 / NEU,NEV,NEW,NEP
      COMMON / D9 / NITR, NCNT
      COMMON / D10 / RO
      COMMON / ARRAY1 / AU3(0:NX0,0:NY0+1,0:NZ0+1),
     $AV3(0:NX0+1,0:NY0,0:NZ0+1),AW3(0:NX0+1,0:NY0+1,0:NZ0),
     $ATP1(NEP0),ATP2(NEP0),ATP3(NEP0),ATP4(NEP0),ATP5(NEP0),
     $ATP6(NEP0),ATP7(NEP0),BP(NEP0),
     $ U(0:NX0,  0:NY0+1,0:NZ0+1),  V(0:NX0+1,0:NY0,  0:NZ0+1),
     $ W(0:NX0+1,0:NY0+1,0:NZ0  ),
     $US(0:NX0,  0:NY0+1,0:NZ0+1), VS(0:NX0+1,0:NY0,  0:NZ0+1),
     $WS(0:NX0+1,0:NY0+1,0:NZ0  ),
     $PS(0:NX0+1,0:NY0+1,0:NZ0+1), PD(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ T(0:NX0+1,0:NY0+1,0:NZ0+1),
     $XX(NEP0)
      COMMON / ARRAY2 /
     $UO(0:NX0  ,0:NY0+1,0:NZ0+1),VO(0:NX0+1,0:NY0,  0:NZ0+1),
     $WO(0:NX0+1,0:NY0+1,0:NZ0  ),TO(0:NX0+1,0:NY0+1,0:NZ0+1)
*     作業用配列
      DIMENSION W1(NEP0),W2(NEP0),W3(NEP0),W4(NEP0),W5(NEP0),W6(NEP0),
     $          W7(NEP0),W8(NEP0),W9(NEP0)
*     バンドガウス消去法のための係数行列配列
      DIMENSION AGB(-NXY0:NXY0,NEP0)
*
*新規計算の場合
      IF ( ICYCLE. EQ. 0 ) THEN
*       U の初期値設定
        DO 10 IX = 0,NX
          DO 20 IY = 0,NY+1
            DO 30 IZ = 0,NZ+1
              U(IX,IY,IZ) = 0.0D0
   30       CONTINUE
   20     CONTINUE
   10   CONTINUE
*       V の初期値設定
        DO 40 IX = 0,NX+1
          DO 50 IY = 0,NY
            DO 60 IZ = 0,NZ+1
              V(IX,IY,IZ) = 0.0D0
   60     CONTINUE
   50     CONTINUE
   40   CONTINUE
*       W の初期値設定
        DO 70 IX = 0,NX+1
          DO 80 IY = 0,NY+1
            DO 90 IZ = 0,NZ
              W(IX,IY,IZ) = 0.0D0
   90     CONTINUE
   80     CONTINUE
   70   CONTINUE
*       P* と P' の初期値設定
        DO 100 IX = 0,NX+1
          DO 110 IY = 0,NY+1
            DO 120 IZ = 0,NZ+1
C             --- SIMPLE ---
              PS(IX,IY,IZ) = 0.0D0
              PD(IX,IY,IZ) = 0.0D0
  120     CONTINUE
  110     CONTINUE
  100   CONTINUE
*----------------------------------------------------------------------*
*（注意）浮力項の計算で温度の配列を使用しているので等温場でもT=0として *
* 初期条件だけは設定する必要がある．ゼロ以外の値を入れると浮力項が計算 *
* される可能性があるので注意．                                         *
*----------------------------------------------------------------------*
*       Tの初期値設定(領域内は高温(+0.5)と低温(-0.5)の中間温度)
        DO 130 IX = 0,NX+1
          DO 140 IY = 0,NY+1
            DO 150 IZ = 0,NZ+1
              T(IX,IY,IZ) = 0.0D0
  150       CONTINUE
  140     CONTINUE
  130   CONTINUE
*       Tの境界：右面（冷却）T=-0.5
        DO 160 IY = 0,NY+1
          DO 170 IZ = 0,NZ+1
            T(NX+1,IY,IZ) = 2.0D0 * ( -0.5D0 ) - T(NX,IY,IZ)
  170     CONTINUE
  160   CONTINUE
*       Tの境界：左面（加熱）T=+0.5
        DO 180 IY = 0,NY+1
          DO 190 IZ = 0,NZ+1
            T(0,IY,IZ) = 2.0D0 * ( +0.5D0 ) - T(1,IY,IZ)
  190     CONTINUE
  180   CONTINUE
*       Tの境界：上面（断熱）
        DO 200 IX = 1,NX
          DO 210 IZ = 1,NZ+1
            T(IX,NY+1,IZ) = T(IX,NY,IZ)
  210     CONTINUE
  200   CONTINUE
*       Tの境界：下面（断熱）
        DO 220 IX = 1,NX
          DO 230 IZ = 1,NZ+1
            T(IX,0,IZ) = T(IX,1,IZ)
  230     CONTINUE
  220   CONTINUE
*       Tの境界：後面（断熱）
        DO 240 IX = 1,NX
          DO 250 IY = 1,NY
            T(IX,IY,NZ+1) = T(IX,IY,NZ)
  250     CONTINUE
  240   CONTINUE
*       Tの境界：前面（断熱）
        DO 260 IX = 1,NX
          DO 270 IY = 1,NY
            T(IX,IY,0) = T(IX,IY,1)
  270     CONTINUE
  260   CONTINUE
*継続計算（すでにある計算結果からスタート）の場合
      ELSE
*       Uデータファイルからの読み込み[Unit No.=16](書式なし形式)
        READ (16) U
*       Vデータファイルからの読み込み[Unit No.=17](書式なし形式)
        READ (17) V
*       Wデータファイルからの読み込み[Unit No.=18](書式なし形式)
        READ (18) W
*       Pデータファイルからの読み込み[Unit No.=19](書式なし形式)
C       --- SIMPLE ---
        READ (19) PS,PD
*---------------------------------------------------------------------*
*    (注意) 等温場の計算でもT(=0)のファイルを読み込む必要がある       *
*---------------------------------------------------------------------*
*       Tデータファイルからの読み込み[Unit No.=20](書式なし形式)
        READ (20) T
        CLOSE (16)
        CLOSE (17)
        CLOSE (18)
        CLOSE (19)
        CLOSE (20)
      END IF
*線形システム解法に関わる配列も初回のみ初期化
*U*, V*, W*, P', Tに関わる線形システム解法の解が配列XXに入る
      DO 280 I=1,NEP
        ATP1(I)=0.0D0
        ATP2(I)=0.0D0
        ATP3(I)=0.0D0
        ATP4(I)=0.0D0
        ATP5(I)=0.0D0
        ATP6(I)=0.0D0
        ATP7(I)=0.0D0
        BP(I)=0.0D0
* METHOD=1,4,5のとき適宜修正(--- SIMPLE ---)
* XX をここでのみ初期化するか線形システム解法の度に直前で初期化するなどする
        XX(I)=0.0D0
        W1(I)=0.0D0
        W2(I)=0.0D0
        W3(I)=0.0D0
        W4(I)=0.0D0
        W5(I)=0.0D0
        W6(I)=0.0D0
        W7(I)=0.0D0
        W8(I)=0.0D0
        W9(I)=0.0D0
        DO 290 II = -NXY0,NXY0
          AGB(II,I)=0.0D0
  290   CONTINUE
  280 CONTINUE
      RETURN
      END
*********************************************************************
*                         時間進行
*********************************************************************
      SUBROUTINE ADV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NX0=20,NY0=20,NZ0=20,NXY0=400,NEU0=7600,
     $           NEV0=7600,NEW0=7600,NEP0=8000)
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP,EPSC,
     $              ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT
      COMMON / D5 / ICYCLE,ITR,IFLG,IFLGC,IRELP,METHOD,NMT
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
C     --- SIMPLE ---
      COMMON / D8 / NEU,NEV,NEW,NEP
      COMMON / D9 / NITR, NCNT
      COMMON / D10 / RO
      COMMON / ARRAY1 / AU3(0:NX0,0:NY0+1,0:NZ0+1),
     $AV3(0:NX0+1,0:NY0,0:NZ0+1),AW3(0:NX0+1,0:NY0+1,0:NZ0),
     $ATP1(NEP0),ATP2(NEP0),ATP3(NEP0),ATP4(NEP0),ATP5(NEP0),
     $ATP6(NEP0),ATP7(NEP0),BP(NEP0),
     $ U(0:NX0,  0:NY0+1,0:NZ0+1),  V(0:NX0+1,0:NY0,  0:NZ0+1),
     $ W(0:NX0+1,0:NY0+1,0:NZ0  ),
     $US(0:NX0,  0:NY0+1,0:NZ0+1), VS(0:NX0+1,0:NY0,  0:NZ0+1),
     $WS(0:NX0+1,0:NY0+1,0:NZ0  ),
     $PS(0:NX0+1,0:NY0+1,0:NZ0+1), PD(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ T(0:NX0+1,0:NY0+1,0:NZ0+1),
     $XX(NEP0)
      COMMON / ARRAY2 /
     $UO(0:NX0  ,0:NY0+1,0:NZ0+1),VO(0:NX0+1,0:NY0,  0:NZ0+1),
     $WO(0:NX0+1,0:NY0+1,0:NZ0  ),TO(0:NX0+1,0:NY0+1,0:NZ0+1)
*
      TIME   = DT*FLOAT(ICYCLE)
      ICYCLE = ICYCLE + 1
*時間進行カウンタ(ICYCLE)を100回毎に表示
      IF (MOD(ICYCLE,100).EQ.0) THEN
        WRITE (6,2000) ICYCLE
 2000   FORMAT ('  CYC = ',2I8)
      END IF
*--------------------------------------------------------------------
* U  -> UO : 必要なら入れ替える前にUとUOから変動量を求める
* U  : 前の時間ステップにおいて最終的に得られた値
*      速度・圧力補正の繰り返しで U が更新されるので，
* UO : 前の時間ステップの値を保存しておく
*--------------------------------------------------------------------
      DO 10 IX = 0,NX
        DO 20 IY = 0,NY+1
          DO 30 IZ = 0,NZ+1
            UO(IX,IY,IZ) = U(IX,IY,IZ)
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
*--------------------------------------------------------------------
* V  -> VO : 必要なら入れ替える前にVとVOから変動量を求める
* V  : 前の時間ステップにおいて最終的に得られた値
*      速度・圧力補正の繰り返しで V が更新されるので，
* VO : 前の時間ステップの値を保存しておく
*--------------------------------------------------------------------
      DO 40 IX = 0,NX+1
        DO 50 IY = 0,NY
          DO 60 IZ = 0,NZ+1
            VO(IX,IY,IZ) = V(IX,IY,IZ)
   60     CONTINUE
   50   CONTINUE
   40 CONTINUE
*--------------------------------------------------------------------
* W  -> WO : 必要なら入れ替える前にWとWOから変動量を求める
* W  : 前の時間ステップにおいて最終的に得られた値
*      速度・圧力補正の繰り返しで W が更新されるので，
* WO : 前の時間ステップの値を保存しておく
*--------------------------------------------------------------------
      DO 70 IX = 0,NX+1
        DO 80 IY = 0,NY+1
          DO 90 IZ = 0,NZ
            WO(IX,IY,IZ) = W(IX,IY,IZ)
   90     CONTINUE
   80   CONTINUE
   70 CONTINUE
*--------------------------------------------------------------------
* T  -> TO : 必要なら入れ替える前にTとTOから変動量を求める
* T  : 前の時間ステップにおいて最終的に得られた値
*      速度・圧力補正の繰り返しで T が更新されるので，
* TO : 前の時間ステップの値を保存しておく
*--------------------------------------------------------------------
      DO 100 IX = 0,NX+1
        DO 110 IY = 0,NY+1
          DO 120 IZ = 0,NZ+1
            TO(IX,IY,IZ) = T(IX,IY,IZ)
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE
      RETURN
      END
*********************************************************************
*                 仮の速度場 (U*, V*, W*) の計算
*********************************************************************
      SUBROUTINE CALVEL (AGB,W1,W2,W3,W4,W5,W6,W7,W8,W9)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NX0=20,NY0=20,NZ0=20,NXY0=400,NEU0=7600,
     $           NEV0=7600,NEW0=7600,NEP0=8000)
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP,EPSC,
     $              ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT
      COMMON / D5 / ICYCLE,ITR,IFLG,IFLGC,IRELP,METHOD,NMT
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
C     --- SIMPLE ---
      COMMON / D8 / NEU,NEV,NEW,NEP
      COMMON / D9 / NITR, NCNT
      COMMON / D10 / RO
      COMMON / ARRAY1 / AU3(0:NX0,0:NY0+1,0:NZ0+1),
     $AV3(0:NX0+1,0:NY0,0:NZ0+1),AW3(0:NX0+1,0:NY0+1,0:NZ0),
     $ATP1(NEP0),ATP2(NEP0),ATP3(NEP0),ATP4(NEP0),ATP5(NEP0),
     $ATP6(NEP0),ATP7(NEP0),BP(NEP0),
     $ U(0:NX0,  0:NY0+1,0:NZ0+1),  V(0:NX0+1,0:NY0,  0:NZ0+1),
     $ W(0:NX0+1,0:NY0+1,0:NZ0  ),
     $US(0:NX0,  0:NY0+1,0:NZ0+1), VS(0:NX0+1,0:NY0,  0:NZ0+1),
     $WS(0:NX0+1,0:NY0+1,0:NZ0  ),
     $PS(0:NX0+1,0:NY0+1,0:NZ0+1), PD(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ T(0:NX0+1,0:NY0+1,0:NZ0+1),
     $XX(NEP0)
      COMMON / ARRAY2 /
     $UO(0:NX0  ,0:NY0+1,0:NZ0+1),VO(0:NX0+1,0:NY0,  0:NZ0+1),
     $WO(0:NX0+1,0:NY0+1,0:NZ0  ),TO(0:NX0+1,0:NY0+1,0:NZ0+1)
*     作業用配列
      DIMENSION W1(NEP0),W2(NEP0),W3(NEP0),W4(NEP0),W5(NEP0),W6(NEP0),
     $          W7(NEP0),W8(NEP0),W9(NEP0)
*     バンドガウス消去法のための係数行列配列
      DIMENSION AGB(-NXY0:NXY0,NEP0)
*
*-------------------------------------------------*
*            U* : US(IX,IY,IZ) の計算             *
*-------------------------------------------------*
* I:線形システム行列用の通し番号を意味する変数
* 3次元配列(IX,IY,IZ):離散化用 -> 1次元配列(I):線形システム用
      I = 1
      DO 10 IZ = 1,NZ
        DO 20 IX = 1,NX-1
          DO 30 IY = 1,NY
*         Uに関するDの計算
          DU5 = VIS * DY*DZ / DX
          DU1 = VIS * DY*DZ / DX
          DU4 = VIS * DX*DZ / DY
          DU2 = VIS * DX*DZ / DY
          DU7 = VIS * DX*DY / DZ
          DU6 = VIS * DX*DY / DZ
*         Uに関するFの計算
*         速度・圧力・温度場の繰り返し計算となっても線形システムの係数となる
*         ここの速度は更新されないように(U,V,W)には前の時刻の収束値(UO,VO,WO)を用いる
          UU5 = ( UO(IX  ,IY  ,IZ  ) + UO(IX+1,IY  ,IZ  ) ) / 2.0D0
          UU1 = ( UO(IX-1,IY  ,IZ  ) + UO(IX  ,IY  ,IZ  ) ) / 2.0D0
          VU4 = ( VO(IX  ,IY  ,IZ  ) + VO(IX+1,IY  ,IZ  ) ) / 2.0D0
          VU2 = ( VO(IX  ,IY-1,IZ  ) + VO(IX+1,IY-1,IZ  ) ) / 2.0D0
          WU7 = ( WO(IX  ,IY  ,IZ  ) + WO(IX+1,IY  ,IZ  ) ) / 2.0D0
          WU6 = ( WO(IX  ,IY  ,IZ-1) + WO(IX+1,IY  ,IZ-1) ) / 2.0D0
          FU5 = RO * UU5 * DY * DZ
          FU1 = RO * UU1 * DY * DZ
          FU4 = RO * VU4 * DX * DZ
          FU2 = RO * VU2 * DX * DZ
          FU7 = RO * WU7 * DX * DY
          FU6 = RO * WU6 * DX * DY
*         Uに関するPeの計算
          PEU5 = FU5 / DU5
          PEU1 = FU1 / DU1
          PEU4 = FU4 / DU4
          PEU2 = FU2 / DU2
          PEU7 = FU7 / DU7
          PEU6 = FU6 / DU6
*         関数A(|Pe|)の計算
*         NMT=1なら 風上法
          IF ( NMT.EQ.1 ) THEN
            APEU5 = 1.0D0
            APEU1 = 1.0D0
            APEU4 = 1.0D0
            APEU2 = 1.0D0
            APEU7 = 1.0D0
            APEU6 = 1.0D0
*         NMT=2なら べき乗法
          ELSE
            APEU5 = DMAX1( 0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEU5))**5 ) )
            APEU1 = DMAX1( 0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEU1))**5 ) )
            APEU4 = DMAX1( 0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEU4))**5 ) )
            APEU2 = DMAX1( 0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEU2))**5 ) )
            APEU7 = DMAX1( 0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEU7))**5 ) )
            APEU6 = DMAX1( 0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEU6))**5 ) )
          END IF
*         U* の線形システムの係数の計算
          AU5 = ( DU5*APEU5 + DMAX1(-FU5,0.0D0) )
          AU1 = ( DU1*APEU1 + DMAX1( FU1,0.0D0) )
          AU4 = ( DU4*APEU4 + DMAX1(-FU4,0.0D0) )
          AU2 = ( DU2*APEU2 + DMAX1( FU2,0.0D0) )
          AU7 = ( DU7*APEU7 + DMAX1(-FU7,0.0D0) )
          AU6 = ( DU6*APEU6 + DMAX1( FU6,0.0D0) )
*         生成項
          SCU = 0.0D0
          SPU = 0.0D0
*         U* の線形システムのBの計算
          BP(I) = SCU*DX*DY*DZ+( PS(IX,IY,IZ)-PS(IX+1,IY,IZ) )*DY*DZ
     $          + RO*DX*DY*DZ/DT*UO(IX,IY,IZ)
*境界の処理
* 1. 境界のU_BND についても変数とみなしてAU3(IX,IY,IZ)を求める
* 2. その後，境界の項 AU(1,2,4,5,6,7)*U_BND をBに吸収させる
* 3. 2に対応して(吸収の後)，境界のAU(1,2,4,5,6,7)をゼロとしてリンクを断つ
          AU3(IX,IY,IZ) = AU1+AU2+AU4+AU5+AU6+AU7 + RO*DX*DY*DZ/DT
     $                  - SPU*DX*DY*DZ
*U（右面）
          IF ( IX.EQ.NX-1 ) THEN
            BP(I) = BP(I)+AU5*U(NX,IY,IZ)
            AU5 = 0.0D0
          END IF
*U（左面）
          IF ( IX.EQ.1 ) THEN
            BP(I) = BP(I)+AU1*U(0,IY,IZ)
            AU1 = 0.0D0
          END IF
*U（上面）
          IF ( IY.EQ.NY ) THEN
            BP(I) = BP(I) + AU4*U(IX,NY+1,IZ)
            AU4 = 0.0D0
          END IF
*U（下面）
          IF ( IY.EQ.1 ) THEN
            BP(I) = BP(I) + AU2*U(IX,0,IZ)
            AU2 = 0.0D0
          END IF
*U（後面）
          IF ( IZ.EQ.NZ ) THEN
            BP(I) = BP(I) + AU7*U(IX,IY,NZ+1)
            AU7 = 0.0D0
          END IF
*U（前面）
          IF ( IZ.EQ.1 ) THEN
            BP(I) = BP(I) + AU6*U(IX,IY,0)
            AU6 = 0.0D0
          END IF
* 線形システム計算用に係数行列を計算
* ALPHAU : U* 計算のための緩和係数
          ATP3(I) = AU3(IX,IY,IZ) / ALPHAU
          ATP5(I) = - AU5
          ATP1(I) = - AU1
          ATP4(I) = - AU4
          ATP2(I) = - AU2
          ATP7(I) = - AU7
          ATP6(I) = - AU6
          BP(I) = BP(I) + (1.0D0-ALPHAU)*ATP3(I)*UO(IX,IY,IZ)
* 線形システム解法前のための初期化(任意) (--- SIMPLE ---)
* METHOD=1,4,5のとき適宜修正
          XX(I) = U(IX,IY,IZ)
          I = I + 1
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
*
*サブルーチンに一般性を持たせるため，NX,NY,NEを引数として渡す
*COMMON文で定義されている値でもあるので，そのまま渡せない！
      NNX = NX-1
      NNY = NY
      NNZ = NZ
      NNE = NEU0
      NNXY = NNX*NNY
*以下の線形システム解法の後，配列XXにUS(U*)が入る
C   1. 直接法 : バンドマトリックスによるガウスの消去法
      IF (METHOD.EQ.1)
     $  CALL GB      (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                AGB)
C   2. 反復法1 : point-SOR 法
      IF (METHOD.EQ.2)
     $  CALL PSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY)
C   3. 反復法2 : line-SOR 法 : OMG=1で十分．大きくしすぎると発散する
      IF (METHOD.EQ.3)
     $  CALL LSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6,W7,W8,W9)
C   4. クリロフ部分空間法1 : 共役残差法
      IF (METHOD.EQ.4)
     $  CALL CRB     (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6)
C   5. クリロフ部分空間法2 : BiCGSTAB
      IF (METHOD.EQ.5)
     $  CALL BICGB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6,W7)
C   METHOD=4,5のときの探索ベクトル計算時のゼロ除算対策
      IF (IFLG.EQ.2)
     $  CALL PSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY)
*線形システム解法で得られた1次元配列の解を3次元配列に置き換える
      DO 40 IX = 1,NX-1
        DO 50 IY = 1,NY
          DO 60 IZ = 1,NZ
            K = IY + (IX-1)*NNY + (IZ-1)*(NNX*NNY)
            US(IX,IY,IZ) = XX(K)
   60     CONTINUE
   50   CONTINUE
   40 CONTINUE
*U*の境界条件の処理
      CALL UBND (US)
*
*-------------------------------------------------*
*            V* : VS(IX,IY,IZ) の計算             *
*-------------------------------------------------*
* I:線形システム行列用の通し番号を意味する変数
* 3次元配列(IX,IY,IZ):離散化用 -> 1次元配列(I):線形システム用
      I = 1
      DO 70 IZ = 1,NZ
        DO 80 IX = 1,NX
          DO 90 IY = 1,NY-1
*         Vに関するDの計算
          DV5 = VIS * DY*DZ / DX
          DV1 = VIS * DY*DZ / DX
          DV4 = VIS * DX*DZ / DY
          DV2 = VIS * DX*DZ / DY
          DV7 = VIS * DX*DY / DZ
          DV6 = VIS * DX*DY / DZ
*         Vに関するFの計算
*         速度・圧力・温度場の繰り返し計算となっても線形システムの係数となる
*         ここの速度は更新されないように(U,V)には前の時刻の収束値(UO,VO)を用いる
          UV5 = ( UO(IX  ,IY  ,IZ  ) + UO(IX  ,IY+1,IZ  ) ) / 2.0D0
          UV1 = ( UO(IX-1,IY  ,IZ  ) + UO(IX-1,IY+1,IZ  ) ) / 2.0D0
          VV4 = ( VO(IX  ,IY  ,IZ  ) + VO(IX  ,IY+1,IZ  ) ) / 2.0D0
          VV2 = ( VO(IX  ,IY-1,IZ  ) + VO(IX  ,IY  ,IZ  ) ) / 2.0D0
          WV7 = ( WO(IX  ,IY  ,IZ  ) + WO(IX  ,IY+1,IZ  ) ) / 2.0D0
          WV6 = ( WO(IX  ,IY  ,IZ-1) + WO(IX  ,IY+1,IZ-1) ) / 2.0D0
          FV5 = RO * UV5 * DY * DZ
          FV1 = RO * UV1 * DY * DZ
          FV4 = RO * VV4 * DX * DZ
          FV2 = RO * VV2 * DX * DZ
          FV7 = RO * WV7 * DX * DY
          FV6 = RO * WV6 * DX * DY
*         Vに関するPeの計算
          PEV5 = FV5 / DV5
          PEV1 = FV1 / DV1
          PEV4 = FV4 / DV4
          PEV2 = FV2 / DV2
          PEV7 = FV7 / DV7
          PEV6 = FV6 / DV6
*         関数A(|Pe|)の計算
*         NMT=1なら 風上法
          IF ( NMT.EQ.1 ) THEN
            APEV5 = 1.0D0
            APEV1 = 1.0D0
            APEV4 = 1.0D0
            APEV2 = 1.0D0
            APEV7 = 1.0D0
            APEV6 = 1.0D0
*         NMT=2なら べき乗法
          ELSE
            APEV5 = DMAX1(0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEV5))**5 ) )
            APEV1 = DMAX1(0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEV1))**5 ) )
            APEV4 = DMAX1(0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEV4))**5 ) )
            APEV2 = DMAX1(0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEV2))**5 ) )
            APEV7 = DMAX1(0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEV7))**5 ) )
            APEV6 = DMAX1(0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEV6))**5 ) )
          END IF
*         V* の線形システムの係数の計算
          AV5 = ( DV5*APEV5 + DMAX1(-FV5,0.0D0) )
          AV1 = ( DV1*APEV1 + DMAX1( FV1,0.0D0) )
          AV4 = ( DV4*APEV4 + DMAX1(-FV4,0.0D0) )
          AV2 = ( DV2*APEV2 + DMAX1( FV2,0.0D0) )
          AV7 = ( DV7*APEV7 + DMAX1(-FV7,0.0D0) )
          AV6 = ( DV6*APEV6 + DMAX1( FV6,0.0D0) )
*         速度・圧力・温度場の繰り返し計算となった場合，更新されたT
*         により浮力も更新される
          TV = ( T(IX,IY,IZ) + T(IX,IY+1,IZ) ) / 2.0D0
*         生成項(浮力項)の計算
          SCV = BUO * TV
          SPV = 0.0D0
*         V* の線形システムのBの計算
          BP(I) = SCV*DX*DY*DZ+( PS(IX,IY,IZ)-PS(IX,IY+1,IZ) )*DX*DZ
     $          + RO*DX*DY*DZ/DT*VO(IX,IY,IZ)

*境界の処理
* 1. 境界のV_BND についても変数とみなしてAV3(IX,IY,IZ)を求める
* 2. その後，境界の項 AV(1,2,4,5,6,7)*V_BND をBに吸収させる
* 3. 2に対応して(吸収の後)，境界のAV(1,2,4,5,6,7)をゼロとしてリンクを断つ
          AV3(IX,IY,IZ) = AV1+AV2+AV4+AV5+AV6+AV7 + RO*DX*DY*DZ/DT
     $                  - SPV*DX*DY
*V（右面）
          IF ( IX.EQ.NX ) THEN
            BP(I) = BP(I) + AV5 * V(NX+1,IY,IZ)
            AV5 = 0.0D0
          END IF
*V（左面）
          IF ( IX.EQ.1 ) THEN
            BP(I) = BP(I) + AV1 * V(0,IY,IZ)
            AV1 = 0.0D0
          END IF
*V（上面）
          IF ( IY.EQ.NY-1 ) THEN
            BP(I) = BP(I) + AV4*V(IX,NY,IZ)
            AV4 = 0.0D0
          END IF
*V（下面）
          IF ( IY.EQ.1 ) THEN
            BP(I) = BP(I) + AV2*V(IX,0,IZ)
            AV2 = 0.0D0
          END IF
*V（後面）
          IF ( IZ.EQ.NZ ) THEN
            BP(I) = BP(I) + AV7*V(IX,IY,NZ+1)
            AV7 = 0.0D0
          END IF
*V（前面）
          IF ( IZ.EQ.1 ) THEN
            BP(I) = BP(I) + AV6*V(IX,IY,0)
            AV6 = 0.0D0
          END IF
* 線形システム計算用に係数行列を計算
* ALPHAV : V* 計算のための緩和係数
          ATP3(I) = AV3(IX,IY,IZ) / ALPHAV
          ATP5(I) = - AV5
          ATP1(I) = - AV1
          ATP4(I) = - AV4
          ATP2(I) = - AV2
          ATP7(I) = - AV7
          ATP6(I) = - AV6
          BP(I) = BP(I) + (1.0D0-ALPHAV)*ATP3(I)*VO(IX,IY,IZ)
* 線形システム解法前のための初期化(任意) (--- SIMPLE ---)
* METHOD=1,4,5のとき適宜修正
          XX(I) = V(IX,IY,IZ)
          I = I + 1
   90     CONTINUE
   80   CONTINUE
   70 CONTINUE
*
*サブルーチンに一般性を持たせるため，NX,NY,NEを引数として渡す
*COMMON文で定義されている値でもあるので，そのまま渡せない！
      NNX = NX
      NNY = NY-1
      NNZ = NZ
      NNE = NEV
      NNXY = NNX*NNY
*以下の線形システム解法の後，配列XXにVS(V*)が入る
C   1. 直接法 : バンドマトリックスによるガウスの消去法
      IF (METHOD.EQ.1)
     $  CALL GB      (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                AGB)
C   2. 反復法1 : point-SOR 法
      IF (METHOD.EQ.2)
     $  CALL PSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY)
C   3. 反復法2 : line-SOR 法 : OMG=1で十分．大きくしすぎると発散する
      IF (METHOD.EQ.3)
     $  CALL LSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6,W7,W8,W9)
C   4. クリロフ部分空間法1 : 共役残差法
      IF (METHOD.EQ.4)
     $  CALL CRB     (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6)
C   5. クリロフ部分空間法2 : BiCGSTAB
      IF (METHOD.EQ.5)
     $  CALL BICGB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6,W7)
C   METHOD=4,5のときの探索ベクトル計算時のゼロ除算対策
      IF (IFLG.EQ.2)
     $  CALL PSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY)
*線形システム解法で得られた1次元配列の解を3次元配列に置き換える
      DO 100 IX = 1,NX
        DO 110 IY = 1,NY-1
          DO 120 IZ = 1,NZ
            K = IY + (IX-1)*NNY + (IZ-1)*(NNX*NNY)
            VS(IX,IY,IZ) = XX(K)
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE
*V*の境界条件の処理
      CALL VBND (VS)
*
*-------------------------------------------------*
*            W* : WS(IX,IY,IZ) の計算             *
*-------------------------------------------------*
* I:線形システム行列用の通し番号を意味する変数
* 3次元配列(IX,IY,IZ):離散化用 -> 1次元配列(I):線形システム用
      I = 1
      DO 130 IZ = 1,NZ-1
        DO 140 IX = 1,NX
          DO 150 IY = 1,NY
*         Wに関するDの計算
          DW5 = VIS * DY*DZ / DX
          DW1 = VIS * DY*DZ / DX
          DW4 = VIS * DX*DZ / DY
          DW2 = VIS * DX*DZ / DY
          DW7 = VIS * DX*DY / DZ
          DW6 = VIS * DX*DY / DZ
*         Wに関するFの計算
*         速度・圧力・温度場の繰り返し計算となっても線形システムの係数となる
*         ここの速度は更新されないように(U,V)には前の時刻の収束値(UO,VO)を用いる
          UW5 = ( UO(IX  ,IY  ,IZ  ) + UO(IX  ,IY  ,IZ+1) ) / 2.0D0
          UW1 = ( UO(IX-1,IY  ,IZ  ) + UO(IX-1,IY  ,IZ+1) ) / 2.0D0
          VW4 = ( VO(IX  ,IY  ,IZ  ) + VO(IX  ,IY  ,IZ+1) ) / 2.0D0
          VW2 = ( VO(IX  ,IY-1,IZ  ) + VO(IX  ,IY-1,IZ+1) ) / 2.0D0
          WW7 = ( WO(IX  ,IY  ,IZ  ) + WO(IX  ,IY  ,IZ+1) ) / 2.0D0
          WW6 = ( WO(IX  ,IY  ,IZ  ) + WO(IX  ,IY  ,IZ-1) ) / 2.0D0
          FW5 = RO * UW5 * DY * DZ
          FW1 = RO * UW1 * DY * DZ
          FW4 = RO * VW4 * DX * DZ
          FW2 = RO * VW2 * DX * DZ
          FW7 = RO * WW7 * DX * DY
          FW6 = RO * WW6 * DX * DY
*         Wに関するPeの計算
          PEW5 = FW5 / DW5
          PEW1 = FW1 / DW1
          PEW4 = FW4 / DW4
          PEW2 = FW2 / DW2
          PEW7 = FW7 / DW7
          PEW6 = FW6 / DW6
*         関数A(|Pe|)の計算
*         NMT=1なら 風上法
          IF ( NMT.EQ.1 ) THEN
            APEW5 = 1.0D0
            APEW1 = 1.0D0
            APEW4 = 1.0D0
            APEW2 = 1.0D0
            APEW7 = 1.0D0
            APEW6 = 1.0D0
*         NMT=2なら べき乗法
          ELSE
            APEW5 = DMAX1(0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEW5))**5 ) )
            APEW1 = DMAX1(0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEW1))**5 ) )
            APEW4 = DMAX1(0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEW4))**5 ) )
            APEW2 = DMAX1(0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEW2))**5 ) )
            APEW7 = DMAX1(0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEW7))**5 ) )
            APEW6 = DMAX1(0.0D0,DBLE( (1.0D0-0.1D0*ABS(PEW6))**5 ) )
          END IF
*         W* の線形システムの係数の計算
          AW5 = ( DW5*APEW5 + DMAX1(-FW5,0.0D0) )
          AW1 = ( DW1*APEW1 + DMAX1( FW1,0.0D0) )
          AW4 = ( DW4*APEW4 + DMAX1(-FW4,0.0D0) )
          AW2 = ( DW2*APEW2 + DMAX1( FW2,0.0D0) )
          AW7 = ( DW7*APEW7 + DMAX1(-FW7,0.0D0) )
          AW6 = ( DW6*APEW6 + DMAX1( FW6,0.0D0) )
*         生成項(浮力項)の計算
          SCW = 0.0D0
          SPW = 0.0D0
*         V* の線形システムのBの計算
          BP(I) = SCW*DX*DY*DZ+( PS(IX,IY,IZ)-PS(IX,IY,IZ+1) )*DX*DY
     $          + RO*DX*DY*DZ/DT*WO(IX,IY,IZ)-SPW*DX*DY*DZ
*境界の処理
* 1. 境界のW_BND についても変数とみなしてAW3(IX,IY,IZ)を求める
* 2. その後，境界の項 AW(1,2,4,5,6,7)*W_BND をBに吸収させる
* 3. 2に対応して(吸収の後)，境界のAW(1,2,4,5,6,7)をゼロとしてリンクを断つ
          AW3(IX,IY,IZ) = AW1+AW2+AW4+AW5+AW6+AW7 + RO*DX*DY*DZ/DT
*W（右面）
          IF ( IX.EQ.NX ) THEN
            BP(I) = BP(I) + AW5 * W(NX+1,IY,IZ)
            AW5 = 0.0D0
          END IF
*W（左面）
          IF ( IX.EQ.1 ) THEN
            BP(I) = BP(I) + AW1 * W(0,IY,IZ)
            AW1 = 0.0D0
          END IF
*W（上面）
          IF ( IY.EQ.NY ) THEN
            BP(I) = BP(I) + AW4*W(IX,NY+1,IZ)
            AW4 = 0.0D0
          END IF
*W（下面）
          IF ( IY.EQ.1 ) THEN
            BP(I) = BP(I) + AW2*W(IX,0,IZ)
            AW2 = 0.0D0
          END IF
*W（後面）
          IF ( IZ.EQ.NZ-1 ) THEN
            BP(I) = BP(I) + AW7*W(IX,IY,NZ)
            AW7 = 0.0D0
          END IF
*W（前面）
          IF ( IZ.EQ.1 ) THEN
            BP(I) = BP(I) + AW6*W(IX,IY,0)
            AW6 = 0.0D0
          END IF
* 線形システム計算用に係数行列を計算
* ALPHAW : W* 計算のための緩和係数
          ATP3(I) = AW3(IX,IY,IZ) / ALPHAW
          ATP5(I) = - AW5
          ATP1(I) = - AW1
          ATP4(I) = - AW4
          ATP2(I) = - AW2
          ATP7(I) = - AW7
          ATP6(I) = - AW6
          BP(I) = BP(I) + (1.0D0-ALPHAW)*ATP3(I)*WO(IX,IY,IZ)
* 線形システム解法前のための初期化(任意) (--- SIMPLE ---)
* METHOD=1,4,5のとき適宜修正
          XX(I) = W(IX,IY,IZ)
          I = I + 1
  150     CONTINUE
  140   CONTINUE
  130 CONTINUE
*
*サブルーチンに一般性を持たせるため，NX,NY,NZ,NEを引数として渡す
*COMMON文で定義されている値でもあるので，そのまま渡せない！
      NNX = NX
      NNY = NY
      NNZ = NZ-1
      NNE = NEW
      NNXY = NNX*NNY
*以下の線形システム解法の後，配列XXにVS(V*)が入る
C   1. 直接法 : バンドマトリックスによるガウスの消去法
      IF (METHOD.EQ.1)
     $  CALL GB      (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                AGB)
C   2. 反復法1 : point-SOR 法
      IF (METHOD.EQ.2)
     $  CALL PSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY)
C   3. 反復法2 : line-SOR 法 : OMG=1で十分．大きくしすぎると発散する
      IF (METHOD.EQ.3)
     $  CALL LSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6,W7,W8,W9)
C   4. クリロフ部分空間法1 : 共役残差法
      IF (METHOD.EQ.4)
     $  CALL CRB     (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6)
C   5. クリロフ部分空間法2 : BiCGSTAB
      IF (METHOD.EQ.5)
     $  CALL BICGB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6,W7)
C   METHOD=4,5のときの探索ベクトル計算時のゼロ除算対策
      IF (IFLG.EQ.2)
     $  CALL PSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY)
*線形システム解法で得られた1次元配列の解を3次元配列に置き換える
      DO 160 IX = 1,NX
        DO 170 IY = 1,NY
          DO 180 IZ = 1,NZ-1
            K = IY + (IX-1)*NNY + (IZ-1)*(NNX*NNY)
            WS(IX,IY,IZ) = XX(K)
  180     CONTINUE
  170   CONTINUE
  160 CONTINUE
*V*の境界条件の処理
      CALL WBND (WS)
*
      RETURN
      END
*********************************************************************
*                        圧力場の計算
*********************************************************************
*                                                                   *
*［圧力補正の線形システムに関する配列の説明］                       *
*                                                                   *
* 圧力補正(PD)の線形システムを ATP_{i,j} PD_{i} = BP_{i} とする．   *
* 有限差分近似を用いて離散化していることから明らかなように，        *
* 1. ATP_{i,j}の大部分はゼロ                                        *
* 2. 非ゼロ要素は疎であり，規則的に並んでいる                       *
* したがってATP_{i,j,K}すべてを記憶させるのは計算容量のムダであるので*
* 本プログラムでは，以下のような規則にしたがって記憶する．          *
* 3次元問題であれば，近隣の7点が関わる．                            *
*                                                                   *
* ●用いる配列(1次元配列に格納)：クリロフ部分空間法を含む反復法     *
*   係数行列     ---> ATP1, ATP2, ATP3, ATP4, ATP5, ATP6, ATP7      *
*   既知ベクトル ---> BP                                            *
*   未知ベクトル ---> XX                                            *
*   未知数の数(未知のPDの数)   ---> NEP = NX * NY * NZ              *
*                                                                   *
*    i=3, j=3, k=3 (m=59): PD(3,3,3) における ATP1,ATP2,ATP3,ATP4   *
*                                             ATP5,ATP6,ATP7        *
*        +------------------/-----+                                 *
* (NY)5  |    |    |    |  / |    |                                 *
*        +---------------ATP7-----+    通し番号の規則               *
*     4  |    |    |ATP4/    |    |    m=(k-1)*(NX*NY)+(i-1)*NY+j   *
*        +-------------/----------+                                 *
*     3  |    |ATP1|ATP3|ATP5|    |    ATP1(m)                      *
*        +---------/--------------+    ATP2(m)                      *
*     2  |    |   /|ATP2|    |    |    ATP3(m)                      *
*        +-----ATP6---------------+    ATP4(m)                      *
*     1  |      /  |    |    |    |    ATP5(m)                      *
*        +-----/------------------+    ATP6(m)                      *
*     j   i=1 /   2    3    4    5(NX) ATP7(m)                      *
*            /                         BP(m)                        *
*         k=1,5(NZ)                    NXY = NX * NY                *
*                                                                   *
*    BP(NE) : 既知ベクトル                                          *
*                                                                   *
*    XX(NE) : 未知ベクトル(各サブルーチンでこれが求まる)            *
*                                                                   *
* ●用いる配列：ガウスの消去法による直接法(バンド行列として格納)    *
*    AGB(-NY:NY,NEP) : 2次元差分近似による規則的非対称行列          *
*    -NY:NY->上の図で考えると K=13 のとき ATP1 から ATP5 までの k は*
*                   ATP1のkは K"-NY"                                *
*                   ATP2のkは K"-1 "                                *
*                   ATP3のkは K"+0 "                                *
*                   ATP4のkは K"+1 "                                *
*                   ATP5のkは K"+NY"                                *
*                   ATP6のkは K"-NXY"                               *
*                   ATP7のkは K"+NXY"                               *
*              このように，kを固定したとき，" "で囲まれた5個の値    *
*              (-NY,-1,0,1,NY,-NXY,NXY)                             *
*     NEP -> 上述のような k は全部で NEP 個定義される               *
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
      SUBROUTINE PRESS (AGB,W1,W2,W3,W4,W5,W6,W7,W8,W9)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NX0=20,NY0=20,NZ0=20,NXY0=400,NEU0=7600,
     $           NEV0=7600,NEW0=7600,NEP0=8000)
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP,EPSC,
     $              ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT
      COMMON / D5 / ICYCLE,ITR,IFLG,IFLGC,IRELP,METHOD,NMT
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
C     --- SIMPLE ---
      COMMON / D8 / NEU,NEV,NEW,NEP
      COMMON / D9 / NITR, NCNT
      COMMON / D10 / RO
      COMMON / ARRAY1 / AU3(0:NX0,0:NY0+1,0:NZ0+1),
     $AV3(0:NX0+1,0:NY0,0:NZ0+1),AW3(0:NX0+1,0:NY0+1,0:NZ0),
     $ATP1(NEP0),ATP2(NEP0),ATP3(NEP0),ATP4(NEP0),ATP5(NEP0),
     $ATP6(NEP0),ATP7(NEP0),BP(NEP0),
     $ U(0:NX0,  0:NY0+1,0:NZ0+1),  V(0:NX0+1,0:NY0,  0:NZ0+1),
     $ W(0:NX0+1,0:NY0+1,0:NZ0  ),
     $US(0:NX0,  0:NY0+1,0:NZ0+1), VS(0:NX0+1,0:NY0,  0:NZ0+1),
     $WS(0:NX0+1,0:NY0+1,0:NZ0  ),
     $PS(0:NX0+1,0:NY0+1,0:NZ0+1), PD(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ T(0:NX0+1,0:NY0+1,0:NZ0+1),
     $XX(NEP0)
      COMMON / ARRAY2 /
     $UO(0:NX0  ,0:NY0+1,0:NZ0+1),VO(0:NX0+1,0:NY0,  0:NZ0+1),
     $WO(0:NX0+1,0:NY0+1,0:NZ0  ),TO(0:NX0+1,0:NY0+1,0:NZ0+1)
*     作業用配列
      DIMENSION W1(NEP0),W2(NEP0),W3(NEP0),W4(NEP0),W5(NEP0),W6(NEP0),
     $          W7(NEP0),W8(NEP0),W9(NEP0)
*     バンドガウス消去法のための係数行列配列
      DIMENSION AGB(-NXY0:NXY0,NEP0)
*
*     速度・圧力・温度を連立させて繰り返し計算を行って，連続の式を満たす
*     ( U*, V* ,W* ) となったかどうかを判定する際の最大値の初期設定
      DMAX = 0.0D0
*     線形システム解法のため，2次元配列を1次元配列にするための変数
      I = 1
      DO 10 IZ = 1,NZ
        DO 20 IX = 1,NX
          DO 30 IY = 1,NY
*速度場(U*,V*,W*)が連続の式を満たすかどうかを計算し，満たしていれば，
*この(U*,V*,W*)をもとに，ここで計算するP'(PD)から得られる(U,V,W,P)を
*この時間ステップの収束解とする．なお，(U,V,W)は常に連続の式を満たす
*ことが保証される．したがって，(U*,V*,W*)が連続の式を満たすということは，
*前回の反復において得られた(U,V,W)を用いると，もはやこれを補正する必要が
*なくなるということを意味する．
*         P'(PD) の線形システムのBの計算
          BP(I)=RO*( ( US(IX-1,IY  ,IZ  )-US(IX,IY,IZ) )*DY*DZ
     $             + ( VS(IX  ,IY-1,IZ  )-VS(IX,IY,IZ) )*DX*DZ
     $             + ( WS(IX  ,IY  ,IZ-1)-WS(IX,IY,IZ) )*DX*DY )
*         (U*, V*)による連続の式の収束性をチェックし，div \VEC{v}の
*         最大値がDMAXとなるようにする
          IF ( DABS(BP(I)). GT. DMAX ) THEN
            DMAX = DABS(BP(I))
          END IF
*P'（右面）
          IF ( IX.EQ.NX ) THEN
            AP5 = 0.0D0
          ELSE
            APU5 = AU3(IX,IY,IZ)
            DP5 =  (DY*DZ) / APU5
            AP5 = RO * DP5 * DY * DZ
          END IF
*P'（左面）
          IF ( IX.EQ.1 ) THEN
            AP1 = 0.0D0
          ELSE
            APU1 = AU3(IX-1,IY,IZ)
            DP1 =  (DY*DZ) / APU1
            AP1 = RO * DP1 * DY * DZ
          END IF
*P'（上面）
          IF ( IY.EQ.NY ) THEN
            AP4 = 0.0D0
          ELSE
            APV4 = AV3(IX,IY,IZ)
            DP4 = (DX*DZ) / APV4
            AP4 = RO * DP4 * DX * DZ
          END IF
*P'（下面）
          IF ( IY.EQ.1 ) THEN
            AP2 = 0.0D0
          ELSE
            APV2 = AV3(IX,IY-1,IZ)
            DP2 = (DX*DZ) / APV2
            AP2 = RO * DP2 * DX * DZ
          END IF
*P'（後面）
          IF ( IZ.EQ.NZ ) THEN
            AP7 = 0.0D0
          ELSE
            APW7 = AW3(IX,IY,IZ)
            DP7 = (DX*DY) / APW7
            AP7 = RO * DP7 * DX * DY
          END IF
*P'（前面）
          IF ( IZ.EQ.1 ) THEN
            AP6 = 0.0D0
          ELSE
            APW6 = AW3(IX,IY,IZ-1)
            DP6 = (DX*DY) / APW6
            AP6 = RO * DP6 * DX * DY
          END IF
          AP3 = AP1 + AP2 + AP4 + AP5 + AP6 + AP7
* 線形システム計算用に係数行列を計算
          ATP3(I) =  AP3
          ATP5(I) = -AP5
          ATP1(I) = -AP1
          ATP4(I) = -AP4
          ATP2(I) = -AP2
          ATP7(I) = -AP7
          ATP6(I) = -AP6
*１次独立な解を得るための処理 : IRELP=1 : 直接法においては必須
*(IX=1,IY=1,IZ=1 ---> K=1を基準点とし，常にPN(1,1,1)=PD(1,1,1)=0とする)
          IF (IRELP.EQ.1.AND.IX.EQ.1.AND.IY.EQ.1.AND.IZ.EQ.1) THEN
            ATP4(1) = 0.0D0
            ATP5(1) = 0.0D0
            ATP7(1) = 0.0D0
            BP(1)   = 0.0D0
*           K=2 の点の処理(K=1とのリンクを断つ)
            ATP2(2) = 0.0D0
*           K=1+NY の点の処理(K=1とのリンクを断つ)
            ATP1(1+NY) = 0.0D0
*           K=1+NY+NX*NY の点の処理(K=1とのリンクを断つ)
            ATP6(1+NY+NX*NY) = 0.0D0
          END IF
* 線形システム解法前のための初期化(任意) (--- SIMPLE ---)
* METHOD=1,4,5のとき適宜修正
          XX(I) = PD(IX,IY,IZ)
          I =  I + 1
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
* --- SIMPLE --- 圧力補正 P'(PD) に関する解法
*サブルーチンに一般性を持たせるため，NX,NY,NEを引数として渡す
*COMMON文で定義されている値でもあるので，そのまま渡せない！
      NNX = NX
      NNY = NY
      NNZ = NZ
      NNE = NEP
      NNXY = NNX*NNY
C   1. 直接法 : バンドマトリックスによるガウスの消去法
      IF (METHOD.EQ.1)
     $  CALL GB      (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                AGB)
C   2. 反復法1 : point-SOR 法
      IF (METHOD.EQ.2)
     $  CALL PSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY)
C   3. 反復法2 : line-SOR 法 : OMG=1で十分．大きくしすぎると発散する
      IF (METHOD.EQ.3)
     $  CALL LSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6,W7,W8,W9)
C   4. クリロフ部分空間法1 : 共役残差法
      IF (METHOD.EQ.4)
     $  CALL CRB     (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6)
C   5. クリロフ部分空間法2 : BiCGSTAB
      IF (METHOD.EQ.5)
     $  CALL BICGB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6,W7)
C   METHOD=4,5のときの探索ベクトル計算時のゼロ除算対策
      IF (IFLG.EQ.2)
     $  CALL PSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY)
*
      DO 40 IX = 1,NX
        DO 50 IY = 1,NY
          DO 60 IZ = 1,NZ
          K = IY + (IX-1)*NNY + (IZ-1)*(NNX*NNY)
*         圧力の相対性の処理 : 基準値の設定
*         1次独立な解を求める場合は圧力の基準点が強制的に設定される
*         圧力の基準点を設けない場合 (IRELP=0,1)
          IF (IRELP.EQ.0. OR. IRELP.EQ.1) PD(IX,IY,IZ) = XX(K)
*         1次従属な解のうちの1つを求めた後，圧力基準を設ける場合 (IRELP=2)
*         P'(1,1,1)=0 ---> P(1,1,1)=0
          IF (IRELP.EQ.2) PD(IX,IY,IZ) = XX(K)-XX(1)
   60     CONTINUE
   50   CONTINUE
   40 CONTINUE
*圧力補正の境界条件の処理
      CALL PDBND (PD)
*速度の修正
      DO 70 IX = 1,NX-1
        DO 80 IY = 1,NY
          DO 90 IZ = 1,NZ
            DU = DY*DZ / AU3(IX,IY,IZ)
            U(IX,IY,IZ)=US(IX,IY,IZ)+DU*(PD(IX,IY,IZ)-PD(IX+1,IY,IZ))
   90     CONTINUE
   80   CONTINUE
   70 CONTINUE
      DO 100 IX = 1,NX
        DO 110 IY = 1,NY-1
          DO 120 IZ = 1,NZ
            DV = DX*DZ / AV3(IX,IY,IZ)
            V(IX,IY,IZ)=VS(IX,IY,IZ)+DV*(PD(IX,IY,IZ)-PD(IX,IY+1,IZ))
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE
      DO 130 IX = 1,NX
        DO 140 IY = 1,NY
          DO 150 IZ = 1,NZ-1
            DW = DX*DY / AW3(IX,IY,IZ)
            W(IX,IY,IZ)=WS(IX,IY,IZ)+DW*(PD(IX,IY,IZ)-PD(IX,IY,IZ+1))
  150     CONTINUE
  140   CONTINUE
  130 CONTINUE
*新たに得られた速度を用いて境界条件を処理する
      CALL UBND (U)
      CALL VBND (V)
      CALL WBND (W)
*P* の修正: P = P* + P': ただしPの値はPSに代入
      DO 160 IX = 1,NX
        DO 170 IY = 1,NY
          DO 180 IZ = 1,NZ
            PS(IX,IY,IZ) = PS(IX,IY,IZ) + ALPHAP*PD(IX,IY,IZ)
  180     CONTINUE
  170   CONTINUE
  160 CONTINUE
*新たに得られたP*を用いて境界条件を処理する
      CALL PSBND (PS)
*速度場 (U*,V*) が EPSC 以下で連続の式を満たすかどうかを判定
      IF ( DMAX. GT. EPSC ) IFLGC=1
      RETURN
      END
*********************************************************************
*                        温度場の計算
*********************************************************************
        SUBROUTINE CALTEM (AGB,W1,W2,W3,W4,W5,W6,W7,W8,W9)
*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NX0=20,NY0=20,NZ0=20,NXY0=400,NEU0=7600,
     $           NEV0=7600,NEW0=7600,NEP0=8000)
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP,EPSC,
     $              ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT
      COMMON / D5 / ICYCLE,ITR,IFLG,IFLGC,IRELP,METHOD,NMT
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
C     --- SIMPLE ---
      COMMON / D8 / NEU,NEV,NEW,NEP
      COMMON / D9 / NITR, NCNT
      COMMON / D10 / RO
      COMMON / ARRAY1 / AU3(0:NX0,0:NY0+1,0:NZ0+1),
     $AV3(0:NX0+1,0:NY0,0:NZ0+1),AW3(0:NX0+1,0:NY0+1,0:NZ0),
     $ATP1(NEP0),ATP2(NEP0),ATP3(NEP0),ATP4(NEP0),ATP5(NEP0),
     $ATP6(NEP0),ATP7(NEP0),BP(NEP0),
     $ U(0:NX0,  0:NY0+1,0:NZ0+1),  V(0:NX0+1,0:NY0,  0:NZ0+1),
     $ W(0:NX0+1,0:NY0+1,0:NZ0  ),
     $US(0:NX0,  0:NY0+1,0:NZ0+1), VS(0:NX0+1,0:NY0,  0:NZ0+1),
     $WS(0:NX0+1,0:NY0+1,0:NZ0  ),
     $PS(0:NX0+1,0:NY0+1,0:NZ0+1), PD(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ T(0:NX0+1,0:NY0+1,0:NZ0+1),
     $XX(NEP0)
      COMMON / ARRAY2 /
     $UO(0:NX0  ,0:NY0+1,0:NZ0+1),VO(0:NX0+1,0:NY0,  0:NZ0+1),
     $WO(0:NX0+1,0:NY0+1,0:NZ0  ),TO(0:NX0+1,0:NY0+1,0:NZ0+1)
*     作業用配列
      DIMENSION W1(NEP0),W2(NEP0),W3(NEP0),W4(NEP0),W5(NEP0),W6(NEP0),
     $          W7(NEP0),W8(NEP0),W9(NEP0)
*     バンドガウス消去法のための係数行列配列
      DIMENSION AGB(-NXY0:NXY0,NEP0)
*
*     線形システム解法のため，3次元配列を1次元配列にするための変数
      I = 1
      DO 10 IZ = 1,NZ
        DO 20 IX = 1,NX
          DO 30 IY = 1,NY
*         Tに関するDの計算
          DT5 = ALP * DY*DZ / DX
          DT1 = ALP * DY*DZ / DX
          DT4 = ALP * DX*DZ / DY
          DT2 = ALP * DX*DZ / DY
          DT7 = ALP * DX*DY / DZ
          DT6 = ALP * DX*DY / DZ
*         Tに関するFの計算
*         速度・圧力・温度場の繰り返し計算の際のUとVを線形システムの係数に反映させる
          UT5 = U(IX  ,IY  ,IZ  )
          UT1 = U(IX-1,IY  ,IZ  )
          VT4 = V(IX  ,IY  ,IZ  )
          VT2 = V(IX  ,IY-1,IZ  )
          WT7 = W(IX  ,IY,  IZ  )
          WT6 = W(IX  ,IY,  IZ-1)
          FT5 = RO * UT5 * DY * DZ
          FT1 = RO * UT1 * DY * DZ
          FT4 = RO * VT4 * DX * DZ
          FT2 = RO * VT2 * DX * DZ
          FT7 = RO * WT7 * DX * DY
          FT6 = RO * WT6 * DX * DY
*         Tに関するPeの計算
          PET5 = FT5 / DT5
          PET1 = FT1 / DT1
          PET4 = FT4 / DT4
          PET2 = FT2 / DT2
          PET7 = FT7 / DT7
          PET6 = FT6 / DT6
*         関数A(|Pe|)の計算
*         NMT=1なら 風上法
          IF ( NMT.EQ.1 ) THEN
            APET5 = 1.0D0
            APET1 = 1.0D0
            APET4 = 1.0D0
            APET2 = 1.0D0
            APET7 = 1.0D0
            APET6 = 1.0D0
*         NMT=2なら べき乗法
          ELSE
            APET5 = DMAX1( 0.0D0,DBLE( (1.0D0-0.1D0*ABS(PET5))**5 ) )
            APET1 = DMAX1( 0.0D0,DBLE( (1.0D0-0.1D0*ABS(PET1))**5 ) )
            APET4 = DMAX1( 0.0D0,DBLE( (1.0D0-0.1D0*ABS(PET4))**5 ) )
            APET2 = DMAX1( 0.0D0,DBLE( (1.0D0-0.1D0*ABS(PET2))**5 ) )
            APET7 = DMAX1( 0.0D0,DBLE( (1.0D0-0.1D0*ABS(PET7))**5 ) )
            APET6 = DMAX1( 0.0D0,DBLE( (1.0D0-0.1D0*ABS(PET6))**5 ) )
          END IF
*         T の線形システムの係数の計算
          AT5 = ( DT5*APET5 + DMAX1(-FT5,0.0D0) )
          AT1 = ( DT1*APET1 + DMAX1( FT1,0.0D0) )
          AT4 = ( DT4*APET4 + DMAX1(-FT4,0.0D0) )
          AT2 = ( DT2*APET2 + DMAX1( FT2,0.0D0) )
          AT7 = ( DT2*APET7 + DMAX1(-FT7,0.0D0) )
          AT6 = ( DT2*APET6 + DMAX1( FT6,0.0D0) )
*         生成項
          SCT = 0.0D0
          SPT = 0.0D0
*         T の線形システムのBの計算
          BP(I) = SCT*DX*DY*DZ + RO*DX*DY*DZ/DT*TO(IX,IY,IZ)
*境界の処理
* 1. 境界のT_BND についても変数とみなしてAT3を求める
* 2. その後，境界の項 AT3*T_BND をBに吸収させる
* 3. 2に対応して(吸収の後)，境界のAT(1,2,4,5,6,7)をゼロとしてリンクを断つ
          AT3 = AT1+AT2+AT4+AT5+AT6+AT7 + RO*DX*DY*DZ/DT - SPT
*T（右面）
          IF (IX.EQ.NX) THEN
            BP(I) = BP(I) + AT5 * T(NX+1,IY,IZ)
            AT5 = 0.0D0
*T（左面）
          ELSE IF (IX.EQ.1) THEN
            BP(I) = BP(I) + AT1 * T(0,IY,IZ)
            AT1 = 0.0D0
          END IF
*T（上面）
          IF (IY.EQ.NY) THEN
            BP(I) = BP(I) + AT4 * T(IX,NY+1,IZ)
            AT4 = 0.0D0
*T（下面）
          ELSE IF (IY.EQ.1) THEN
            BP(I) = BP(I) + AT2 * T(IX,0,IZ)
            AT2 = 0.0D0
          END IF
*T（後面）
          IF (IZ.EQ.NZ) THEN
            BP(I) = BP(I) + AT7 * T(IX,IY,NZ+1)
            AT7 = 0.0D0
*T（前面）
          ELSE IF (IZ.EQ.1) THEN
            BP(I) = BP(I) + AT6 * T(IX,IY,0)
            AT6 = 0.0D0
          END IF
* 線形システム計算用に係数行列を計算
* ALPHAT : T 計算のための緩和係数
          ATP3(I) =   AT3 / ALPHAT
          ATP5(I) = - AT5
          ATP1(I) = - AT1
          ATP4(I) = - AT4
          ATP2(I) = - AT2
          ATP7(I) = - AT7
          ATP6(I) = - AT6
          BP(I) = BP(I) + (1.0D0-ALPHAT)*ATP3(I)*TO(IX,IY,IZ)
* 線形システム解法前のための初期化(任意) (--- SIMPLE ---)
* METHOD=1,4,5のとき適宜修正
          XX(I) = T(IX,IY,IZ)
          I = I + 1
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
*
*サブルーチンに一般性を持たせるため，NX,NY,NEを引数として渡す
*COMMON文で定義されている値でもあるので，そのまま渡せない！
      NNX = NX
      NNY = NY
      NNZ = NZ
      NNE = NEP
      NNXY = NNX*NNY
C   1. 直接法 : バンドマトリックスによるガウスの消去法
      IF (METHOD.EQ.1)
     $  CALL GB      (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                AGB)
C   2. 反復法1 : point-SOR 法
      IF (METHOD.EQ.2)
     $  CALL PSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY)
C   3. 反復法2 : line-SOR 法 : OMG=1で十分．大きくしすぎると発散する
      IF (METHOD.EQ.3)
     $  CALL LSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6,W7,W8,W9)
C   4. クリロフ部分空間法1 : 共役残差法
      IF (METHOD.EQ.4)
     $  CALL CRB     (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6)
C   5. クリロフ部分空間法2 : BiCGSTAB
      IF (METHOD.EQ.5)
     $  CALL BICGB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY,
     $                W1,W2,W3,W4,W5,W6,W7)
C   METHOD=4,5のときの探索ベクトル計算時のゼロ除算対策
      IF (IFLG.EQ.2)
     $  CALL PSORB   (ATP1,ATP2,ATP3,ATP4,ATP5,ATP6,ATP7,BP,XX,
     $                NNX,NNY,NNZ,NNE,NNXY)
*線形システム解法で得られた1次元配列の解を2次元配列に置き換える
      DO 40 IX = 1,NX
        DO 50 IY = 1,NY
          DO 60 IZ = 1,NZ
            K = IY + (IX-1)*NNY + (IZ-1)*(NNX*NNY)
            T(IX,IY,IZ) = XX(K)
   60     CONTINUE
   50   CONTINUE
   40 CONTINUE
*境界条件の処理
      CALL TBND (T)
      RETURN
      END
*********************************************************************
*                速度 ( U* と U ) の境界条件の処理
*********************************************************************
      SUBROUTINE UBND (XU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D1 / NX,NY,NZ
      DIMENSION XU(0:NX,0:NY+1,0:NZ+1)
*U（左面）
      DO 10 IY = 1,NY
        DO 20 IZ = 1,NZ
          XU(0,IY,IZ) = 0.0D0
   20   CONTINUE
   10 CONTINUE
*U（右面）
      DO 30 IY = 1,NY
        DO 40 IZ = 1,NZ
          XU(NX,IY,IZ) = 0.0D0
   40   CONTINUE
   30 CONTINUE
*U（下面）
      DO 50 IX = 0,NX
        DO 60 IZ = 1,NZ
          XU(IX,0,IZ) = -XU(IX,1,IZ)
   60   CONTINUE
   50 CONTINUE
*U（上面）
      DO 70 IX = 0,NX
        DO 80 IZ = 1,NZ
          XU(IX,NY+1,IZ) = -XU(IX,NY,IZ)
   80   CONTINUE
   70 CONTINUE
*U（前面）
      DO 90 IX = 0,NX
        DO 100 IY = 0,NY+1
          XU(IX,IY,0) = -XU(IX,IY,1)
  100   CONTINUE
   90 CONTINUE
*U（後面）
      DO 110 IX = 0,NX
        DO 120 IY = 0,NY+1
          XU(IX,IY,NZ+1) = -XU(IX,IY,NZ)
  120   CONTINUE
  110 CONTINUE
*
      RETURN
      END
*********************************************************************
*                速度 ( V* と V ) の境界条件の処理
*********************************************************************
      SUBROUTINE VBND (XV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D1 / NX,NY,NZ
      DIMENSION XV(0:NX+1,0:NY,0:NZ+1)
*V（左面）
      DO 10 IY = 1,NY-1
        DO 20 IZ = 1,NZ
          XV(0,IY,IZ) = -XV(1,IY,IZ)
   20   CONTINUE
   10 CONTINUE
*V（右面）
      DO 30 IY = 1,NY-1
        DO 40 IZ = 1,NZ
          XV(NX+1,IY,IZ) = -XV(NX,IY,IZ)
   40   CONTINUE
   30 CONTINUE
*V（下面）
      DO 50 IX = 0,NX+1
        DO 60 IZ = 1,NZ
          XV(IX,0,IZ) = 0.0D0
   60 CONTINUE
   50 CONTINUE
*V（上面）
      DO 70 IX = 0,NX+1
        DO 80 IZ = 1,NZ
          XV(IX,NY,IZ) = 0.0D0
   80   CONTINUE
   70 CONTINUE
*V（前面）
      DO  90 IX = 0,NX+1
        DO 100 IY = 0,NY
          XV(IX,IY,0) = -XV(IX,IY,1)
  100   CONTINUE
   90 CONTINUE
*V（後面）
      DO 110 IX = 0,NX+1
        DO 120 IY = 0,NY
          XV(IX,IY,NZ+1) = -XV(IX,IY,NZ)
  120   CONTINUE
  110 CONTINUE
*
      RETURN
      END
*********************************************************************
*                速度 ( W* と W ) の境界条件の処理
*********************************************************************
      SUBROUTINE WBND (XW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D1 / NX,NY,NZ
      DIMENSION XW(0:NX+1,0:NY+1,0:NZ)
*W（左面）
      DO 10 IY = 1,NY
        DO 20 IZ = 1,NZ-1
          XW(0,IY,IZ) = -XW(1,IY,IZ)
   20   CONTINUE
   10 CONTINUE
*W（右面）
      DO 30 IY = 1,NY
        DO 40 IZ = 1,NZ-1
          XW(NX+1,IY,IZ) = -XW(NX,IY,IZ)
   40   CONTINUE
   30 CONTINUE
*W（下面）
      DO 50 IX = 0,NX+1
        DO 60 IZ = 1,NZ-1
          XW(IX,0,IZ) = -XW(IX,1,IZ)
   60   CONTINUE
   50 CONTINUE
*W（上面）
      DO 70 IX = 0,NX+1
        DO 80 IZ = 1,NZ-1
          XW(IX,NY+1,IZ) = -XW(IX,NY,IZ)
   80   CONTINUE
   70 CONTINUE
*W（前面）
      DO 90 IX = 0,NX+1
        DO 100 IY = 0,NY+1
          XW(IX,IY,0) = 0.0D0
  100   CONTINUE
   90 CONTINUE
*W（後面）
      DO 110 IX = 0,NX+1
        DO 120 IY = 0,NY+1
          XW(IX,IY,NZ) = 0.0D0
  120   CONTINUE
  110 CONTINUE
      RETURN
      END
*********************************************************************
*                  温度の境界条件の処理
*********************************************************************
      SUBROUTINE TBND (XT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D1 / NX,NY,NZ
      DIMENSION XT(0:NX+1,0:NY+1,0:NZ+1)
* 左面
      DO 10 IY = 1,NY
        DO 20 IZ = 1,NZ
          XT(0,IY,IZ) = 2.0D0 * ( +0.5D0 ) - XT(1,IY,IZ)
   20   CONTINUE
   10 CONTINUE
* 右面
      DO 30 IY = 1,NY
        DO 40 IZ = 1,NZ
          XT(NX+1,IY,IZ) = 2.0D0 * ( -0.5D0 ) - XT(NX,IY,IZ)
   40   CONTINUE
   30 CONTINUE
* 下面
      DO 50 IX = 0,NX+1
        DO 60 IZ = 1,NZ
          XT(IX,0,IZ) = XT(IX,1,IZ)
   60   CONTINUE
   50 CONTINUE
* 上面
      DO 70 IX = 0,NX+1
        DO 80 IZ = 1,NZ
          XT(IX,NY+1,IZ) = XT(IX,NY,IZ)
   80   CONTINUE
   70 CONTINUE
*前面
      DO 90 IX = 0,NX+1
        DO 100 IY = 0,NY+1
           XT(IX,IY,0) = XT(IX,IY,1)
  100   CONTINUE
   90 CONTINUE
*後面
      DO 110 IX = 0,NX+1
        DO 120 IY = 0,NY+1
          XT(IX,IY,NZ+1) = XT(IX,IY,NZ)
  120   CONTINUE
  110 CONTINUE
*
      RETURN
      END
*********************************************************************
*                      圧力補正の境界条件                           *
*********************************************************************
      SUBROUTINE PDBND (XP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D1 / NX,NY,NZ
      DIMENSION XP(0:NX+1,0:NY+1,0:NZ+1)
* 左面
      DO 10 IY = 1,NY
        DO 20 IZ = 1,NZ
          XP(0,IY,IZ) = XP(1,IY,IZ)
   20   CONTINUE
   10 CONTINUE
* 右面
      DO 30 IY = 1,NY
        DO 40 IZ = 1,NZ
          XP(NX+1,IY,IZ) = XP(NX,IY,IZ)
   40   CONTINUE
   30 CONTINUE
* 下面
      DO 50 IX = 0,NX+1
        DO 60 IZ = 1,NZ
          XP(IX,0,IZ) = XP(IX,1,IZ)
   60   CONTINUE
   50 CONTINUE
* 上面
      DO 70 IX = 0,NX+1
        DO 80 IZ = 1,NZ
          XP(IX,NY+1,IZ) = XP(IX,NY,IZ)
   80   CONTINUE
   70 CONTINUE
*後面
      DO 90 IX = 0,NX+1
        DO 100 IY = 0,NY+1
          XP(IX,IY,NZ+1) = XP(IX,IY,NZ)
  100   CONTINUE
   90 CONTINUE
*前面
      DO 110 IX = 0,NX+1
        DO 120 IY = 0,NY+1
          XP(IX,IY,0) = XP(IX,IY,1)
  120   CONTINUE
  110 CONTINUE
      RETURN
      END
*********************************************************************
*                      圧力の境界条件                               *
*********************************************************************
      SUBROUTINE PSBND (XP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D1 / NX,NY,NZ
      DIMENSION XP(0:NX+1,0:NY+1,0:NZ+1)
* 左面
      DO 10 IY = 1,NY
        DO 20 IZ = 1,NZ
          XP(0,IY,IZ) = XP(1,IY,IZ)
   20   CONTINUE
   10 CONTINUE
* 右面
      DO 30 IY = 1,NY
        DO 40 IZ = 1,NZ
          XP(NX+1,IY,IZ) = XP(NX,IY,IZ)
   40   CONTINUE
   30 CONTINUE
* 下面
      DO 50 IX = 0,NX+1
        DO 60 IZ = 1,NZ
          XP(IX,0,IZ) = XP(IX,1,IZ)
   60   CONTINUE
   50 CONTINUE
* 上面
      DO 70 IX = 0,NX+1
        DO 80 IZ = 1,NZ
          XP(IX,NY+1,IZ) = XP(IX,NY,IZ)
   80   CONTINUE
   70 CONTINUE
*後面
      DO 90 IX = 0,NX+1
        DO 100 IY = 0,NY+1
          XP(IX,IY,NZ+1) = XP(IX,IY,NZ)
  100   CONTINUE
   90 CONTINUE
*前面
      DO 110 IX = 0,NX+1
        DO 120 IY = 0,NY+1
          XP(IX,IY,0) = XP(IX,IY,1)
  120   CONTINUE
  110 CONTINUE
*
      RETURN
      END
*********************************************************************
*                       データ出力
*********************************************************************
      SUBROUTINE PROUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NX0=20,NY0=20,NZ0=20,NXY0=400,NEU0=7600,
     $           NEV0=7600,NEW0=7600,NEP0=8000)
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP,EPSC,
     $              ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT
      COMMON / D5 / ICYCLE,ITR,IFLG,IFLGC,IRELP,METHOD,NMT
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
C     --- SIMPLE ---
      COMMON / D8 / NEU,NEV,NEW,NEP
      COMMON / D9 / NITR, NCNT
      COMMON / D10 / RO
      COMMON / ARRAY1 / AU3(0:NX0,0:NY0+1,0:NZ0+1),
     $AV3(0:NX0+1,0:NY0,0:NZ0+1),AW3(0:NX0+1,0:NY0+1,0:NZ0),
     $ATP1(NEP0),ATP2(NEP0),ATP3(NEP0),ATP4(NEP0),ATP5(NEP0),
     $ATP6(NEP0),ATP7(NEP0),BP(NEP0),
     $ U(0:NX0,  0:NY0+1,0:NZ0+1),  V(0:NX0+1,0:NY0,  0:NZ0+1),
     $ W(0:NX0+1,0:NY0+1,0:NZ0  ),
     $US(0:NX0,  0:NY0+1,0:NZ0+1), VS(0:NX0+1,0:NY0,  0:NZ0+1),
     $WS(0:NX0+1,0:NY0+1,0:NZ0  ),
     $PS(0:NX0+1,0:NY0+1,0:NZ0+1), PD(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ T(0:NX0+1,0:NY0+1,0:NZ0+1),
     $XX(NEP0)
      COMMON / ARRAY2 /
     $UO(0:NX0  ,0:NY0+1,0:NZ0+1),VO(0:NX0+1,0:NY0,  0:NZ0+1),
     $WO(0:NX0+1,0:NY0+1,0:NZ0  ),TO(0:NX0+1,0:NY0+1,0:NZ0+1)
*
      WRITE (11) U
      WRITE (12) V
      WRITE (13) W
      WRITE (14) PS,PD
      WRITE (15) T
      RETURN
      END
*********************************************************************
*                    Tecplot用データ出力
*********************************************************************
      SUBROUTINE TECPLT (FNAME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NX0=20,NY0=20,NZ0=20,NXY0=400,NEU0=7600,
     $           NEV0=7600,NEW0=7600,NEP0=8000)
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP,EPSC,
     $              ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT
      COMMON / D5 / ICYCLE,ITR,IFLG,IFLGC,IRELP,METHOD,NMT
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
C     --- SIMPLE ---
      COMMON / D8 / NEU,NEV,NEW,NEP
      COMMON / D9 / NITR, NCNT
      COMMON / D10 / RO
      COMMON / ARRAY1 / AU3(0:NX0,0:NY0+1,0:NZ0+1),
     $AV3(0:NX0+1,0:NY0,0:NZ0+1),AW3(0:NX0+1,0:NY0+1,0:NZ0),
     $ATP1(NEP0),ATP2(NEP0),ATP3(NEP0),ATP4(NEP0),ATP5(NEP0),
     $ATP6(NEP0),ATP7(NEP0),BP(NEP0),
     $ U(0:NX0,  0:NY0+1,0:NZ0+1),  V(0:NX0+1,0:NY0,  0:NZ0+1),
     $ W(0:NX0+1,0:NY0+1,0:NZ0  ),
     $US(0:NX0,  0:NY0+1,0:NZ0+1), VS(0:NX0+1,0:NY0,  0:NZ0+1),
     $WS(0:NX0+1,0:NY0+1,0:NZ0  ),
     $PS(0:NX0+1,0:NY0+1,0:NZ0+1), PD(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ T(0:NX0+1,0:NY0+1,0:NZ0+1),
     $XX(NEP0)
      COMMON / ARRAY2 /
     $UO(0:NX0  ,0:NY0+1,0:NZ0+1),VO(0:NX0+1,0:NY0,  0:NZ0+1),
     $WO(0:NX0+1,0:NY0+1,0:NZ0  ),TO(0:NX0+1,0:NY0+1,0:NZ0+1)
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
             UU = ( U(IX  ,IY,  IZ  )+U(IX  ,IY+1,IZ  )
     $             +U(IX  ,IY,  IZ+1)+U(IX  ,IY+1,IZ+1) )/4.0D0
             VV = ( V(IX  ,IY,  IZ  )+V(IX+1,IY  ,IZ  )
     $             +V(IX  ,IY,  IZ+1)+V(IX+1,IY,  IZ+1) )/4.0D0
             WW = ( W(IX  ,IY  ,IZ  )+W(IX  ,IY+1,IZ  )
     $             +W(IX+1,IY  ,IZ  )+W(IX+1,IY+1,IZ  ) )/4.0D0
             TT = ( T(IX  ,IY  ,IZ  )+T(IX+1,IY  ,IZ  )
     $             +T(IX  ,IY+1,IZ  )+T(IX+1,IY+1,IZ  )
     $             +T(IX  ,IY  ,IZ+1)+T(IX+1,IY  ,IZ+1)
     $             +T(IX  ,IY+1,IZ+1)+T(IX+1,IY+1,IZ+1) )/8.0D0
            WRITE (21,4010) X,Y,Z,UU,VV,WW,TT
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
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP,EPSC,
     $              ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT
      COMMON / D5 / ICYCLE,ITR,IFLG,IFLGC,IRELP,METHOD,NMT
C     --- SIMPLE ---
      COMMON / D9 / NITR, NCNT
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
*.....係数行列が特異なら計算終了する (--- SIMPLE ---)
      IF ( DABS(A(0,NE)).LE.1.0D-50 ) THEN
        WRITE (6,*) ' Matrix is singular : |A(0,NE)| < 1E-50 '
        IFLG = 1
        RETURN
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
*.....線形システムの総反復回数出力のために追加 (--- SIMPLE ---)
      ITR = 1
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
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP,EPSC,
     $              ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT
      COMMON / D5 / ICYCLE,ITR,IFLG,IFLGC,IRELP,METHOD,NMT
      COMMON / D9 / NITR, NCNT
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
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP,EPSC,
     $              ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT
      COMMON / D5 / ICYCLE,ITR,IFLG,IFLGC,IRELP,METHOD,NMT
      COMMON / D9 / NITR, NCNT
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
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP,EPSC,
     $              ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT
      COMMON / D5 / ICYCLE,ITR,IFLG,IFLGC,IRELP,METHOD,NMT
      COMMON / D9 / NITR, NCNT
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
*.......探索方向のための計算が 0除算 なら計算終了(--- SIMPLE ---)
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
*.......探索方向のための計算が 0除算 なら計算終了(--- SIMPLE ---)
        IF (DABS(APAP).LT.1.0D-50) THEN
          WRITE (6,*) ' 0 division : BETA_{K} in Conjugate Residual '
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
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP,EPSC,
     $              ALPHAP,ALPHAU,ALPHAV,ALPHAW,ALPHAT
      COMMON / D5 / ICYCLE,ITR,IFLG,IFLGC,IRELP,METHOD,NMT
      COMMON / D9 / NITR, NCNT
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
*.......探索方向のための計算が 0除算 なら計算終了(--- SIMPLE ---)
        IF (DABS(SAP).LT.1.0D-50) THEN
          WRITE (6,*) ' 0 division : ALPHA_{K} in Bi-CGSTAB '
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
*.......探索方向のための計算が 0除算 なら計算終了(--- SIMPLE ---)
        IF (DABS(ATAT).LT.1.0D-50) THEN
          WRITE (6,*) ' 0 division : XI_{K} in Bi-CGSTAB '
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
