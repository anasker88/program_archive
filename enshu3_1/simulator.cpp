// simulatotion of 2D fluid dynamics
#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

template <typename T>
using vec2d = vector<vector<T>>;

FILE *output;
vector<char> buf(1 << 20);

// 問題設定
const int NX = 100; // x方向の格子点数
const int NY = 100; // y方向の格子点数
int MAX_CYCLE;
double EPSP, OMG, DT, L, G, BETA, THETA_H, THETA_C, A, NU;
double DX = 1.0 / NX; // x方向の格子間隔
double DY = 1.0 / NY; // y方向の格子間隔
double PR, GR, RA;    // Prandtl数, Grashof数, Rayleigh数
int mymax = 0;

vec2d<double> U(NX + 1, vector<double>(NY + 2));   // x方向速度
vec2d<double> V(NX + 2, vector<double>(NY + 1));   // y方向速度
vec2d<double> U_T(NX + 1, vector<double>(NY + 2)); // x方向速度の仮値
vec2d<double> V_T(NX + 2, vector<double>(NY + 1)); // y方向速度の仮値
vec2d<double> P(NX + 2, vector<double>(NY + 2));   // 圧力
vec2d<double> T(NX + 2, vector<double>(NY + 2));   // 温度
// vec2d<double> U_D(NX + 1, vector<double>(NY + 2)); // x方向速度の修正値
// vec2d<double> V_D(NX + 2, vector<double>(NY + 1)); // y方向速度の修正値
vec2d<double> P_D(NX + 2, vector<double>(NY + 2));     // 圧力の修正値
vec2d<double> P_D_new(NX + 2, vector<double>(NY + 2)); // 圧力の修正値(更新用)
// 各項
vec2d<double> CNVU(NX + 1, vector<double>(NY + 2));
vec2d<double> CNVV(NX + 2, vector<double>(NY + 1));
vec2d<double> CNVT(NX + 2, vector<double>(NY + 2));
vec2d<double> DIFU(NX + 1, vector<double>(NY + 2));
vec2d<double> DIFV(NX + 2, vector<double>(NY + 1));
vec2d<double> DIFT(NX + 2, vector<double>(NY + 2));
vec2d<double> BUOV(NX + 2, vector<double>(NY + 1));
vec2d<double> DIV(NX + 2, vector<double>(NY + 2));

void read_input();
void init();
void calc_terms();
void calc_UV_T();
void calc_DP();
void calc_UVPT();
void set_boundary();
void write_result(int cycle);

int main()
{
    read_input();
    init();
    for (int cycle = 0; cycle < MAX_CYCLE; cycle++)
    {
        calc_terms();
        calc_UV_T();
        calc_DP();
        calc_UVPT();
        set_boundary();
        if ((cycle & 0x1f) == 0)
        {
            write_result(cycle);
        }
    }
    // ファイルを閉じる
    fclose(output);
    return 0;
}

// パラメータ読み込み
void read_input()
{
    // parameter.datからパラメータを読み込む
    ifstream ifs("parameters.dat");
    // 1行ずつ読み込む(2つ目がパラメータの値)
    string str;
    getline(ifs, str);
    MAX_CYCLE = stoi(str.substr(str.find(" ") + 1));
    getline(ifs, str);
    EPSP = stod(str.substr(str.find(" ") + 1));
    getline(ifs, str);
    OMG = stod(str.substr(str.find(" ") + 1));
    getline(ifs, str);
    DT = stod(str.substr(str.find(" ") + 1));
    getline(ifs, str);
    L = stod(str.substr(str.find(" ") + 1));
    getline(ifs, str);
    G = stod(str.substr(str.find(" ") + 1));
    getline(ifs, str);
    BETA = stod(str.substr(str.find(" ") + 1));
    getline(ifs, str);
    THETA_H = stod(str.substr(str.find(" ") + 1));
    getline(ifs, str);
    THETA_C = stod(str.substr(str.find(" ") + 1));
    getline(ifs, str);
    A = stod(str.substr(str.find(" ") + 1));
    getline(ifs, str);
    NU = stod(str.substr(str.find(" ") + 1));
    ifs.close();
    PR = NU / A;
    GR = G * BETA * (THETA_H - THETA_C) * L * L * L / (NU * NU);
    RA = GR * PR;
    cout << "PR=" << PR << endl;
    cout << "RA=" << RA << endl;
    cout << "GR=" << GR << endl;
    cout << "DX=" << DX << endl;
    cout << "DY=" << DY << endl;
}

void init()
{
    // 出力ファイルの設定
    output = fopen("result", "wb");
    // NX, NY
    fwrite(&NX, sizeof(int), 1, output);
    fwrite(&NY, sizeof(int), 1, output);
    // 初期条件
    // 速度は0
    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 2; j++)
        {
            U.at(i).at(j) = 0.0;
        }
    }
    for (int i = 0; i < NX + 2; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            V.at(i).at(j) = 0.0;
        }
    }
    // 圧力は0
    for (int i = 0; i < NX + 2; i++)
    {
        for (int j = 0; j < NY + 2; j++)
        {
            P.at(i).at(j) = 0.0;
        }
    }
    // 境界以外で温度は0
    // 上下の境界も0
    for (int i = 1; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 2; j++)
        {
            T.at(i).at(j) = 0.0;
        }
    }
    for (int i = 0; i < NY + 2; i++)
    {
        T.at(0).at(i) = 0.5;       // 左側:加熱
        T.at(NX + 1).at(i) = -0.5; // 右側:冷却
    }
}

// x方向対流項
double calc_CNVU(int i, int j)
{
    double V_av = 0.25 * (V.at(i).at(j) + V.at(i).at(j - 1) + V.at(i + 1).at(j) + V.at(i + 1).at(j - 1));
    // return (U.at(i).at(j) * (U.at(i + 1).at(j) - U.at(i - 1).at(j)) - abs(U.at(i).at(j)) * (U.at(i + 1).at(j) + U.at(i - 1).at(j) - 2 * U.at(i).at(j))) / (2.0 * DX) + (V_av * (U.at(i).at(j + 1) - U.at(i).at(j - 1)) - abs(V_av) * (U.at(i).at(j + 1) + U.at(i).at(j - 1) - 2 * U.at(i).at(j))) / (2.0 * DY);
    return U.at(i).at(j) * (U.at(i + 1).at(j) - U.at(i - 1).at(j)) / (2.0 * DX) - abs(U.at(i).at(j)) * (U.at(i + 1).at(j) + U.at(i - 1).at(j) - 2 * U.at(i).at(j)) / (2.0 * DX) + V_av * (U.at(i).at(j + 1) - U.at(i).at(j - 1)) / (2.0 * DY) - abs(V_av) * (U.at(i).at(j + 1) + U.at(i).at(j - 1) - 2 * U.at(i).at(j)) / (2.0 * DY);
}

// y方向対流項
double calc_CNVV(int i, int j)
{
    double U_av = 0.25 * (U.at(i).at(j) + U.at(i - 1).at(j) + U.at(i).at(j + 1) + U.at(i - 1).at(j + 1));
    // return (U_av * (V.at(i + 1).at(j) - V.at(i - 1).at(j)) - abs(U_av) * (V.at(i + 1).at(j) + V.at(i - 1).at(j) - 2 * V.at(i).at(j))) / (2.0 * DX) + (V.at(i).at(j) * (V.at(i).at(j + 1) - V.at(i).at(j - 1)) - abs(V.at(i).at(j)) * (V.at(i).at(j + 1) + V.at(i).at(j - 1) - 2 * V.at(i).at(j))) / (2.0 * DY);
    return U_av * (V.at(i + 1).at(j) - V.at(i - 1).at(j)) / (2.0 * DX) - abs(U_av) * (V.at(i + 1).at(j) + V.at(i - 1).at(j) - 2 * V.at(i).at(j)) / (2.0 * DX) + V.at(i).at(j) * (V.at(i).at(j + 1) - V.at(i).at(j - 1)) / (2.0 * DY) - abs(V.at(i).at(j)) * (V.at(i).at(j + 1) + V.at(i).at(j - 1) - 2 * V.at(i).at(j)) / (2.0 * DY);
}

// 温度対流項
double calc_CNVT(int i, int j)
{
    double U_av = 0.5 * (U.at(i).at(j) + U.at(i - 1).at(j));
    double V_av = 0.5 * (V.at(i).at(j) + V.at(i).at(j - 1));
    // return (U_av * (T.at(i + 1).at(j) - T.at(i - 1).at(j)) - abs(U_av) * (T.at(i + 1).at(j) + T.at(i - 1).at(j) - 2 * T.at(i).at(j))) / (2.0 * DX) + (V_av * (T.at(i).at(j + 1) - T.at(i).at(j - 1)) - abs(V_av) * (T.at(i).at(j + 1) + T.at(i).at(j - 1) - 2 * T.at(i).at(j))) / (2.0 * DY);
    return U_av * (T.at(i + 1).at(j) - T.at(i - 1).at(j)) / (2.0 * DX) - abs(U_av) * (T.at(i + 1).at(j) + T.at(i - 1).at(j) - 2 * T.at(i).at(j)) / (2.0 * DX) + V_av * (T.at(i).at(j + 1) - T.at(i).at(j - 1)) / (2.0 * DY) - abs(V_av) * (T.at(i).at(j + 1) + T.at(i).at(j - 1) - 2 * T.at(i).at(j)) / (2.0 * DY);
}

// x方向拡散項
double calc_DIFU(int i, int j)
{
    return PR * ((U.at(i + 1).at(j) - 2.0 * U.at(i).at(j) + U.at(i - 1).at(j)) / (DX * DX) + (U.at(i).at(j + 1) - 2.0 * U.at(i).at(j) + U.at(i).at(j - 1)) / (DY * DY));
}

// y方向拡散項
double calc_DIFV(int i, int j)
{
    return PR * ((V.at(i + 1).at(j) - 2.0 * V.at(i).at(j) + V.at(i - 1).at(j)) / (DX * DX) + (V.at(i).at(j + 1) - 2.0 * V.at(i).at(j) + V.at(i).at(j - 1)) / (DY * DY));
}

// 温度拡散項
double calc_DIFT(int i, int j)
{
    return PR * ((T.at(i + 1).at(j) - 2.0 * T.at(i).at(j) + T.at(i - 1).at(j)) / (DX * DX) + (T.at(i).at(j + 1) - 2.0 * T.at(i).at(j) + T.at(i).at(j - 1)) / (DY * DY));
}

// 重力項
double calc_BUOV(int i, int j)
{
    return RA * PR * (T.at(i).at(j) + T.at(i).at(j + 1)) / 2.0;
}

// 各項を計算
void calc_terms()
{
    // x方向速度
    for (int i = 1; i < NX; i++)
    {
        for (int j = 1; j < NY + 1; j++)
        {
            CNVU.at(i).at(j) = calc_CNVU(i, j);
            DIFU.at(i).at(j) = calc_DIFU(i, j);
        }
    }
    // y方向速度
    for (int i = 1; i < NX + 1; i++)
    {
        for (int j = 1; j < NY; j++)
        {
            CNVV.at(i).at(j) = calc_CNVV(i, j);
            DIFV.at(i).at(j) = calc_DIFV(i, j);
            BUOV.at(i).at(j) = calc_BUOV(i, j);
        }
    }
    // 温度
    for (int i = 1; i < NX + 1; i++)
    {
        for (int j = 1; j < NY + 1; j++)
        {
            CNVT.at(i).at(j) = calc_CNVT(i, j);
            DIFT.at(i).at(j) = calc_DIFT(i, j);
        }
    }
}

// 仮の速度場を求める
void calc_UV_T()
{
    // x方向速度
    for (int i = 1; i < NX; i++)
    {
        for (int j = 1; j < NY + 1; j++)
        {
            U_T.at(i).at(j) = U.at(i).at(j) + DT * ((P.at(i).at(j) - P.at(i + 1).at(j)) / DX + DIFU.at(i).at(j) - CNVU.at(i).at(j));
        }
    }
    // y方向速度
    for (int i = 1; i < NX + 1; i++)
    {
        for (int j = 1; j < NY; j++)
        {
            V_T.at(i).at(j) = V.at(i).at(j) + DT * ((P.at(i).at(j) - P.at(i).at(j + 1)) / DY + DIFV.at(i).at(j) + BUOV.at(i).at(j) - CNVV.at(i).at(j));
        }
    }
}

inline double poison_succeed(int i, int j)
{
    return (1.0 - OMG) * P_D_new.at(i).at(j) + OMG * (-DIV.at(i).at(j) / DT + (P_D_new.at(i + 1).at(j) + P_D_new.at(i - 1).at(j)) / DX / DX + (P_D_new.at(i).at(j + 1) + P_D_new.at(i).at(j - 1)) / DY / DY) / (2.0 * (1.0 / DX / DX + 1.0 / DY / DY));
}

// ポアソン方程式を解く
void calc_DP()
{
    for (int i = 1; i < NX + 1; i++)
    {
        for (int j = 1; j < NY + 1; j++)
        {
            DIV.at(i).at(j) = (U_T.at(i).at(j) - U_T.at(i - 1).at(j)) / DX + (V_T.at(i).at(j) - V_T.at(i).at(j - 1)) / DY;
        }
    }
    double err = 1.0;
    // P_D, P_D_newを初期化
    memset(&P_D[0][0], 0, sizeof(double) * (NX + 2) * (NY + 2));
    memset(&P_D_new[0][0], 0, sizeof(double) * (NX + 2) * (NY + 2));
    int count = 0;
    while (err > EPSP)
    {
        count++;
        // red-black法
        for (int i = 1; i < NX + 1; i += 2)
        {
            for (int j = 1; j < NY + 1; j += 2)
            {
                P_D_new.at(i).at(j) = poison_succeed(i, j);
            }
        }
        for (int i = 2; i < NX + 1; i += 2)
        {
            for (int j = 2; j < NY + 1; j += 2)
            {
                P_D_new.at(i).at(j) = poison_succeed(i, j);
            }
        }
        for (int i = 2; i < NX + 1; i += 2)
        {
            for (int j = 1; j < NY + 1; j += 2)
            {
                P_D_new.at(i).at(j) = poison_succeed(i, j);
            }
        }
        for (int i = 1; i < NX + 1; i += 2)
        {
            for (int j = 2; j < NY + 1; j += 2)
            {
                P_D_new.at(i).at(j) = poison_succeed(i, j);
            }
        }
        // 収束判定
        err = 0.0;
        // errorは差の二乗の合計
        for (int i = 1; i < NX + 1; i++)
        {
            for (int j = 1; j < NY + 1; j++)
            {
                err += (P_D_new.at(i).at(j) - P_D.at(i).at(j)) * (P_D_new.at(i).at(j) - P_D.at(i).at(j));
            }
        }
        // P_Dを更新
        memcpy(&P_D[0][0], &P_D_new[0][0], sizeof(double) * (NX + 2) * (NY + 2));
        if (count > 100000)
        {
            cout << "count=" << count << endl;
            break;
        }
    }
    if (count > mymax)
    {
        mymax = count;
        cout << "max=" << mymax << endl;
    }
}

// 解をもとに各値を修正
void calc_UVPT()
{
    // x方向速度
    for (int i = 1; i < NX; i++)
    {
        for (int j = 1; j < NY + 1; j++)
        {
            U.at(i).at(j) = U_T.at(i).at(j) - DT * (P_D.at(i + 1).at(j) - P_D.at(i).at(j)) / DX;
        }
    }
    // y方向速度
    for (int i = 1; i < NX + 1; i++)
    {
        for (int j = 1; j < NY; j++)
        {
            V.at(i).at(j) = V_T.at(i).at(j) - DT * (P_D.at(i).at(j + 1) - P_D.at(i).at(j)) / DY;
        }
    }
    // 圧力
    for (int i = 1; i < NX + 1; i++)
    {
        for (int j = 1; j < NY + 1; j++)
        {
            P.at(i).at(j) = P_D.at(i).at(j) + P.at(i).at(j);
        }
    }
    // 温度
    for (int i = 1; i < NX + 1; i++)
    {
        for (int j = 1; j < NY + 1; j++)
        {
            T.at(i).at(j) = T.at(i).at(j) + DT * (DIFT.at(i).at(j) - CNVT.at(i).at(j));
        }
    }
}

// 境界条件
void set_boundary()
{
    // 速度
    // 境界上の速度は0
    for (int i = 0; i < NY + 2; i++)
    {
        U.at(0).at(i) = 0.0;
        U.at(NX).at(i) = 0.0;
    }
    for (int i = 0; i < NX + 2; i++)
    {
        V.at(i).at(0) = 0.0;
        V.at(i).at(NY) = 0.0;
    }
    // 境界上に代表点がないなら平均を0に
    for (int i = 1; i < NX; i++)
    {
        U.at(i).at(0) = -U.at(i).at(1);
        U.at(i).at(NY + 1) = -U.at(i).at(NY);
    }
    for (int i = 1; i < NY; i++)
    {
        V.at(0).at(i) = -V.at(1).at(i);
        V.at(NX + 1).at(i) = -V.at(NX).at(i);
    }
    // 圧力(境界上は0)
    for (int i = 0; i < NY + 2; i++)
    {
        P.at(0).at(i) = 0.0;
        P.at(NX + 1).at(i) = 0.0;
    }
    for (int i = 0; i < NX + 2; i++)
    {
        P.at(i).at(0) = 0.0;
        P.at(i).at(NY + 1) = 0.0;
    }
    // 温度
    // 左側:加熱
    for (int i = 0; i < NY + 2; i++)
    {
        T.at(0).at(i) = 1.0 - T.at(1).at(i);
    }
    // 右側:冷却
    for (int i = 0; i < NY + 2; i++)
    {
        T.at(NX + 1).at(i) = -1.0 - T.at(NX).at(i);
    }
    // 上下の境界は断熱
    for (int i = 1; i < NX + 1; i++)
    {
        T.at(i).at(0) = T.at(i).at(1);
        T.at(i).at(NY + 1) = T.at(i).at(NY);
    }
}

// 結果をファイルに書き込む
void write_result(int cycle)
{
    // ファイルに書き込む
    fwrite(&cycle, sizeof(int), 1, output);
    fwrite(&U[0][0], sizeof(double), (NX + 1) * (NY + 2), output);
    fwrite(&V[0][0], sizeof(double), (NX + 2) * (NY + 1), output);
    fwrite(&P[0][0], sizeof(double), (NX + 2) * (NY + 2), output);
    fwrite(&T[0][0], sizeof(double), (NX + 2) * (NY + 2), output);
}
