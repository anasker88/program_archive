#include <stdio.h>
#include <omp.h>
#include <vector>
#include <iostream>
#include <stdio.h>

using namespace std;

// 単純に２乗を計算して返す関数
double test_square(double x)
{
    return x * x;
}

int main()
{
    printf("使用可能な最大スレッド数：%d\n", omp_get_max_threads());

    double test;

    // 最後に結果を格納する配列
    vector<double> box;
    box.resize(10);

    // 各for文で呼ばれる配列
    vector<double> chocolates;
    for (int i = 0; i < 10; i++)
    {
        chocolates.push_back(i);
    }

// for文だよ〜
// chocolateの各要素を２乗してboxの配列に格納するよ〜
#pragma omp parallel for
    for (int i = 0; i < 10; i++)
    {

        test = test_square(chocolates[i]);
        box[i] = test;
        printf("thread = %d, i = %2d\n", omp_get_thread_num(), i);
    }

    // 最後に結果をプリントするよ〜
    for (int i = 0; i < 10; i++)
    {
        cout << box[i] << endl;
    }
    return 0;
}
