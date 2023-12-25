#include <bits/stdc++.h>

using namespace std;

int main()
{
    int n, m;
    int zero_x, zero_y;
    cin >> n >> m;
    vector<vector<int>> a(n, vector<int>(m));
    for (auto &row : a)
        for (auto &x : row)
        {
            cin >> x;
            if (x == 0)
            {
                zero_x = &row - &a[0];
                zero_y = &x - &row[0];
            }
        }
    vector<vector<int>> end(n, vector<int>(m));
    int num = 1;
    for (auto &row : end)
        for (auto &x : row)
            x = num++;
    end[n - 1][m - 1] = 0;

    return 0;
}
