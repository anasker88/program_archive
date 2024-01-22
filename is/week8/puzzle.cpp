#include <bits/stdc++.h>

using namespace std;

vector<pair<int, int>> dir = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}};

int n, m;
int board_num_end;
pair<bool, stack<int>>
solve(vector<vector<int>> a, vector<vector<int>> end, int zero_x, int zero_y, int last, int depth, unordered_map<int, int> &checked);

int board_to_num(vector<vector<int>> &a);

int main()
{
    int zero_x = 0, zero_y = 0;
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
    board_num_end = board_to_num(end);
    int already_checked = 0;
    unordered_map<int, int> checked = {};
    for (int depth = 3; true; depth += 3)
    {
        auto [ok, s] = solve(a, end, zero_x, zero_y, -1, depth, checked);
        if (ok)
        {
            cout << s.size() << "steps" << endl
                 << endl;
            for (auto &row : a)
            {
                for (auto &x : row)
                    cout << setw(2) << x << " ";
                cout << endl;
            }
            cout << endl;
            int i = 0;
            while (!s.empty())
            {
                i = s.top();
                // cout << i << endl;
                swap(a[zero_x][zero_y], a[zero_x + dir[i].first][zero_y + dir[i].second]);
                zero_x += dir[i].first;
                zero_y += dir[i].second;
                for (auto &row : a)
                {
                    for (auto &x : row)
                        cout << setw(2) << x << " ";
                    cout << endl;
                }
                cout << endl;
                s.pop();
            }
            return 0;
        }
        cout << "depth " << depth << " failed" << endl;
        // checkedのtrueの数を数える
        int count = 0;
        for (auto &p : checked)
            if (p.second)
                count++;
        cout << "checked: " << count << endl;
        if (count == already_checked)
        {
            cout << "no answer" << endl;
            break;
        }
        already_checked = count;
    }
    return 0;
}

int board_to_num(vector<vector<int>> &a)
{
    // boardを(n*m)!-1以下の数値に変換する
    int num = 0;
    int fact = n * m;
    unordered_map<int, bool> used;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            int cur = a[i][j];
            if (used[cur])
            {
                cout << "include same number" << endl;
                exit(1);
            }
            used[cur] = true;
            int index = cur;
            for (int k = 0; k < cur; k++)
            {
                if (used[k])
                {
                    index--;
                }
            }
            num *= fact;
            num += index;
            fact--;
        }
    }
    return num;
}

pair<bool, stack<int>> solve(vector<vector<int>> a, vector<vector<int>> end, int zero_x, int zero_y, int last, int depth, unordered_map<int, int> &checked)
{
    stack<int> s = {};
    int board_num = board_to_num(a);
    if (checked[board_num] < depth)
        checked[board_num] = depth;
    else
        return {false, s};
    if (board_num == board_num_end)
    {
        return {true, s};
    }
    if (depth == 0)
    {
        return {false, s};
    }
    for (int i = 0; i < 4; i++)
    {
        // cout << i << endl;
        pair<int, int> d = dir[i];
        int x = zero_x + d.first;
        int y = zero_y + d.second;
        vector<vector<int>> tmp_a = a;
        int tmp_zero_x = zero_x;
        int tmp_zero_y = zero_y;
        if (x >= 0 && x < n && y >= 0 && y < m)
        {
            tmp_a = a;
            swap(tmp_a[zero_x][zero_y], tmp_a[x][y]);
            tmp_zero_x = x;
            tmp_zero_y = y;
            auto [ok, s] = solve(tmp_a, end, tmp_zero_x, tmp_zero_y, i, depth - 1, checked);
            if (ok)
            {
                s.push(i);
                return {true, s};
            }
        }
    }
    return {false, s};
}
