#include <bits/stdc++.h>

using namespace std;

class job
{
public:
    // type: それぞれの仕事の種類
    // 1: ALU(1クロック)
    // 2: FPU(2クロック)
    // 3: MEM(2クロック)
    int type;
    // レイテンシ
    int latency;
    // depend: それぞれの仕事の依存関係
    vector<int> depend;
    // priority: それぞれの仕事の優先度
    int priority;
    // 実行可能になる時間
    int start_time;
    // コンストラクタ
    job()
    {
        type = 0;
        priority = 0;
        depend = {};
        start_time = 0;
    };
    void set_latency()
    {
        switch (type)
        {
        case 1:
            latency = 1;
            break;
        case 2:
            latency = 2;
            break;
        case 3:
            latency = 2;
            break;
        default:
            cout << "error: unknown type" << endl;
            exit(1);
        }
    };
    bool delete_depend(int index)
    {
        for (int i = 0; i < depend.size(); i++)
        {
            if (depend[i] == index)
            {
                depend.erase(depend.begin() + i);
                return true;
            }
        }
        return false;
    };
};

tuple<vector<int>, int> time_priority(int n, vector<job> jobs)
{
    // Ready Setを作成
    vector<int> ready_set;
    for (int i = 0; i < n; i++)
    {
        if (jobs[i].depend.size() == 0)
        {
            ready_set.push_back(i);
            // 先行命令がないので優先度は0
            jobs[i].priority = 0;
        }
    }
    // 各資源の使用状況を記録
    vector<int> alu(1);
    vector<int> fpu(2);
    vector<int> mem(2);
    // 空なら-1,使用中ならその仕事の番号を記録
    alu[0] = -1;
    fpu[0] = -1;
    fpu[1] = -1;
    mem[0] = -1;
    mem[1] = -1;
    int time = 0;
    // 終わった順に仕事の番号を記録
    vector<int> done_order;
    // Ready Setが空になるまで繰り返す
    while (ready_set.size() > 0)
    {
        // Ready Setの中で最も優先度が高い仕事を探す
        int max_priority = -1;
        int index;
        for (int i = 0; i < ready_set.size(); i++)
        {
            if (jobs[ready_set[i]].priority > max_priority && jobs[ready_set[i]].start_time <= time)
            {
                index = ready_set[i];
                max_priority = jobs[ready_set[i]].priority;
            }
        }
        if (max_priority != -1)
        {
            // その仕事を実行
            switch (jobs[index].type)
            {
            case 1:
                // ALU
                // ALUが空いていれば実行
                if (alu[0] == -1)
                {
                    alu[0] = index;
                    done_order.push_back(index);
                    // ready_setから削除
                    for (int i = 0; i < ready_set.size(); i++)
                    {
                        if (ready_set[i] == index)
                        {
                            ready_set.erase(ready_set.begin() + i);
                        }
                    }
                }
                break;
            case 2:
                // FPU
                // FPUが空いていれば実行
                if (fpu[0] == -1)
                {
                    fpu[0] = index;
                    done_order.push_back(index);
                    // ready_setから削除
                    for (int i = 0; i < ready_set.size(); i++)
                    {
                        if (ready_set[i] == index)
                        {
                            ready_set.erase(ready_set.begin() + i);
                        }
                    }
                }
                break;
            case 3:
                // MEM
                // MEMが空いていれば実行
                if (mem[0] == -1)
                {
                    mem[0] = index;
                    done_order.push_back(index);
                    // ready_setから削除
                    for (int i = 0; i < ready_set.size(); i++)
                    {
                        if (ready_set[i] == index)
                        {
                            ready_set.erase(ready_set.begin() + i);
                        }
                    }
                }
                break;
            default:
                cout << "error: unknown type" << endl;
                exit(1);
            }
        }
        // 依存グラフの更新
        if (done_order[done_order.size() - 1] == index)
        {
            // 依存グラフの更新
            for (int i = 0; i < n; i++)
            {
                if (jobs[i].delete_depend(index))
                {
                    // 依存していた仕事がなくなったらReady Setに追加
                    if (jobs[i].depend.size() == 0)
                    {
                        ready_set.push_back(i);
                        // 優先度は先行命令のレイテンシ
                        jobs[i].priority = jobs[index].latency + jobs[i].priority;
                        // 実行可能になる時間を更新
                        if (jobs[i].start_time < time + jobs[index].latency)
                        {
                            jobs[i].start_time = time + jobs[index].latency;
                        }
                    }
                }
            }
        }
        // 時間更新
        ++time;
        // 使用中の資源の時間を更新
        fpu[1] = fpu[0];
        mem[1] = mem[0];
        alu[0] = -1;
        fpu[0] = -1;
        mem[0] = -1;
    }
    time += jobs[done_order[done_order.size() - 1]].latency - 1;
    return make_tuple(done_order, time);
}

tuple<vector<int>, int> resource_priority(int n, vector<job> jobs)
{
    // Ready Setを作成
    vector<int> ready_set;
    for (int i = 0; i < n; i++)
    {
        if (jobs[i].depend.size() == 0)
        {
            ready_set.push_back(i);
        }
    }
    // 各資源の使用状況を記録
    vector<int> alu(1);
    vector<int> fpu(2);
    vector<int> mem(2);
    // 空なら-1,使用中ならその仕事の番号を記録
    alu[0] = -1;
    fpu[0] = -1;
    fpu[1] = -1;
    mem[0] = -1;
    mem[1] = -1;
    int time = 0;
    // 終わった順に仕事の番号を記録
    vector<int> done_order;
    // Ready Setが空になるまで繰り返す
    while (ready_set.size() > 0)
    {
        // Ready Setの中で最も優先度が高い仕事を探す
        int max_priority = -1;
        int index;
        for (int i = 0; i < ready_set.size(); i++)
        {
            if (jobs[ready_set[i]].priority > max_priority)
            {
                index = ready_set[i];
                max_priority = jobs[ready_set[i]].priority;
            }
        }
        if (max_priority != -1 && jobs[index].start_time <= time)
        {
            // その仕事を実行
            switch (jobs[index].type)
            {
            case 1:
                // ALU
                // ALUが空いていれば実行
                if (alu[0] == -1)
                {
                    alu[0] = index;
                    done_order.push_back(index);
                    // ready_setから削除
                    for (int i = 0; i < ready_set.size(); i++)
                    {
                        if (ready_set[i] == index)
                        {
                            ready_set.erase(ready_set.begin() + i);
                        }
                    }
                }
                break;
            case 2:
                // FPU
                // FPUが空いていれば実行
                if (fpu[0] == -1)
                {
                    fpu[0] = index;
                    done_order.push_back(index);
                    // ready_setから削除
                    for (int i = 0; i < ready_set.size(); i++)
                    {
                        if (ready_set[i] == index)
                        {
                            ready_set.erase(ready_set.begin() + i);
                        }
                    }
                }
                break;
            case 3:
                // MEM
                // MEMが空いていれば実行
                if (mem[0] == -1)
                {
                    mem[0] = index;
                    done_order.push_back(index);
                    // ready_setから削除
                    for (int i = 0; i < ready_set.size(); i++)
                    {
                        if (ready_set[i] == index)
                        {
                            ready_set.erase(ready_set.begin() + i);
                        }
                    }
                }
                break;
            default:
                cout << "error: unknown type" << endl;
                exit(1);
            }
        }
        // 依存グラフの更新
        if (done_order[done_order.size() - 1] == index)
        {
            // 依存グラフの更新
            for (int i = 0; i < n; i++)
            {
                if (jobs[i].delete_depend(index))
                {
                    // 依存していた仕事がなくなったらReady Setに追加
                    if (jobs[i].depend.size() == 0)
                    {
                        ready_set.push_back(i);
                        // 優先度は先行命令の数
                        jobs[i].priority = max(jobs[i].priority, jobs[index].priority + 1);
                        // 実行可能になる時間を更新
                        if (jobs[i].start_time < time + jobs[index].latency)
                        {
                            jobs[i].start_time = time + jobs[index].latency;
                        }
                    }
                }
            }
        }
        // 時間更新
        ++time;
        // 使用中の資源の時間を更新
        fpu[1] = fpu[0];
        mem[1] = mem[0];
        alu[0] = -1;
        fpu[0] = -1;
        mem[0] = -1;
    }
    time += jobs[done_order[done_order.size() - 1]].latency - 1;
    return make_tuple(done_order, time);
}

int main()
{
    // 条件入力
    // n: 仕事の数
    // m: 仕事の依存関係の数
    int n, m;
    cin >> n >> m;
    vector<job> jobs(n);
    // jobsのtypeを入力
    for (int i = 0; i < n; i++)
    {
        cin >> jobs[i].type;
    }
    // jobsのdependを入力
    for (int i = 0; i < m; i++)
    {
        int x, y;
        cin >> x >> y;
        jobs[y].depend.push_back(x);
    }
    // jobsのレイテンシを入力
    for (int i = 0; i < n; i++)
    {
        jobs[i].set_latency();
    }
    int time;
    vector<int> done_order;
    tie(done_order, time) = time_priority(n, jobs);
    cout << "time_priority" << endl;
    cout << "   time: " << time << endl;
    cout << "   done_order: ";
    for (int i = 0; i < done_order.size(); i++)
    {
        cout << done_order[i] << " ";
    }
    cout << endl;
    tie(done_order, time) = resource_priority(n, jobs);
    cout << "resource_priority" << endl;
    cout << "   time: " << time << endl;
    cout << "   done_order: ";
    for (int i = 0; i < done_order.size(); i++)
    {
        cout << done_order[i] << " ";
    }
    cout << endl;
    return 0;
}
