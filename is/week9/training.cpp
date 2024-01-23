#include <bits/stdc++.h>

using namespace std;

int main()
{
    // パラメータの最尤推定
    map<pair<string, string>, int> transition;
    map<pair<string, string>, int> emission;
    map<string, int> word_count;

    ifstream ifs("wsj_sample_train.pos");
    string line;
    string buf;
    while (ifs && getline(ifs, line))
    {
        istringstream iss(line);
        string prev = "BOS";
        while (iss >> buf)
        {
            // buf: "word/POS"
            // wordは最後の'/'より前の文字列
            // POSは最後の'/'より後の文字列
            int pos = buf.rfind('/');
            string word = buf.substr(0, pos);
            string POS = buf.substr(pos + 1);
            // wordを小文字に変換
            transform(word.begin(), word.end(), word.begin(), ::tolower);
            // transitionのカウント
            transition[make_pair(prev, POS)]++;
            // emissionのカウント
            emission[make_pair(POS, word)]++;
            // word_countのカウント
            word_count[word]++;
            prev = POS;
        }
        // EOSのカウント
        transition[make_pair(prev, "EOS")]++;
    }
    ifs.close();

    // 未知語の変換
    for (auto itr = emission.begin(); itr != emission.end(); itr++)
    {
        if (word_count[itr->first.second] <= 1)
        {
            emission[make_pair(itr->first.first, "UNK")] += itr->second;
            emission.erase(itr);
        }
    }
    // transitionとemissionのカウント
    unordered_map<string, int> tran_count;
    unordered_map<string, int> emi_count;
    for (auto itr = transition.begin(); itr != transition.end(); itr++)
    {
        tran_count[itr->first.first] += itr->second;
    }
    for (auto itr = emission.begin(); itr != emission.end(); itr++)
    {
        emi_count[itr->first.first] += itr->second;
    }
    // a,bを推定
    map<pair<string, string>, double> a;
    map<pair<string, string>, double> b;
    for (auto itr = transition.begin(); itr != transition.end(); itr++)
    {
        a[itr->first] = (double)itr->second / tran_count[itr->first.first];
    }
    for (auto itr = emission.begin(); itr != emission.end(); itr++)
    {
        b[itr->first] = (double)itr->second / emi_count[itr->first.first];
    }
    // a,bの出力
    ofstream ofs1("a.txt");
    for (auto itr = a.begin(); itr != a.end(); itr++)
    {
        ofs1 << itr->first.first << " " << itr->first.second << " " << itr->second << endl;
    }
    ofstream ofs2("b.txt");
    for (auto itr = b.begin(); itr != b.end(); itr++)
    {
        ofs2 << itr->first.first << " " << itr->first.second << " " << itr->second << endl;
    }
    // テストデータの予測
    ifstream ifs2("wsj_sample_test.pos");
    ofstream ofs("result.pos");
    while (ifs2 && getline(ifs2, line))
    {
        istringstream iss(line);
        while (iss >> buf)
        {
        }
    }

    return 0;
}
