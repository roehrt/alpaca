#include <bits/stdc++.h>
#include "alpaca.cpp"
using namespace std;
#ifndef NDEBUG
extern "C" int __lsan_is_turned_off() { return 1; }
#endif

int sigints_received = 0;
void handle_sigint(int) { if (++sigints_received == 3) exit(SIGINT); }

int main() {
    signal(SIGINT, handle_sigint);
    srand(0);

    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    while (cin.peek() == 'c') cin.ignore(numeric_limits<streamsize>::max(), '\n');
    string p, task;
    cin >> p >> task;
    if (p != "p") exit(42);

    int n;
    vector<vector<int32_t> > clauses;

    if (task == "ds") {
        int m;
        cin >> n >> m;
        clauses.resize(n);
        for (int i = 1; i <= n; ++i) clauses[i - 1].push_back(i);
        for (int i = 0; i < m; ++i) {
            int u, v;
            cin >> u >> v;
            clauses[u - 1].push_back(v);
            clauses[v - 1].push_back(u);
        }
    } else if (task == "hs") {
        int m;
        cin >> n >> m;
        clauses.resize(m);
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
        for (auto &clause: clauses) {
            string line;
            getline(cin, line);
            istringstream iss(line);
            for (int x; iss >> x;) clause.push_back(x);
            assert(!clause.empty());
        }
    }

    Alpaca::Settings settings;
    // settings.should_terminate = [&]() -> bool { return sigints_received; };
    settings.verbosity = 0;

    if (task == "hs") {
        size_t total_lits = 0;
        for (const auto &clause: clauses) total_lits += clause.size();
        if (2 * total_lits < 7 * clauses.size()) settings.preprocessing_techniques = "[bus]#[[[[vu]b]sr]lc]";
        settings.ilp_enabled = true;
    }

    Alpaca::Solver solver(settings);
    for (const auto &clause: clauses)
        solver.add(clause);
    solver.solve();
    cout << solver.objective_value() << "\n";
    for (int i = 1; i <= n; ++i)
        if (solver.value(i) > 0)
            cout << i << "\n";

    return 0;
}
