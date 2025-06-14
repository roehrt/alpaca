#include <bits/extc++.h>
#include <fmt/format.h>
#include "scippp/model.hpp"
#include "preprocessorinterface.hpp"
#include "cadical.hpp"

extern "C" {
#include "kissat.h"
}
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "objscip/objscip.h"

namespace Alpaca {
    using namespace std;

    namespace SetPacking {
        void pack(vector<vector<int32_t> > &sets) {
            if (sets.size() >= 100000) {
                ranges::sort(sets, {}, [](const vector<int32_t> &s) { return s.size(); });
                unordered_set<int32_t> seen;
                erase_if(sets, [&](const vector<int32_t>& s) {
                    for (int32_t x : s) if (seen.contains(x)) return true;
                    for (int32_t x : s) seen.insert(x);
                    return false;
                });
                return;
            }
            unordered_map<int32_t, vector<int32_t> > literal_to_set;
            for (size_t i = 0; i < sets.size(); ++i)
                for (int32_t lit: sets[i])
                    literal_to_set[abs(lit)].emplace_back(i);
            vector g(sets.size(), vector(sets.size(), false));
            for (auto const &s: views::values(literal_to_set))
                for (auto i : s)
                    for (auto j : s)
                        g[i][j] = i != j;

            vector<size_t> degree(sets.size(), 0);
            for (size_t i = 0; i < sets.size(); ++i) degree[i] = ranges::count(g[i], true);
            vector blocked(sets.size(), false);
            for (size_t i = 0; i < sets.size(); ++i) {
                size_t u = ranges::min_element(degree) - begin(degree);
                degree[u] = -1;
                if (blocked[u]) {
                    sets[u].clear();
                    continue;
                }
                for (size_t v = 0; v < sets.size(); ++v)
                    if (g[u][v]) {
                        --degree[v];
                        blocked[v] = true;
                    }
            }
            erase_if(sets, [](auto &s) { return s.empty(); });
        }
    }

    enum Status {
        ERROR = 0,
        UNKNOWN = 1,
        UNSATISFIABLE = 2,
        SATISFIABLE = 3,
        OPTIMAL = 4,
    };

    struct Settings {
        double preprocessing_seconds = 120;
        string preprocessing_techniques = "[bus]#[[[[[[vu]b]sr]lc]G]ea]";
        int sat_coreminimization_budget = 5;
        int sat_coreminimization_conflicts = 1000;
        int sat_core_rounds = 8;
        bool ilp_enabled = false;
        int verbosity = 1;
        uint64_t seed = 0;
        function<bool()> should_terminate = [] {
            return false;
        };
    };

    struct SolverContext {
        Settings settings;

        int32_t vars = 0;
        vector<int32_t> soft;
        vector<vector<int32_t> > clauses;

        uint64_t offset = 0;
        uint64_t bound = 0;
        uint64_t objective_value = numeric_limits<uint64_t>::max();
        vector<int32_t> model;

        bool presolve = true;
        mt19937 rng;
        chrono::steady_clock::time_point start;

        explicit SolverContext(const Settings &settings) : settings(settings), rng(settings.seed),
                                                           start(chrono::steady_clock::now()) {
        }

        [[nodiscard]] string prefix() const {
            auto elapsed_seconds = static_cast<double>(chrono::duration_cast<chrono::milliseconds>(
                                       chrono::steady_clock::now() - start).count()) / 1000.0;
            if (presolve)
                return fmt::format(
                    "{:.1f}s [presolve] ",
                    elapsed_seconds
                );
            string objective_string = objective_value == numeric_limits<uint64_t>::max()
                                          ? "inf"
                                          : to_string(offset + objective_value);
            return fmt::format(
                "{:.1f}s [{}-{}] ",
                elapsed_seconds,
                offset + bound,
                objective_string
            );
        }

        void log(const string &msg, int level = 2) const {
            if (settings.verbosity < level) return;
            cerr << prefix() << msg << "\n";
        }
    };

    class SCIPSolutionHandler : public scip::ObjEventhdlr {
    public:
        SolverContext &context;

        SCIPSolutionHandler(SCIP *scip, SolverContext &context) : ObjEventhdlr(scip, "solutions",
                                                                               "pass solutions to alpaca"),
                                                                  context(context) {
        }

        ~SCIPSolutionHandler() override = default;

        SCIP_DECL_EVENTINITSOL(scip_initsol) override {
            SCIP_CALL(SCIPcatchEvent( scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, nullptr, nullptr));
            return SCIP_OKAY;
        }

        SCIP_DECL_EVENTEXEC(scip_exec) override {
            SCIP_SOL *sol = SCIPeventGetSol(event);
            if (sol == nullptr) return SCIP_OKAY;
            auto obj = llround(SCIPgetSolOrigObj(scip, sol)) - context.offset;
            if (obj < context.objective_value) {
                context.objective_value = obj;
                SCIP_VAR **vars = SCIPgetVars(scip);
                for (int32_t i = 1; i <= context.vars; ++i)
                    context.model[i] = SCIPisZero(scip, SCIPgetSolVal(scip, sol, vars[i - 1])) ? -i : i;
                context.log("SCIP improved model");
            }
            return SCIP_OKAY;
        }
    };

    class SCIPBoundHandler : public scip::ObjEventhdlr {
    public:
        SolverContext &context;

        SCIPBoundHandler(SCIP *scip, SolverContext &context) : ObjEventhdlr(scip, "bounds",
                                                                            "pass bounds to alpaca"),
                                                               context(context) {
        }

        ~SCIPBoundHandler() override = default;

        SCIP_DECL_EVENTINIT(scip_init) override {
            SCIP_CALL(SCIPcatchEvent( scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, nullptr, nullptr));
            return SCIP_OKAY;
        }

        SCIP_DECL_EVENTEXEC(scip_exec) override {
            auto bound = llround(SCIPceil(scip, SCIPgetDualbound(scip))) - context.offset;
            if (bound > context.bound) {
                context.bound = bound;
                context.log("SCIP improved bound");
            }
            return SCIP_OKAY;
        }
    };

    class SCIPTerminator : public scip::ObjEventhdlr {
    public:
        SolverContext &context;

        SCIPTerminator(SCIP *scip, SolverContext &context) : ObjEventhdlr(scip, "terminator",
                                                                          "terminates SCIP solve"),
                                                             context(context) {
        }

        ~SCIPTerminator() override = default;

        SCIP_DECL_EVENTINIT(scip_init) override {
            SCIP_CALL(SCIPcatchEvent( scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, nullptr, nullptr));
            return SCIP_OKAY;
        }

        SCIP_DECL_EVENTEXEC(scip_exec) override {
            if (context.settings.should_terminate())
                SCIP_CALL(SCIPinterruptSolve(scip));
            return SCIP_OKAY;
        }
    };

    class SCIPMessageHandler : public scip::ObjMessagehdlr {
        SolverContext &context;

    public:
        SCIPMessageHandler(SolverContext &context, bool bufferedoutput) : ObjMessagehdlr(bufferedoutput),
                                                                          context(context) {
        }

        ~SCIPMessageHandler() override = default;

        static void trim_right(string &s) {
            while (!s.empty() && isspace(s.back()))
                s.pop_back();
        }

        SCIP_DECL_MESSAGEWARNING(scip_warning) override {
            string s = msg;
            trim_right(s);
            if (!empty(s)) context.log(s, 2);
        }

        SCIP_DECL_MESSAGEDIALOG(scip_dialog) override {
            string s = msg;
            trim_right(s);
            if (!empty(s)) context.log(s, 3);
        }

        SCIP_DECL_MESSAGEINFO(scip_info) override {
            string s = msg;
            trim_right(s);
            if (!empty(s)) context.log(s, 3);
        }
    };

    struct CaDiCaLTerminator : CaDiCaL::Terminator {
        SolverContext& context;
        chrono::steady_clock::time_point start;
        explicit CaDiCaLTerminator(SolverContext& context) : context(context), start(chrono::steady_clock::now()) {}
        bool terminate () override {
            if (context.settings.ilp_enabled && chrono::steady_clock::now() - start > chrono::seconds(30)) return true;
            return context.settings.should_terminate();
        }
    };

    class PrefixBuf : public streambuf {
    public:
        PrefixBuf(streambuf* buf, SolverContext& context)
            : dest(buf), context(context), atLineStart(true) {

        }

    protected:
        int overflow(int ch) override {
            if (ch == EOF) {
                return !EOF;
            }

            if (atLineStart && ch != '\n') {
                string prefix = context.prefix();
                dest->sputn(prefix.c_str(), prefix.size());
            }

            atLineStart = (ch == '\n');
            return dest->sputc(ch);
        }

        int sync() override {
            return dest->pubsync();
        }

    private:
        streambuf* dest;
        SolverContext& context;
        bool atLineStart;
    };

    struct Preprocessor {
        static constexpr uint64_t TOP_WEIGHT = numeric_limits<uint64_t>::max();
        SolverContext &context;
        unique_ptr<maxPreprocessor::PreprocessorInterface> maxpre;

        explicit Preprocessor(SolverContext &context) : context(context) {
        }

        vector<int32_t> labels;
        vector<vector<int32_t> > clauses;
        vector<uint64_t> weights;

        int32_t original_size = 0; // TODO check if there is a better/simpler solution
        vector<int32_t> decompress;

        void compress() {
            for (const auto &clause: context.clauses) for (auto lit: clause) decompress.push_back(abs(lit));
            ranges::sort(decompress);
            decompress.erase(ranges::unique(decompress).begin(), decompress.end());
            decompress.insert(decompress.begin(), 0);
            for (auto &clause: context.clauses)
                for (auto &lit: clause) {
                    auto var = ranges::lower_bound(decompress, abs(lit)) - begin(decompress);
                    lit = lit > 0 ? var : -var;
                }
            for (auto &i: context.soft) {
                auto var = ranges::lower_bound(decompress, abs(i)) - begin(decompress);
                i = i > 0 ? var : -var;
            }
            context.vars = static_cast<int32_t>(decompress.size() - 1);
        }

        void preprocess(const double time_limit = numeric_limits<double>::infinity()) {
            original_size = context.vars;
            for (int32_t i = 1; i <= context.vars; ++i) clauses.push_back({-i}), weights.push_back(1);
            for (const auto &clause: context.clauses) clauses.push_back(clause), weights.push_back(TOP_WEIGHT);
            maxpre = make_unique<maxPreprocessor::PreprocessorInterface>(clauses, weights, TOP_WEIGHT, false);
            string techniques = context.settings.preprocessing_techniques;
            context.log(fmt::format("Running maxpre with techniques {}", techniques));
            context.soft.clear();
            streambuf* originalBuf = cerr.rdbuf();
            PrefixBuf prefixBuf(cerr.rdbuf(), context);
            cerr.rdbuf(&prefixBuf);
            maxpre->preprocess(techniques, context.settings.verbosity - 1, time_limit);
            cerr.rdbuf(originalBuf);
            maxpre->getInstance(clauses, weights, labels, false, false);
            vector<vector<int32_t> > preprocessed_collection;
            for (size_t i = 0; i < clauses.size(); ++i) {
                const auto &clause = clauses[i];
                const auto &weight = weights[i];
                if (weight == TOP_WEIGHT) {
                    preprocessed_collection.push_back(clause);
                } else {
                    assert(weight == 1 && clause.size() == 1);
                    context.soft.push_back(clause[0]);
                }
            }
            context.offset = maxpre->getRemovedWeight()[0];
            if (maxpre->getUpperBound() < TOP_WEIGHT) {
                context.log(fmt::format("bound {}", maxpre->getUpperBound() - context.offset));
                // context.model = maxpre->bestModel;
                // TODO use this
            }
            context.clauses = preprocessed_collection;
            context.log(fmt::format("soft {}", context.soft.size()));
            compress();
            context.model.resize(context.vars + 1);
            context.presolve = false;
            context.log("maxpre presolve");
        }

        vector<int32_t> reconstruct(const vector<int32_t> &model) {
            vector<int32_t> internal_model(original_size + 1);
            for (int32_t i = 1; i < internal_model.size(); ++i) internal_model[i] = -i;
            for (int32_t i: model) {
                if (i == 0) continue;
                int32_t j = decompress[abs(i)];
                if (j > original_size) continue; // BVA variables
                internal_model[j] = i > 0 ? j : -j;
            }
            erase(internal_model, 0);
            return maxpre->reconstruct(internal_model);
        }
    };

    struct Export final : CaDiCaL::ClauseIterator {
        function<bool(const vector<int> &)> fn;

        explicit Export(const function<bool(const vector<int> &)> &fn)
            : fn(fn) {
        }

        bool clause(const vector<int> &c) override {
            return fn(c);
        }
    };

    struct Internal {
        // TODO state machine
        // TODO factor out ilp into own class

        struct Node {
            shared_ptr<Node> left, right;
            size_t size;
            vector<int32_t> cnt, children;

            explicit Node(const vector<int32_t> &children) : children(children) {
                if (children.size() == 1) {
                    left = nullptr, right = nullptr;
                    size = 1;
                    cnt = {children[0]};
                    return;
                }
                auto mid = children.begin() + ssize(children) / 2;
                left = make_shared<Node>(vector(children.begin(), mid));
                right = make_shared<Node>(vector(mid, children.end()));
                size = left->size + right->size;
            }
        };

        void extend(const shared_ptr<Node> &node) {
            if (node->cnt.size() >= node->size) return;
            extend(node->left);
            extend(node->right);
            const size_t k = node->cnt.size();
            node->cnt.push_back(new_var());
            if (k < node->left->cnt.size()) oracle.clause(-node->left->cnt[k], node->cnt[k]);
            if (k < node->right->cnt.size()) oracle.clause(-node->right->cnt[k], node->cnt[k]);
            for (size_t i = 0; i < min(k, node->left->cnt.size()); ++i) {
                size_t j = k - (i + 1);
                if (j < node->right->cnt.size())
                    oracle.clause(-node->left->cnt[i], -node->right->cnt[j], node->cnt[k]);
            }
        }

        void relax(const shared_ptr<Node> &node) {
            if (node->cnt.size() >= node->size) return;
            extend(node);
            assumptions.insert(-node->cnt.back());
            assumptions.erase(-end(node->cnt)[-2]);
            depth[abs(node->cnt.back())] = depth[abs(end(node->cnt)[-2])];
            forest[-node->cnt.back()] = node;
            forest.erase(-end(node->cnt)[-2]);
        }

        SolverContext &context;
        CaDiCaL::Solver oracle;
        set<int32_t> assumptions;
        map<int32_t, shared_ptr<Node> > forest;
        map<int32_t, int64_t > depth;

        explicit Internal(SolverContext &context) : context(context) {
            for (const auto &clause: context.clauses) oracle.clause(clause);
        }

        void handle_model() {
            if (const uint64_t value = ranges::count_if(context.soft,
                                                        [&](const int32_t i) {
                                                            return (oracle.val(abs(i)) > 0) != (i > 0);
                                                        });
                value < context.objective_value) {
                context.objective_value = value;
                for (int32_t i = 1; i <= context.vars; ++i)
                    context.model[abs(i)] = (oracle.val(abs(i)) > 0) == (i > 0) ? i : -i;
                context.log("model improved");
            }
        }

        int32_t new_var() {
            static int32_t var = 0;
            return var = max(var, oracle.vars()) + 1;
        }

        void handle_core(const vector<int32_t> &core) {
            int64_t max_depth = 0;
            for (int32_t lit: core) {
                assumptions.erase(lit);
                max_depth = max(max_depth, heu(lit));
                if (forest.contains(lit)) relax(forest[lit]);
            }
            vector<int32_t> relaxed = core;
            for (auto &lit : relaxed) lit = -lit;
            ranges::shuffle(relaxed, context.rng);
            auto node = make_shared<Node>(relaxed);
            extend(node);
            depth[abs(node->cnt.back())] = max_depth + node->size;
            relax(node);
            oracle.clause(node->cnt[0]);
        }

        CaDiCaL::Status sat_solve(const auto &assumed) {
            for (const int32_t lit: assumed) if (lit) oracle.assume(lit);
            auto status = oracle.solve();
            if (status == CaDiCaL::SATISFIABLE) handle_model();
            return static_cast<CaDiCaL::Status>(status);
        }

        int64_t heu(int32_t lit) {
            assert(depth.contains(abs(lit)));
            return depth[abs(lit)];
        }

        void minimize(vector<int32_t> &core) {
            ranges::shuffle(core, context.rng);
            ranges::stable_sort(core, {}, [&](int32_t u) { return -heu(u); });
            int budget = context.settings.sat_coreminimization_budget;
            for (int32_t &u: core) {
                oracle.limit("conflicts", context.settings.sat_coreminimization_conflicts);
                auto core_without_u = core;
                erase(core_without_u, u);
                if (sat_solve(core_without_u) == CaDiCaL::UNSATISFIABLE) u = 0;
                else if (--budget == 0) break;
            }
            erase(core, 0);
        }

        size_t aggressiveness = 1;

        vector<vector<int> > find_cores() {
            vector<vector<int> > cores;
            for (int r = 0; r < context.settings.sat_core_rounds || cores.empty(); ++r) {
                set<int> assumed = assumptions;
                for (int _i = 0; _i < 8; ++_i) {
                    if (context.settings.should_terminate()) return cores;
                    vector a(begin(assumed), end(assumed));
                    ranges::shuffle(a, context.rng);
                    a.resize(a.size() - min(aggressiveness, a.size()));
                    if (sat_solve(a) == CaDiCaL::UNSATISFIABLE) {
                        vector<int> core;
                        for (int lit: a) if (oracle.failed(lit)) core.push_back(lit);
                        minimize(core);
                        for (int lit: core) assumed.erase(lit);
                        cores.emplace_back(core);
                        if (aggressiveness >= 5) aggressiveness += 2;
                    } else {
                        aggressiveness *= 3;
                        aggressiveness /= 4;
                        break;
                    }
                }
            }
            return cores;
        }

        [[nodiscard]] Status zero_weight() const {
            kissat* solver = kissat_init();
            for (const auto &clause: context.clauses) {
                for (int32_t lit : clause) kissat_add(solver, lit);
                kissat_add(solver, 0);
            }
            for (int32_t lit: context.soft) {
                kissat_add(solver, lit);
                kissat_add(solver, 0);
            }
            auto status = kissat_solve(solver);
            if (status == CaDiCaL::SATISFIABLE) {
                context.objective_value = 0;
                context.model.resize(context.vars + 1);
                for (int32_t i = 1; i <= context.vars; ++i)
                    context.model[i] = kissat_value(solver, i) > 0 ? i : -i;
                context.log("trivial sat");
                return OPTIMAL;
            }
            return UNKNOWN;
        }

        void packing_bound() {
            auto trivial_cores = context.clauses;
            erase_if(trivial_cores, [&](const auto &clause) {
                return !ranges::all_of(clause, [&](const int32_t lit) {
                    return assumptions.contains(-lit);
                });
            });
            SetPacking::pack(trivial_cores);
            for (auto &clause: trivial_cores) {
                ranges::transform(clause, clause.begin(), negate());
                handle_core(clause);
                ++context.bound;
            }
            context.log("set packing bound");
        }

        double score(const vector<int32_t> &core) {
            double sum = core.size();
            for (int32_t lit: core) sum += pow(heu(lit), 1.4);
            return sum;
        }

        SCIP_RETCODE runilp() {
            using namespace scippp;
            Model model("model");

            vector<int> objective(context.vars);
            uint64_t offset = context.offset;
            for (int32_t lit: context.soft) {
                objective[abs(lit) - 1] = lit > 0 ? -1 : 1;
                offset += lit > 0;
            }
            auto x = model.addVars("x_", context.vars, objective, VarType::BINARY);
            auto offset_var = model.addVar("offset", offset, VarType::BINARY, 1, 1);

            for (const vector<int> &c : context.clauses) {
                LinExpr clause;
                for (int lit: c)
                    clause += lit > 0 ? x[abs(lit) - 1] : (1 - x[abs(lit) - 1]);
                model.addConstr(clause >= 1, "c");
            }

            auto scip = model.scip();
            SCIP_SOL *sol;
            SCIP_CALL(SCIPcreatePartialSol(scip, &sol, nullptr));
            for (int32_t i = 1; i <= context.vars; ++i)
                SCIP_CALL(SCIPsetSolVal(scip, sol, x[i-1].getVar(), context.model[i] > 0));
            SCIP_CALL(SCIPsetSolVal(scip, sol, offset_var.getVar(), 1));
            SCIP_Bool stored = false;
            SCIP_CALL(SCIPaddSolFree(scip, &sol, &stored));

            model.setParam(params::Param<bool>("misc/catchctrlc"), false);
            model.setParam(params::Param<char>("estimation/restarts/restartpolicy"), 'c');
            model.setParam(params::Param<bool>("separating/filtercutpoolrel"), true);
            SCIP_CALL(SCIPsetHeuristics(scip, SCIP_PARAMSETTING_AGGRESSIVE, true));
            SCIP_CALL(SCIPincludeObjEventhdlr(scip, new SCIPSolutionHandler(scip, context), true));
            SCIP_CALL(SCIPincludeObjEventhdlr(scip, new SCIPBoundHandler(scip, context), true));
            SCIP_CALL(SCIPincludeObjEventhdlr(scip, new SCIPTerminator(scip, context), true));
            SCIP_MESSAGEHDLR *messagehdlr = nullptr;
            SCIP_CALL(SCIPcreateObjMessagehdlr(&messagehdlr, new SCIPMessageHandler(context, true), true));
            SCIP_CALL(SCIPsetMessagehdlr(scip, messagehdlr));
            model.setObjsense(Sense::MINIMIZE);
            model.solve();
            return SCIP_OKAY;
        }

        Status solve() {
            if (zero_weight() == OPTIMAL) return OPTIMAL;
            if (context.settings.should_terminate()) return UNKNOWN;
            size_t sum = 0;
            for (const auto &clause: context.clauses) sum += clause.size();
            context.settings.ilp_enabled &= 2 * sum < 7 * context.clauses.size();
            oracle.connect_terminator(new CaDiCaLTerminator(context));
            if (sat_solve(vector<int32_t>{}) == CaDiCaL::UNSATISFIABLE) {
                context.log("unsatisfiable");
                return UNSATISFIABLE;
            }
            if (context.settings.should_terminate()) return SATISFIABLE;
            for (int32_t lit: context.soft) depth[abs(lit)] = 1;
            assumptions = set(begin(context.soft), end(context.soft));
            packing_bound();
            // deque<double> history;
            // for (int i = 0; i < 10; ++i) history.push_back(50);
            while (sat_solve(assumptions) == CaDiCaL::UNSATISFIABLE) {
                ++context.bound;
                context.log("unsat bound");

                auto cores = find_cores();
                if (context.settings.should_terminate()) return SATISFIABLE;
                if (context.bound == context.objective_value) break;
                auto all_cores = cores;
                SetPacking::pack(all_cores);
                if (all_cores.size() - 1 + context.bound == context.objective_value) {
                    context.bound = context.objective_value;
                    return OPTIMAL;
                }

                ranges::sort(cores, {}, [&](const auto &core) { return score(core); });
                // double avg = accumulate(begin(history), end(history), 0.0) / history.size();

                auto accepted_cores = cores;
                erase_if(accepted_cores, [&](const auto &core) { return score(cores[0]) * 4 < score(core) * 3; });
                /*vector<vector<int>> accepted_cores = {cores[0]};
                for (const auto& core : cores | views::drop(1))
                    if (score(core) < avg * context.settings.sat_coverthreshold)
                        accepted_cores.emplace_back(core);*/
                // history.pop_front();
                // history.push_back(score(cores[0]));
                SetPacking::pack(accepted_cores);

                for (const auto &core: accepted_cores) handle_core(core);

                context.bound += accepted_cores.size() - 1;
                if (accepted_cores.size() > 1) context.log("disjoint cores bound");
                context.log(fmt::format("core quality {:.2f}", score(accepted_cores.back())));

                if (context.settings.should_terminate()) return SATISFIABLE;
                if (context.bound + 1 == context.objective_value) {
                    context.log("hardening");
                    for (const int32_t a: assumptions) oracle.clause(a);
                    if (sat_solve(vector<int32_t>{}) == CaDiCaL::UNSATISFIABLE)
                        context.bound = context.objective_value;
                    return OPTIMAL;
                }
            }
            if (context.bound != context.objective_value && context.settings.ilp_enabled) {
                if (runilp() != SCIP_OKAY) return ERROR;
                return context.bound == context.objective_value ? OPTIMAL : SATISFIABLE;
            }
            if (context.settings.should_terminate()) return SATISFIABLE;
            assert(context.bound == context.objective_value);
            return OPTIMAL;
        }
    };

    class Solver {
        SolverContext context;
        vector<int32_t> model;
    public:
        explicit Solver(const Settings &settings = {}) : context(settings) {}
        void add(const vector<int32_t> &clause) { context.clauses.emplace_back(clause); }
        Status solve() {
            context.start = chrono::steady_clock::now();
            context.vars = 0;
            for (const auto& clause : context.clauses)
                for (int32_t lit : clause)
                    context.vars = max(context.vars, abs(lit));
            Preprocessor pre(context);
            pre.preprocess(context.settings.preprocessing_seconds);
            Internal internal(context);
            Status status = internal.solve();
            if (status != OPTIMAL) exit(1);
            if (status == SATISFIABLE || status == OPTIMAL)
                model = pre.reconstruct(context.model), model.insert(model.begin(), 0);
            return status;
        }

        [[nodiscard]] uint64_t objective_value() const {
            return context.objective_value + context.offset;
        }

        [[nodiscard]] uint64_t bound() const {
            return context.bound + context.offset;
        }

        [[nodiscard]] bool value(const int32_t lit) const {
            assert(1 <= lit && lit <= model.size());
            return model[lit] > 0;
        }
    };
} // namespace Alpaca
