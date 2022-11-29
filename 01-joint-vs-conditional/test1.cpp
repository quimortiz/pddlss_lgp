#include "KOMO/local_utils.hpp"
#include <boost/test/unit_test_suite.hpp>
#define BOOST_TEST_DYN_LINK
#include "Core/util.h"
#include "KOMO/conv_komo_graph.h"
#include "KOMO/komo.h"
#include "LGP/LGP_symconflict.h"
#include "LGP/LGP_tree.h"
#include <boost/config.hpp> // put this first to suppress some VC++ warnings
#include <boost/test/unit_test.hpp>
#include <filesystem>
#include <queue>

std::ofstream out_file;
int num_threads;

std::vector<double> PANDA_Q0 = {0.0, -.5, 0., -2, -0., 2., -.5};

using Match = std::vector<std::pair<int, int>>;

template <typename T> bool in_vector(const T &e, const std::vector<T> &v) {
  return std::find(v.begin(), v.end(), e) != v.end();
}

std::vector<std::string> string_to_tokens(const std::string &in,
                                          char sep = '\n') {
  std::stringstream test(in);
  std::string segment;
  std::vector<std::string> seglist;

  while (std::getline(test, segment, sep)) {
    seglist.push_back(segment);
  }
  return seglist;
}

std::string int_to_string(int number, int len = 4) {
  std::stringstream ss;
  ss << std::setw(len) << std::setfill('0') << number;
  return ss.str();
}

void set_a_var(BGraph &graph, const std::string &name, int time,
               const std::vector<double> &x) {
  // continue here!!

  auto [beg_ref, end_ref] = boost::vertices(graph);
  bool found = false;
  for (auto it = beg_ref; it != end_ref && !found; it++) {
    auto &v = graph[*it];
    if (__startsWith(v.name, name) && v.time == time) {
      v.x = x;
      found = true;
    }
  }
  if (!found) {
    throw std::runtime_error("warning: could not propagate");
    //
  }
}

void propagate_fixed_vars(BGraph &graph) {

  auto [beg_ref, end_ref] = boost::vertices(graph);

  size_t max_time =
      *std::max_element(beg_ref, end_ref, [&](const auto &a, const auto &b) {
        return graph[a].time < graph[b].time;
      });
  size_t start_time = 0;

  for (size_t t = start_time; t <= max_time; t++) {

    for (auto it = beg_ref; it != end_ref; it++) {
      auto &v = graph[*it];
      if (!v.is_var)
        continue;
      if (v.time != t)
        continue;
      if (v.x.size())
        continue;
      auto [nb, ne] = boost::adjacent_vertices(*it, graph);
      for (auto nit = nb; nit != ne; nit++) {
        auto &n = graph[*nit];
        if (__startsWith(n.name, "F_qZeroVel")) {
          auto [nbx, nex] = boost::adjacent_vertices(*nit, graph);

          Qassert(std::distance(nbx, nex) == 2);
          auto other = nbx;
          if (*other == *it) {
            other++;
          }
          // i can only propagate forward
          //
          if (graph[*other].x.size() && graph[*other].time < graph[*it].time) {
            std::cout << "I can propagate, lets copy the value" << std::endl;
            std::cout << "source is: " << graph[*other] << std::endl;
            std::cout << "target is: " << graph[*it] << std::endl;
            graph[*it].x = graph[*other].x;
          }
        }
      }
    }
  }
}

void warmstart_nlp_from_graph(BGraph &graph,
                              shared_ptr<Conv_KOMO_FineStructuredProblem> &mp,
                              bool fix_vars = false) {

  std::cout << "doing graph based warmstart" << std::endl;
  auto [b, e] = boost::vertices(graph);
  for (auto it = b; it != e; it++) {
    auto &v = graph[*it];
    if (v.is_var) {
      if (v.x.size()) {
        auto &vi = mp->__variableIndex(v.id);
        vi.setValue(arr(v.x, false));
      }
    }
  }

  if (fix_vars) {
    uintA conditional_vars = {};
    auto [b, e] = boost::vertices(graph);
    for (auto it = b; it != e; it++) {
      auto &v = graph[*it];
      if (v.is_var) {
        if (v.x.size()) {
          conditional_vars.append(v.id);
        }
      }
    }
    mp->subSelect({}, conditional_vars);
  }
}

void subgraph_with_distance(const BGraph &graph_small, const BGraph &graph_big,
                            std::vector<Match> &matches, bool &match_vars,
                            bool &match_fix_vars, bool verbose = false) {
  double tol = 1e-6;
  match_fix_vars = false;
  Vf2_print_store_callback<BGraph, BGraph> callback(graph_small, graph_big,
                                                    matches);
  match_vars = vf2_subgraph_iso(
      graph_small, graph_big, callback, get(boost::vertex_index, graph_small),
      get(boost::vertex_index, graph_big),
      boost::vertex_order_by_mult(graph_small), boost::always_equivalent(),
      Equivalentq<BGraph>(graph_small, graph_big));

  std::cout << "match vars is " << match_vars << std::endl;
  if (match_vars) {
    auto &match = matches.front();
    double distance = 0.;
    for (auto &[a, b] : match) {
      if (graph_small[a].is_var) {
        if (verbose) {
          std::cout << "***\n";
          std::cout << "***\n";
          std::cout << graph_small[a] << std::endl;
          std::cout << "--\n";
          std::cout << graph_big[b] << std::endl;
          std::cout << "***\n";
          std::cout << "***\n";
        }

        if (graph_big[b].x.size()) {
          distance += euclideanDistance(arr(graph_big[b].x, false),
                                        arr(graph_small[a].x, false));
        }
      }
    }
    std::cout << "distance between fixed variables" << std::endl;
    std::cout << distance << std::endl;
    if (distance < tol) {
      match_fix_vars = true;
      std::cout << "match fix vars is " << match_fix_vars << std::endl;
    }
  }
}

auto extract_interesting_subgraphs(BGraph &graph) {
  const bool write_intermeidate_to_file = false;
  std::vector<BGraph> out;
  BGraph clean;
  boost::copy_graph(make_filtered_graph(
                        graph, boost::keep_all(),
                        Keep_FF<BGraph>(
                            &graph, [](auto &s) { return true; },
                            [](auto &s) {
                              std::cout << s.name << std::endl;
                              std::vector<std::string> del = {
                                  "F_qQuaternionNorms",
                                  "F_PositionRel/0-l_gripper-l_panda_base",
                                  "F_PositionRel/0-r_gripper-r_panda_base",
                                  "F_Pose/1-block", "F_Pose/2-block"};
                              return std::find_if(
                                         del.begin(), del.end(), [&](auto &j) {
                                           return _startsWith(s.name, j);
                                         }) == del.end();
                            },
                            false)),
                    clean);

  if (write_intermeidate_to_file) {
    std::cout << "writing clean " << std::endl;
    std::ofstream file("tmp/g.dot");
    write_graphviz_easy(file, clean, make_label_writer3(clean));
    file.close();
    system_s("dot -Tpdf tmp/g.dot -o tmp/g.pdf");
    // system_s("zathura tmp/g.pdf & ");
    std::cin.get();
  }

  std::vector<std::vector<std::string>> name_refs = {
      {"F_PoseRel/1-block1-block1_ref"},
      {"F_PoseRel/1-block2-block2_ref"},

      {"F_PoseRel/1-block1-r_gripper"},
      {"F_PoseRel/1-block2-r_gripper"},

      {"F_PoseRel/1-block1-l_gripper"},
      {"F_PoseRel/1-block2-l_gripper"},

      {"F_PoseRel/1-block1-block1_ref", "F_PoseRel/1-block1-r_gripper"},
      {"F_PoseRel/1-block2-block2_ref", "F_PoseRel/1-block2-r_gripper"},

      {"F_PoseRel/1-block1-block1_ref", "F_PoseRel/1-block1-l_gripper"},
      {"F_PoseRel/1-block2-block2_ref", "F_PoseRel/1-block2-l_gripper"},

  };

  std::vector<std::vector<int>> vars_each_graph;

  for (auto &name_ref : name_refs) {
    auto &g = clean;
    auto [beg, end] = boost::vertices(g);
    std::vector<int> idxs_id;
    std::vector<std::pair<int, std::vector<int>>> data;

    for (auto it = beg; it != end; it++) {
      auto &v = g[*it];
      if (!v.is_var && in_vector(v.name, name_ref)) {
        std::cout << "found constraint" << std::endl;
        idxs_id.push_back(v.id);
        std::vector<int> neighs;

        auto [nb, ne] = boost::adjacent_vertices(*it, g);
        for (auto nit = nb; nit != ne; nit++) {
          auto &ve = g[*nit];
          neighs.push_back(ve.id);
        }
        data.push_back({v.id, neighs});
      }
    }
    printContainer(idxs_id, ",", true);
    for (auto &d : data) {
      std::cout << d.first << ": ";
      printContainer(d.second);
      std::cout << std::endl;
    }

    BGraph small;
    boost::copy_graph(
        make_filtered_graph(
            graph, boost::keep_all(),
            Keep_FF<BGraph>(
                &graph, [](auto &s) { return true; },
                [&](auto &s) { return in_vector(s.name, name_ref); }, false,
                true)),
        small);

    // get the vars

    std::vector<int> vars;
    std::vector<int> cons;
    vars_and_cons_from_graph(small, vars, cons);

    // check that I don't put the same graph twice

    if (write_intermeidate_to_file) {
      std::cout << "writing small " << std::endl;
      std::ofstream file("tmp/g.dot");
      write_graphviz_easy(file, small, make_label_writer3(small));
      file.close();
      system_s("dot -Tpdf tmp/g.dot -o tmp/g.pdf");
      std::cin.get();
    }

    BGraph small2;
    boost::copy_graph(
        make_filtered_graph(clean, boost::keep_all(),
                            Keep_FF<BGraph>(
                                &clean,
                                [&](auto &s) {
                                  return std::find(vars.begin(), vars.end(),
                                                   s.id) != vars.end();
                                },
                                [&](auto &s) { return true; }, true, false)),
        small2);

    if (boost::num_vertices(small2)) {

      bool show_added = false;
      if (show_added) {
        std::cout << "writing found" << std::endl;
        std::ofstream file("tmp/g.dot");
        write_graphviz_easy(file, small2, make_label_writer3(small2));
        file.close();
        system_s("dot -Tpdf tmp/g.dot -o tmp/g.pdf");
        std::cin.get();
      }

      std::vector<int> varsx;
      std::vector<int> consx;
      vars_and_cons_from_graph(small2, varsx, consx);

      bool is_repeated = false;
      for (auto &vv : vars_each_graph) {
        if (vv == varsx) {
          std::cout << "this graph is the same as a previous one" << std::endl;
          is_repeated = true;
          break;
        }
      }

      if (is_repeated)
        continue;

      vars_each_graph.push_back(varsx);

      out.push_back(small2);
    }

    // auto &graph_small = small2;
    // auto &graph_big = graph;
    // Qvf2_print_callback<BGraph, BGraph> callback(graph_small, graph_big);

    // std::cout << "checking if it is a subgraph " << std::endl;
    // bool match = vf2_subgraph_iso(
    //     graph_small, graph_big, callback, get(boost::vertex_index,
    //     graph_small), get(boost::vertex_index, graph_big),
    //     boost::vertex_order_by_mult(graph_small),
    //     boost::always_equivalent(), Equivalentq<BGraph>(graph_small,
    //     graph_big));

    // BOOST_ASSERT(match);
  }
  return out;
};

std::vector<std::string> conflict_files;

int global_node_id = 0;
int global_conflict_counter = 0;
int global_id_interesting_subgraphs = 0;

auto get_nlps_and_graph(rai::LGP_Node *node) {

  auto bound = rai::BD_seq;
  bool verbose = false;

  StringAA collisions_list = {};
  StringAA collisions_pairwise = {};
  node->prepareProblem(bound, verbose, collisions_list, collisions_pairwise);
  auto factored_nlp = node->problem(bound).komo->mp_FineStructured();
  auto new_variableIndex =
      break_chain(factored_nlp->__variableIndex, separate_by_name,
                  [](auto &s) { return true; });

  factored_nlp->__variableIndex = new_variableIndex;
  factored_nlp->recompute_frame2ids_and_features();
  factored_nlp->komo.run_prepare(0);
  factored_nlp->report(cout, 3);
  arr default_x = factored_nlp->komo.x.copy();
  BGraph graph = create_boost_graph(*factored_nlp, false);
  return std::make_pair(graph, factored_nlp);
}

void fix_subgraph(BGraph &graph, BGraph &subgraph) {

  auto [beg, end] = boost::vertices(subgraph);
  auto [beg_ref, end_ref] = boost::vertices(graph);

  for (auto it = beg; it != end; it++) {
    auto &v = subgraph[*it];
    if (!v.is_var)
      continue;
    auto it_in_ref = std::find_if(beg_ref, end_ref, [&](const auto &s) {
      return graph[s].is_var && graph[s].id == v.id;
    });
    Qassert(it_in_ref != end_ref);

    auto &x = v.x;
    // arr x = mp->__variableIndex(v.id).getValue(true);
    // arr x2 = mp->__variableIndex(v.id).getValue(true);
    graph[*it_in_ref].x = std::vector(x.begin(), x.end());
  }
}

std::string from_sas_to_komo(std::string s) {

  // sas: Atom on(block2_ref, block2)
  // komo: (on block2_ref block2)

  std::replace(s.begin(), s.end(), '(', ' ');
  std::replace(s.begin(), s.end(), ')', ' ');
  std::replace(s.begin(), s.end(), ',', ' ');

  std::stringstream test(s);
  std::string segment;
  std::vector<std::string> seglist;
  std::cout << "s " << s << std::endl;

  while (std::getline(test, segment, ' ')) {
    if (segment.size())
      seglist.push_back(segment);
  }

  Qassert(seglist.size());
  Qassert(seglist.at(0) == "Atom");

  for (auto &s : seglist) {
    std::cout << s << std::endl;
  }
  if (seglist.size() == 4) {
    std::string komo =
        "(" + seglist.at(1) + " " + seglist.at(2) + " " + seglist.at(3) + ")";
    return komo;

  } else {
    throw std::runtime_error("not implemented");
  }
}

bool is_subset(const std::vector<std::string> &ref,
               const std::vector<std::string> &small) {

  return all_of(small.begin(), small.end(), [&](auto &s) {
    return find_if(ref.begin(), ref.end(), [&](auto &r) { return s == r; }) !=
           ref.end();
  });
}

bool is_subsequence(const std::vector<std::vector<std::string>> &ref,
                    const std::vector<std::vector<std::string>> &small) {

  for (size_t i = 0; i < ref.size() - small.size() + 1; i++) {
    bool subseq = true;
    for (size_t j = 0; j < small.size(); j++) {
      if (!is_subset(ref.at(i + j), small.at(j))) {
        subseq = false;
        break;
      }
    }
    if (subseq)
      return true;
  }
}

std::vector<std::vector<std::string>> get_state_sequence(rai::LGP_Node *node) {

  if (!node) {
    return {};
  }

  auto node_path = node->getTreePath();

  std::vector<std::vector<std::string>> state_sequence;
  for (auto &n : node_path) {
    std::vector<std::string> state;
    for (uint i = 0; i < n->folState->N; i++) {
      std::stringstream ss;
      {
        if (n->folState->elem(i))
          n->folState->elem(i)->write(ss);
        state.push_back(ss.str());
      }
    }
    state_sequence.push_back(state);
  }
  return state_sequence;
}

std::vector<std::vector<std::string>>
pddl_conflict(const std::string &filename) {

  std::cout << filename << std::endl;
  std::ifstream file(filename);
  Qassert(file.is_open());

  std::string line;
  std::vector<std::vector<std::string>> block;
  std::vector<std::string> lines;

  std::cout << "hello" << std::endl;
  while (std::getline(file, line)) {
    std::cout << "line " << line << "*" << std::endl;
    if (line == "") {
      if (lines.size()) {
        block.push_back(lines);
        lines.clear();
      }
    } else {
      lines.push_back(line);
    }
  }
  if (lines.size())
    block.push_back(lines);

  return block;
}

auto komo_state_sequence_from_PDDL(const std::string &filename) {

  auto out = pddl_conflict(filename);

  for (auto &v : out) {
    std::cout << "***" << std::endl;
    for (auto &vv : v) {
      std::cout << vv << std::endl;
    }
  }

  std::vector<std::vector<std::string>> komo_style;

  for (auto &v : out) {
    std::vector<std::string> neww;
    for (auto &vv : v) {
      neww.push_back(from_sas_to_komo(vv));
    }
    komo_style.push_back(neww);
  }

  std::cout << "translation " << std::endl;
  for (auto &v : komo_style) {
    std::cout << "***" << std::endl;
    for (auto &vv : v) {
      std::cout << vv << std::endl;
    }
  }
  return komo_style;
}

bool replace(std::string &str, const std::string &from, const std::string &to) {
  size_t start_pos = str.find(from);
  if (start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;
}

bool files_in_folder(const std::string &folder) {
  for (const auto &entry : std::filesystem::directory_iterator(folder)) {
    return true;
  }
  return false;
}

BOOST_AUTO_TEST_CASE(hello) { std::cout << "hello world " << std::endl; }

BOOST_AUTO_TEST_CASE(store_graph) {
  // BOOST_CHECK(false);
}

BOOST_AUTO_TEST_CASE(easy_minimal_conflict) {}

BOOST_AUTO_TEST_CASE(plan_solve_join) {

  auto argv = boost::unit_test::framework::master_test_suite().argv;
  auto argc = boost::unit_test::framework::master_test_suite().argc;

  rai::initCmdLine(argc, argv);

  int seed = rai::getParameter<int>("seed", 0);
  rnd.seed(seed);

  rai::String folFile =
      rai::getParameter<rai::String>("folFile", STRING("none"));
  rai::String confFile =
      rai::getParameter<rai::String>("confFile", STRING("none"));
  rai::String plan = rai::getParameter<rai::String>("plan", STRING("none"));

  bool visualize = rai::getParameter<bool>("vis", false);

  rai::Configuration C;
  C.addFile(confFile);
  if (visualize)
    C.watch(true);

  rai::LGP_Tree lgp(C, folFile);

  auto node = lgp.walkToNode(plan);

  rai::BoundType bound = rai::BoundType::BD_seq;
  bool collisions = true;
  int verbose = 1;
  const rai::String folder = "";
  const rai::String &conflict = "";
  bool relaxation = false;
  const StringAA &collisions_list = {};
  const StringAA &collisions_pairwise = {};
  std::map<std::string, double> *data = nullptr;
  OptConflict opt = OptConflict();
  bool extra_try = false;

  bool check_komo = false;
  if (check_komo) {
    node->optBound(bound, collisions, verbose, folder, conflict, relaxation,
                   collisions_list, collisions_pairwise, data, opt, extra_try);
    node->problem(bound).komo->view(true);
    node->problem(bound).komo->view_play(true);
  }

  node->prepareProblem(bound, verbose, collisions_list, collisions_pairwise);

  auto factored_nlp = node->problem(bound).komo->mp_FineStructured();

  // solve now using the factored nlp
  //

  auto new_variableIndex =
      break_chain(factored_nlp->__variableIndex, separate_by_name,
                  [](auto &s) { return true; });

  factored_nlp->__variableIndex = new_variableIndex;
  // std::cout << factored_nlp->__variableIndex << std::endl;
  factored_nlp->recompute_frame2ids_and_features();
  factored_nlp->komo.run_prepare(0);
  factored_nlp->report(cout, 3);

  arr default_x = factored_nlp->komo.x.copy();

  BGraph graph = create_boost_graph(*factored_nlp, false);
  std::cout << "writing graph " << std::endl;
  std::ofstream file("tmp/g.dot");
  write_graphviz_easy(file, graph, make_label_writer3(graph));
  factored_nlp->subSelect({}, {});

  bool solve_graph_full = false;
  if (solve_graph_full) {

    std::cout << std::endl;

    // solve
    auto out = MP_Solver().setProblem(factored_nlp).solve();

    out->write(std::cout);
  }
  bool solve_forward = false;

  if (solve_forward) {
    Solve_nlp_with_counter solver;

    for (size_t i = 0; i < 4; i++) {
      // int i = 0;
      std::cout << "i " << i << std::endl;
      BGraph graph_slice;
      boost::copy_graph(
          make_filtered_graph(graph, boost::keep_all(),
                              Keep_FF<BGraph>(
                                  &graph,
                                  [&](const GBNode &s) {
                                    return s.time >= i && s.time <= i + 1;
                                  },
                                  [&](auto &s) { return true; }, true)),
          graph_slice);

      bool flag = solver.solve_nlp(factored_nlp, default_x, graph_slice,
                                   visualize, true);

      BOOST_CHECK_EQUAL(flag, true);
      factored_nlp->komo.view(true);
      factored_nlp->komo.view_play(true);

      for (auto &v : factored_nlp->__variableIndex) {
        std::cout << "***" << std::endl;
        std::cout << v << std::endl;
        std::cout << "value " << v.getValue() << std::endl;
        std::cout << "***" << std::endl;
      }

      auto [vars, cons] = vars_cons_from_graph(graph_slice);
      for (auto &j : vars) {
        auto &v = graph[j];
        factored_nlp->__variableIndex(v.id).fixed = true;
      }
    }
  }

  bool solve_sequential_backward = true;
  std::cout << "solving backward" << std::endl;
  if (solve_sequential_backward) {
    // backward is a little bit strange...
    // i have to solve for x_{i-1}, x_i, x_{i+1}, where x_{i+1} is fixed.
    // and then only fix x_{i}

    Solve_nlp_with_counter solver;

    int num_actions = 4;
    for (size_t i = 0; i < num_actions + 1; i++) {
      BGraph graph_slice;
      std::cout << "i " << i << std::endl;
      boost::copy_graph(
          make_filtered_graph(graph, boost::keep_all(),
                              Keep_FF<BGraph>(
                                  &graph,
                                  [&](const GBNode &s) {
                                    return s.time >= num_actions - i - 1 &&
                                           s.time <= num_actions - i + 1;
                                  },
                                  [&](auto &s) { return true; }, true)),
          graph_slice);

      auto [vi, ci] = vars_cons_from_graph(graph_slice);

      std::cout << "VARS: " << vi << std::endl;
      std::cout << "CONS: " << ci << std::endl;

      bool flag = solver.solve_nlp(factored_nlp, default_x, graph_slice,
                                   visualize, true);

      BOOST_CHECK_EQUAL(flag, true);
      factored_nlp->komo.view(true);
      factored_nlp->komo.view_play(true);

      for (auto &v : factored_nlp->__variableIndex) {
        std::cout << "***" << std::endl;
        std::cout << v << std::endl;
        std::cout << "value " << v.getValue() << std::endl;
        std::cout << "***" << std::endl;
      }

      auto [vars, cons] = vars_cons_from_graph(graph_slice);
      std::cout << "check what to fix " << std::endl;
      for (auto &j : vars) {
        auto &v = graph[j];
        std::cout << "v " << v << std::endl;
        auto &vv = factored_nlp->__variableIndex(v.id);
        std::cout << "vv " << vv << std::endl;
        // i only have to fix the variables of the TOP.
        if (vv.time == 4 - i) {
          std::cout << "fixing " << v << std::endl;
          vv.fixed = true;
        }
      }
    }
    // TODO: double check that all constraints are fulfiled!
  }
}

enum class TYPE {
  with_values,
  with_variables,
};

bool path_to_node_has_conflict(rai::LGP_Node *node,
                               const std::string &conflict_folder) {
  auto state_sequence = get_state_sequence(node);
  for (const auto &entry :
       std::filesystem::directory_iterator(conflict_folder)) {
    auto komo_style = komo_state_sequence_from_PDDL(entry.path());
    if (is_subsequence(state_sequence, komo_style))
      return true;
  }
  return false;
}

std::pair<int, rai::String>
compute_pddl_heuristic(rai::LGP_Node *node,
                       const std::string &path_to_downward_folder,
                       const std::string &conflict_folder) {

  using cs = const std::string;
  auto &tree = node->tree;

  cs sas = "sas";
  auto plan_c = "plan";
  auto sas_c = "sas";
  auto states_c = "states";

  std::string id = rand_id(6);

  tree.fol.setState(node->folState);
  tree.fol.writePDDLfiles(("tmp/child_" + id).c_str(), true);
  cs plan_file = string_format("tmp/%s.%s.txt", plan_c, id.c_str());
  std::string sas_file = string_format("tmp/%s.%s.txt", sas_c, id.c_str());
  cs states_file = string_format("tmp/%s.%s.txt", states_c, id.c_str());

  cs generate_sas_cmd = "/usr/bin/python3 " + path_to_downward_folder +
                        "/src/translate/translate.py " + "tmp/child_" + id +
                        ".domain.pddl " + "tmp/child_" + id + ".problem.pddl";

  cs modify_sas_cmd = "/usr/bin/python3 " + path_to_downward_folder +
                      "/src/translate/modify_sas.py ";

  cs down_cmd = path_to_downward_folder + "/builds/release/bin/"
                                          "downward --search 'astar(ff())' ";

  system_s(generate_sas_cmd + " --sas-file " + sas_file);

  bool modify_sas_with_conflicts = true;

  if (modify_sas_with_conflicts) {

    std::cout << "check if plan until here contains the action" << std::endl;

    std::cout << "MODIFYING  heuristic" << std::endl;

    int counter = 0;
    for (const auto &entry :
         std::filesystem::directory_iterator(conflict_folder)) {
      counter++;
    }

    if (counter > 0) {
      std::cout << "there are conflicts in the folder" << std::endl;
      cs new_sas_file = sas_file + ".c.txt";
      cs modify_sas_args =
          string_format("%s %s %s", conflict_folder.c_str(), sas_file.c_str(),
                        new_sas_file.c_str());

      cs path_to_downward_folder = "/home/quim/stg/lgp-pddl/downward/";
      cs modify_sas_cmd = "/usr/bin/python3 " + path_to_downward_folder +
                          "/src/translate/modify_sas.py ";
      system_s(modify_sas_cmd + modify_sas_args);
      sas_file = new_sas_file;
    }
  }

  cs down_args =
      string_format("--internal-plan-file %s < "
                    "%s --internal-state-file  %s",
                    plan_file.c_str(), sas_file.c_str(), states_file.c_str());

  system_s(down_cmd + down_args);

  if (!std::filesystem::exists(plan_file)) {
    std::cout << "Warning: "
              << "i can not solve planning task" << std::endl;
    std::cout << "Assign BIG_NUMBER_COST_TO_GO as heuristic" << std::endl;
    return {false, ""};
  }

  rai::String plan = rai::String(FILE(plan_file.c_str()));

  // cut the last line comment with ';'
  uint i = plan.N;
  for (; i--;)
    if (plan(i) == ';')
      break;
  plan.resize(i, true);

  cout << "FOUND PLAN: \n" << plan << endl;

  cout << "plan until here " << std::endl;
  return {true, plan};
}

void find_conflict(rai::LGP_Node *node,
                   shared_ptr<Conv_KOMO_FineStructuredProblem> factored_nlp,
                   const arr &default_x, const std::string &conflict_folder) {

  using cs = const std::string;
  cs komo_conflict = "komo_conflict";
  cs conflict_sas = "conflict_sas";

  Solve_nlp_with_counter solver;
  solver.init_database_with_init_guess(*factored_nlp, default_x);
  solver.use_local_skip = true;
  solver.use_global_skip = true;
  solver.max_it_solve = 1;

  OptConflict opt = OptConflict();
  int num_conflicts = 1;
  bool visualize = false;
  factored_nlp->reset_subselect(default_x);
  auto graph = create_boost_graph(*factored_nlp, false);

  std::vector<BGraph> conflicts = conflict_extraction(
      graph, factored_nlp, solver, default_x, visualize, opt, num_conflicts);

  Qassert(conflicts.size());
  auto &cconflict = conflicts.front();
  uintA vars_sorted, cons_sorted;
  boost::tie(vars_sorted, cons_sorted) = vars_cons_from_graph(cconflict);

  std::cout << "found infeasible relaxation" << std::endl;
  Relaxation infeas_relaxation = {0, 0, {false, vars_sorted, cons_sorted, {}}};

  auto &r = infeas_relaxation;
  CHECK(r.r.cons.N && r.r.vars.N, "infeas can not be empty");

  std::vector<int> vars(r.r.vars.p, r.r.vars.p + r.r.vars.N);
  std::vector<int> cons(r.r.cons.p, r.r.cons.p + r.r.cons.N);

  std::cout << "conflict vars are: " << std::endl;

  for (auto &v : vars)
    std::cout << factored_nlp->__variableIndex(v) << std::endl;

  int offset = 0;
  auto sym_states_conflict =
      get_sym_conflict(vars, node->getTreePath(), *factored_nlp, offset);

  print_states_conflict(std::cout, sym_states_conflict, true);

  std::string path_to_utils = "../rai/utils";
  using cs = const std::string;
  std::ofstream file_conflicts(
      komo_conflict + "." + std::to_string(global_conflict_counter) + ".txt");
  print_states_conflict(file_conflicts, sym_states_conflict, true);
  file_conflicts.close();
  system_s(
      string_format("bash " + path_to_utils + "/cmd_conflict.sh %s %d %s %s",
                    komo_conflict.c_str(), global_conflict_counter,
                    conflict_sas.c_str(), conflict_folder.c_str()));
  global_conflict_counter++;
}

struct Compu_node {
  // lets try to order: minimize the time to reach the goal.
  // with values: heuristic, never added again
  // with variables: heuristic + expands, added again. I only do the
  // symbolic expansion once.
  // root is of type: with values

  Compu_node() { id = ++global_node_id; }

  const int BIG_NUMBER_COST_TO_GO = 1000;
  const double FEASIBLE_THRESHOLD = .1;
  bool feasible = true;
  int id;
  const std::string conflict_folder = "__conflicts";
  rai::LGP_Node *node; // a Compu Node is a LGP Node with some
                       // assigned values!
  rai::Configuration *C;
  rai::LGP_Tree *tree; // or fol?
  rai::String plan;
  std::vector<std::string> logic_states; // other representation?
  BGraph graph;
  std::map<rai::String, arr> values;
  TYPE type = TYPE::with_values;
  Compu_node *parent;
  std::vector<Compu_node *> childs;
  bool symbolic_expansion_done = false; // it can only be done once
  int depth =
      0; // takes into account the num of expands operations in the parent
  int tree_depth = 0; // does not consider this
  int heuristic = 0;
  int path_from_root = 0;
  bool terminal = false;
  int max_time_fixed = 0;
  std::vector<int> matches_feasible_subgraphs = {};

  std::string path_feasible_samples = "tmp_feasible_graphs";

  StringAA collisions_list = {};
  StringAA collisions_pairwise = {};

  using cs = const std::string;
  cs path_to_downward_folder = "/home/quim/stg/lgp-pddl/downward/";

  std::vector<rai::String> tabu_plans;

  bool highlight = false; // you can use this to highlight nodes in
                          // graphviz

  std::vector<int> expands; // id's of expands for debugging, set from
                            // outside
  int num_expands = 0;

  void compute_heuristic() {

    const bool invalidate_childs = true;
    const bool check_if_path_to_node_is_conflict = true;

    Qassert(tree);
    Qassert(node);

    // the parent can set the heuristic of this to
    // BIG_NUMBER_COST_TO_GO
    if (heuristic == BIG_NUMBER_COST_TO_GO)
      return;

    if (check_if_path_to_node_is_conflict &&
        path_to_node_has_conflict(node, conflict_folder)) {

      std::cout << "it contains an infeasible subsequence!" << std::endl;
      heuristic = BIG_NUMBER_COST_TO_GO;
      feasible = false;

      if (invalidate_childs) {
        std::cout << "invalidating all childs" << std::endl;
        std::queue<Compu_node *> queue;
        queue.push(this);
        while (queue.size()) {
          auto n = queue.front();
          n->feasible = false;
          n->heuristic = BIG_NUMBER_COST_TO_GO;
          for (auto &nn : n->childs)
            queue.push(nn);
          queue.pop();
        }
      }
      return;
    }

    auto [flag, plan] =
        compute_pddl_heuristic(node, path_to_downward_folder, conflict_folder);

    if (flag)
      heuristic = std::count(plan.p, plan.p + plan.N, '(');
    else
      heuristic = BIG_NUMBER_COST_TO_GO;

    // block plans that are in the the tabu list
    rai::String plan_until_here = node->getTreePathString('\n');
    rai::String plan_from_root = plan_until_here + plan;
    std::cout << "plan from root" << plan_from_root << std::endl;
    std::cout << "blocking plan if it is in the tabu list" << std::endl;
    for (auto &p : tabu_plans) {
      if (p == plan_from_root) {
        // std::cout << p << std::endl;
        // std::cout << plan_from_root << std::endl;
        // std::cout << plan_until_here << std::endl;
        std::cout << "plan is in tabu list!" << std::endl;
        std::cout << "plan_from_root:" << plan_from_root << std::endl;
        heuristic = BIG_NUMBER_COST_TO_GO;
        break;
      }
    }
  }

  std::vector<Compu_node *> expand(bool *recompute_heuristic = nullptr) {
    const bool check_conflicts = true;
    std::vector<Compu_node *> new_nodes;
    num_expands++;
    bool expand_symbolic =
        true; // this will be set to false if optimization fails
    if (num_expands > 1 && type == TYPE::with_values) {
      throw std::runtime_error("node with values can only be expanded once");
    }
    if (type == TYPE::with_variables) {

      auto [graph, factored_nlp] = get_nlps_and_graph(node);
      factored_nlp->report(cout, 3);

      for (auto &v : factored_nlp->__variableIndex)
        std::cout << v << std::endl;

      arr default_x = factored_nlp->komo.x.copy();

      // example: id 0 dim 7 dofs [block1.202 ]  name block1 nameID block1.202
      // time 0
      std::cout << "variables of the problems are " << std::endl;
      for (auto &v : factored_nlp->__variableIndex)
        std::cout << v << std::endl;

      auto &vi = factored_nlp->__variableIndex;

      auto [beg, end] = boost::vertices(graph);

      uintA condition_variable = {};
      for (auto &[k, v] : values) {
        std::cout << "fixed vars are " << std::endl;
        std::cout << "[k,v] " << k << " " << v << std::endl;
        auto it = std::find_if(vi.begin(), vi.end(),
                               [&](auto &s) { return s.nameID == k; });
        Qassert(it != vi.end());
        it->setValue(v);
        condition_variable.append(it->id);

        auto it_g = std::find_if(beg, end, [&](const auto &s) {
          return graph[s].is_var && graph[s].name == std::string(k.p);
        });
        Qassert(it_g != end);
        graph[*it_g].x = std::vector(v.begin(), v.end());
      }

      auto it = std::max_element(beg, end, [&](const auto &a, const auto &b) {
        return graph[a].time < graph[b].time;
      });

      int max_time = graph[*it].time;

      auto it2 = std::max_element(beg, end, [&](const auto &a, const auto &b) {
        return bool(graph[a].x.size()) * graph[a].time <
               bool(graph[b].x.size()) * graph[b].time;
      });

      Qassert(graph[*it2].x.size());
      int max_time_conditional = graph[*it2].time;

      std::cout << "max time all and conditional " << max_time << " "
                << max_time_conditional << std::endl;

      bool use_database_samples = true;

      if (use_database_samples) {
        // auto &vii = factored_nlp->__variableIndex;
        //

        BGraph graph_slice;
        boost::copy_graph(
            make_filtered_graph(graph, boost::keep_all(),
                                Keep_FF<BGraph>(
                                    &graph,
                                    [&](const GBNode &s) {
                                      return s.time >= max_time_conditional &&
                                             s.time <= max_time;
                                    },
                                    [&](auto &s) { return true; }, true)),
            graph_slice);

        std::vector<GraphStruct> graphs =
            load_graphs(path_feasible_samples.c_str());

        // we have to take only the graph of: [ last fixed, max_time ]

        // i want the same max-min interval
        int it_counter = 0;

        std::vector<std::vector<std::pair<int, int>>> matches;
        BGraph chosen;
        bool match_vars, match_fix_vars;

        int match_index = 0;
        for (auto &g_candidate : graphs) {

          if (in_vector(match_index, matches_feasible_subgraphs)) {
            std::cout << "this sugraphs has already been matched for this "
                         "expansion -- continueing"
                      << std::endl;
            continue;
          }

          std::cout << it_counter << std::endl;
          auto &graph_small = g_candidate.g;
          auto &graph_big = graph_slice;

          auto [gbeg, gend] = boost::vertices(graph_small);
          auto [itmin, itmax] = std::minmax_element(
              gbeg, gend, [&](const auto &a, const auto &b) {
                return graph_small[a].time < graph_small[b].time;
              });

          int min_time_candidate = graph_small[*itmin].time;
          int max_time_candidate = graph_small[*itmax].time;
          int dif = max_time_candidate - min_time_candidate;
          if (dif != max_time - max_time_conditional) {
            std::cout << "skipping, different times" << std::endl;
            continue;
          }
          subgraph_with_distance(graph_small, graph_big, matches, match_vars,
                                 match_fix_vars);
          std::cout << "result is " << std::endl;
          std::cout << match_vars << " " << match_fix_vars << std::endl;

          if (match_vars && match_fix_vars) {
            chosen = graph_small;
            std::cout << "info: we break at first match" << std::endl;
            matches_feasible_subgraphs.push_back(match_index);
            break;
          } else {
            // I clear the vector because maybe there was a symbolic
            // match only.
            matches.clear();
          }
          match_index++;
        }

        //  Now lets do the warmstart

        // Qassert(matches.size());
        if (matches.size()) {
          std::cout << "there is a match" << std::endl;
          for (auto &[a, b] : matches.front()) {
            if (chosen[a].is_var && !graph_slice[b].x.size())
              graph_slice[b].x = chosen[a].x;
          }

          // I have to copy from graph_slice back to graph
          fix_subgraph(graph, graph_slice);
        }

        // TODO: think about what to do with the other robot.
        // propagate_a_var(graph_big, "l_panda", 2);

        propagate_fixed_vars(graph);

        bool fixed_non_used_robot = true;

        if (fixed_non_used_robot) {

          std::cout << "fixing only last decision" << std::endl;
          std::stringstream ss;
          ss << *node->decision;
          std::string ss_str = ss.str();

          std::vector<std::string> actions = string_to_tokens(plan.p);

          Qassert(actions.size() == max_time);

          for (size_t t = max_time_conditional + 1; t <= max_time; t++) {

            std::stringstream ss(actions.at(t));
            auto ss_str = ss.str();
            std::cout << "checking action at time " << ss_str << " " << t
                      << std::endl;

            if (ss_str.find("l_") == std::string::npos) {
              std::cout << "fixing l_ at " << t << std::endl;
              set_a_var(graph, "l_", t, PANDA_Q0);
            }
            if (ss_str.find("r_") == std::string::npos) {
              std::cout << "fixing r_ at " << t << std::endl;
              set_a_var(graph, "r_", t, PANDA_Q0);
            }
          }
        }
      }

      bool fix_vars = true;
      warmstart_nlp_from_graph(graph, factored_nlp, fix_vars);

      shared_ptr<SolverReturn> out;
      std::cout << " before optimization " << std::endl;
      factored_nlp->report(cout, 5);
      out = MP_Solver().setProblem(factored_nlp).solve();
      out->write(std::cout);
      std::cout << std::endl;
      std::cout << " after optimization " << std::endl;
      factored_nlp->report(cout, 5);

      bool visualize_komo = true;
      if (visualize_komo) {
        factored_nlp->komo.view_play(true);
      }

      for (auto &v : factored_nlp->__variableIndex)
        std::cout << v << std::endl;

      bool assignment_feasible =
          fabs(out->ineq) + fabs(out->eq) < FEASIBLE_THRESHOLD;
      if (!assignment_feasible) {
        expand_symbolic = false;
        std::cout << "CAUTION, optimizaton is infeasible" << std::endl;
        std::cout << "setting expand symbolic to false" << std::endl;
        if (check_conflicts && max_time_conditional == 0) {
          std::cout << "Warning: check conflicts is only working  for no "
                       "condition variables"
                    << std::endl;
          find_conflict(node, factored_nlp, default_x, conflict_folder);
          if (recompute_heuristic) {
            *recompute_heuristic = true;
          }
        }

      } else {

        // I fix all the graph
        get_values_from_nlp(*factored_nlp, graph);

        auto interesting_subgraphs = extract_interesting_subgraphs(graph);

        // but I should only consider the newly assigned variables!
        // save to folder

        std::string command2 = "mkdir -p " + path_feasible_samples;
        system_s(command2.c_str());

        std::cout << "writing intersesting subgraphs" << std::endl;
        for (auto &g : interesting_subgraphs) {
          std::cout << global_id_interesting_subgraphs << std::endl;
          ofstream file_graphviz(
              path_feasible_samples + "/" +
              int_to_string(global_id_interesting_subgraphs) + ".dot");
          write_graphviz_easy(file_graphviz, g, make_label_writer3(g));
          ofstream file_graphboost(
              path_feasible_samples + "/" +
              int_to_string(global_id_interesting_subgraphs) + ".dat");
          file_graphboost << boost::write(g);
          std::cin.get();
          global_id_interesting_subgraphs++;
        }

        // extract samples and add to the database
      }

      // I have to fix the variables for the child!
      std::map<rai::String, arr> new_values;
      for (auto &v : factored_nlp->__variableIndex) {
        std::cout << v << std::endl;
        new_values[v.nameID] = v.getValue(true);
      }

      std::cout << "fixed vars are" << std::endl;
      for (auto &[k, v] : new_values)
        std::cout << "k,v: " << k << " " << v << std::endl;

      Compu_node *cc = new Compu_node();
      cc->node = node; // it points to the same node
      cc->tree_depth = tree_depth;
      cc->feasible = assignment_feasible;
      cc->max_time_fixed = max_time;
      cc->depth = depth + num_expands - 1; //
                                           // first time it will be = depth,
                                           // then depth+1,...
      cc->tree = tree;
      cc->heuristic = heuristic;    // it has the same heuristic
      cc->type = TYPE::with_values; // it only has values
      cc->values = new_values;
      cc->terminal = terminal;
      childs.push_back(cc);
      new_nodes.push_back(cc); // is it a goal?
    }

    if (!symbolic_expansion_done && expand_symbolic && feasible) {
      Qassert(node);
      node->expand(); // if it is goal, it does not have children
      for (rai::LGP_Node *ch : node->children) {
        // I sould create new nodes
        Compu_node *cc = new Compu_node();
        cc->node = ch;
        cc->terminal = ch->isTerminal;
        cc->depth = depth + 1;
        cc->tree_depth = tree_depth + 1;
        cc->max_time_fixed = max_time_fixed;
        cc->tree = tree;
        cc->compute_heuristic();
        cc->type = TYPE::with_variables;
        cc->values = values;
        childs.push_back(cc);
        new_nodes.push_back(cc);
      };
      symbolic_expansion_done = true;
    }
    return new_nodes;
  };

  ~Compu_node() {
    for (auto &c : childs) {
      delete c;
    }
  }

  void write(std::ostream &os) {

    os << "\n***NODE***" << std::endl;
    os << "symbolic_expansion_done:" << symbolic_expansion_done << std::endl;
    os << "tree depth" << tree_depth << std::endl;
    os << "depth:" << depth << std::endl;
    os << "heuristic:" << heuristic << std::endl;
    os << "path_from_root:" << path_from_root << std::endl;
    os << "id:" << id << std::endl;
    os << "type:" << int(type) << std::endl;
    os << "node in LGP:" << std::endl;
    node->write(os);
    os << "\nfixed values " << std::endl;
    for (auto &[k, v] : values) {
      os << " k,v: " << k << " " << v << std::endl;
    }

    os << "***NODE***\n" << std::endl;
  }

  std::string to_graphviz() {
    std::stringstream ss;
    write(ss);
    std::string shape;
    std::string label;
    if (type == TYPE::with_values)
      shape = "box";
    else {
      shape = "box";
      {
        std::stringstream ss;
        ss << *node->decision;
        std::string ss_str = ss.str();
        replace(ss_str, "pick", "I");
        replace(ss_str, "place", "A");
        replace(ss_str, "r_gripper", "R");
        replace(ss_str, "l_gripper", "L");
        replace(ss_str, "block", "B");
        replace(ss_str, "block", "B"); // again, to replace second
                                       // occurrence xd
        replace(ss_str, "goal", "g");
        replace(ss_str, "table", "T");
        label = "a:" + ss_str;
      }
    }

    label += " h:" + std::to_string(heuristic);
    label += ",d:" + std::to_string(depth);
    label += ",td:" + std::to_string(tree_depth);
    label += ",ne:" + std::to_string(num_expands);
    label += ",Tf:" + std::to_string(max_time_fixed);
    label += ",exp:";
    std::string vec = "";
    for (auto &e : expands)
      vec += std::to_string(e) + " ";
    label += vec;

    std::string color = type == TYPE::with_values ? "lightgray" : "white";

    if (terminal)
      color = type == TYPE::with_values ? "lightgreen" : "lightblue";

    if (highlight)
      color = "yellow";

    if (!feasible)
      color = "red";

    std::string out = "[shape=" + shape + " style=filled fillcolor=" + color +
                      " label=\"" + label + "\" ]";
    return out;
  };
};

void tree_to_graphviz(std::ostream &os, Compu_node *root) {

  if (!root) {
    os << "digraph D {\n";
    os << "}\n";
    return;
  }

  os << "digraph D {\n";

  Compu_node *node = root;
  int global_id = 0;

  std::map<int, std::vector<int>> childs_map;

  std::queue<std::pair<int, Compu_node *>> queue;
  queue.push(std::make_pair(global_id, node));

  // add nodes
  while (queue.size()) {
    auto [id, node] = queue.front();
    os << id << " " << node->to_graphviz() << std::endl;
    queue.pop();
    std::vector<int> childs;
    for (auto &c : node->childs) {
      global_id++;
      queue.push(std::make_pair(global_id, c));
      childs.push_back(global_id);
    }
    childs_map.insert({id, childs});
  }
  std::cout << "done " << std::endl;

  // add edges
  std::string arrow = " -> ";
  for (auto &[k, v] : childs_map) {
    if (v.size()) {
      for (auto &c : v) {
        os << k << arrow << c << std::endl;
      }
    }
  }

  os << "}\n";
}

struct Compu_tree {

  Compu_node *root;
  int max_it = 1000;
  bool interactive = false;
  bool open_viewer_first_time = true;
  bool block_found_plans = false;

  std::vector<Compu_node *> open_nodes; // queue for search
  std::vector<Compu_node *> goals;
  Compu_node *pop_and_erase_best() {

    std::cout << "printing nodes " << std::endl;
    for (auto &n : open_nodes) {
      n->write(std::cout);
    }

    auto score_with_expands = [](auto &a) {
      return a->heuristic + a->depth + a->num_expands;
    };

    // heuristic distance to goal + cost to arrive
    auto fun = [](auto &a, auto &b) {
      return a->heuristic + a->depth < b->heuristic + b->depth;
    };

    // heuristic distance to goal + cost to arrive + expands
    auto fun2 = [&](auto &a, auto &b) {
      return score_with_expands(a) < score_with_expands(b);
    };

    // we minimize only the cost to generate one sample.
    auto fun3 = [&](auto &a, auto &b) {
      return score_with_expands(a) - a->max_time_fixed <
             score_with_expands(b) - b->max_time_fixed;
    };

    // we use computational cost to break ties.
    auto fun4 = [&](auto &a, auto &b) {
      auto sa = score_with_expands(a);
      auto sb = score_with_expands(b);
      if (sa == sb) {
        return -a->max_time_fixed < -b->max_time_fixed;
      } else {
        return sa < sb;
      }
    };

    // fifo best first queue
    auto it = std::min_element(open_nodes.begin(), open_nodes.end(), fun2);

    auto min = *it;
    open_nodes.erase(it);
    return min;
  }

  void search() {

    int it = 0;
    int expand_counter = 0;
    while (it++ < max_it) {

      Compu_node *node = pop_and_erase_best();

      if (interactive) {
        node->highlight = true;
        std::ofstream file("g.dot");
        tree_to_graphviz(file, root);
        file.close();
        system_s("dot -Tpdf g.dot -o g.pdf");
        if (open_viewer_first_time && it == 1)
          system_s("zathura g.pdf & ");
        std::cin.get();
        node->highlight = false;
      }

      std::cout << "Chosen node is" << std::endl;
      node->write(std::cout);
      node->expands.push_back(expand_counter++);

      bool recompute_heuristic = false;
      auto new_nodes = node->expand(&recompute_heuristic);

      if (node->type == TYPE::with_variables) {
        std::cout << "node has variables, so I add it again" << std::endl;
        open_nodes.push_back(node);
      }

      for (auto &n : new_nodes) {
        if (n->terminal && n->type == TYPE::with_values) {
          goals.push_back(n);
          std::cout << "We have found a solution" << std::endl;
          rai::String plan = n->node->getTreePathString('\n');
          std::cout << "plan is: " << plan << std::endl;

          if (block_found_plans) {
            std::cout << "lets block found plan -- recomputing the heuristic";
            for (auto &n : open_nodes) {
              n->tabu_plans.push_back(plan);
              n->compute_heuristic();
            }
            std::cout << "I don't add the other child nodes to the graph"
                      << std::endl;
            break;
          }
        }
        open_nodes.push_back(n);
      }

      if (recompute_heuristic) {
        std::cout << "recomputing heuristics on all nodes" << std::endl;
        for (auto &n : open_nodes) {
          n->compute_heuristic();
        }
      }

      if (interactive) {
        node->highlight = true;
        std::ofstream file("g.dot");
        tree_to_graphviz(file, root);
        file.close();
        system_s("dot -Tpdf g.dot -o g.pdf");
        std::cin.get();
        node->highlight = false;
      }
    }
  }
};

BOOST_AUTO_TEST_CASE(compu_node) {

  rai::String folFile = "fol_lab_bench_only_one.g";
  rai::String confFile = "models/lab_setting_bench_only_one.g";

  rai::String goal = "(on goal1_table block1) (on block1 block2)";
  rai::String plan = "(pick block1 block1_ref r_gripper)\n(place block1 "
                     "r_gripper goal1_table)\n(pick block2 block2_ref "
                     "r_gripper)\n(place block2 r_gripper block1)";

  bool visualize = true;
  rai::Configuration C;
  C.addFile(confFile);
  if (visualize)
    C.watch(true);

  rai::LGP_Tree lgp(C, folFile);
  lgp.fol.addTerminalRule(goal);

  // String getTreePathString(char sep=' ') const;

  // rai::String pp = "(pick block1 block1_ref r_gripper)\n(place block1 "
  //                  "r_gripper goal1_table)\n(pick block2 block2_ref
  //                  r_gripper)";

  // auto node = lgp.walkToNode(plan);
  // std::cout << "path is " << std::endl;
  // std::cout << node->getTreePathString('\n') << std::endl;

  // rai::String plan_past_goal =
  //     "(pick block1 block1_ref r_gripper)\n(place block1 "
  //     "r_gripper goal1_table)\n(pick block2 block2_ref "
  //     "r_gripper)\n(place block2 r_gripper block1)\n"
  //     "(pick block2 block1 r_gripper)";

  // std::cout << "walk past goal" << std::endl;
  // auto node = lgp.walkToNode(plan);
  // node->expand();
  // for (auto &c : node->children) {
  //   c->write(std::cout);
  // }
  // std::cout << node->isTerminal << std::endl;
  // throw -1;

  auto root = lgp.root;
  std::cout << "lgp node " << std::endl;
  root->write(std::cout);

  Compu_node croot;
  croot.node = root;
  croot.tree = &lgp;

  std::cout << "heuristic for root" << std::endl;
  croot.compute_heuristic();

  // croot.expand();

  Compu_tree compu_tree;
  compu_tree.max_it = 50;
  compu_tree.root = &croot;
  compu_tree.open_nodes = {&croot};
  compu_tree.interactive = true;
  compu_tree.search();
  std::ofstream file("g.dot");
  tree_to_graphviz(file, &croot);
}

BOOST_AUTO_TEST_CASE(compu_node2) {
  // TODO: double check wether handover is possible or not!

  rai::String folFile = "fol_lab_bench_two_one_object.g";
  rai::String confFile = "./lab_setting_two_two_objects.g";

  rai::String goal = "(on goal2_table block1) (on r_gripper block2)";

  bool visualize = true;
  rai::Configuration C;
  C.addFile(confFile);
  if (visualize)
    C.watch(true);

  rai::LGP_Tree lgp(C, folFile);
  lgp.fol.addTerminalRule(goal);

  // String getTreePathString(char sep=' ') const;

  // rai::String pp = "(pick block1 block1_ref r_gripper)\n(place block1 "
  //                  "r_gripper goal1_table)\n(pick block2 block2_ref
  //                  r_gripper)";

  // auto node = lgp.walkToNode(plan);
  // std::cout << "path is " << std::endl;
  // std::cout << node->getTreePathString('\n') << std::endl;

  // rai::String plan_past_goal =
  //     "(pick block1 block1_ref r_gripper)\n(place block1 "
  //     "r_gripper goal1_table)\n(pick block2 block2_ref "
  //     "r_gripper)\n(place block2 r_gripper block1)\n"
  //     "(pick block2 block1 r_gripper)";

  // std::cout << "walk past goal" << std::endl;
  // auto node = lgp.walkToNode(plan);
  // node->expand();
  // for (auto &c : node->children) {
  //   c->write(std::cout);
  // }
  // std::cout << node->isTerminal << std::endl;
  // throw -1;

  auto root = lgp.root;
  std::cout << "lgp node " << std::endl;
  root->write(std::cout);

  Compu_node croot;
  croot.node = root;
  croot.tree = &lgp;

  system_s("rm -rf " + croot.conflict_folder + " && mkdir " +
           croot.conflict_folder);
  system_s("rm -rf " + croot.path_feasible_samples + " && mkdir " +
           croot.path_feasible_samples);

  croot.compute_heuristic();

  // fixing values for the root.
  //
  auto [graph, factored_nlp] = get_nlps_and_graph(croot.node);
  factored_nlp->report(cout, 3);

  for (auto &v : factored_nlp->__variableIndex)
    std::cout << v << std::endl;

#if 0
  {
  std::cout << "writing clean " << std::endl;
  std::ofstream filex("tmp/g.dot");
  write_graphviz_easy(filex, graph, make_label_writer3(graph));
  filex.close();
  system_s("dot -Tpdf tmp/g.dot -o tmp/g.pdf");
  system_s("zathura tmp/g.pdf & ");
  std::cin.get();
  }
#endif

  factored_nlp->subSelect({}, {});
  auto out = MP_Solver().setProblem(factored_nlp).solve();

  std::map<rai::String, arr> new_values;
  for (auto &v : factored_nlp->__variableIndex) {
    std::cout << v << std::endl;
    new_values[v.nameID] = v.getValue(true);
  }

  croot.values = new_values;

  Compu_tree compu_tree;
  compu_tree.max_it = 50;
  compu_tree.root = &croot;
  compu_tree.open_nodes = {&croot};
  compu_tree.interactive = true;
  compu_tree.search();
  std::ofstream file("g.dot");
  tree_to_graphviz(file, &croot);
}

BOOST_AUTO_TEST_CASE(plan_contains_states) {

  rai::String folFile = "fol_lab_bench_two_one_object.g";
  rai::String confFile = "./lab_setting_two_two_objects.g";

  rai::String goal = "(on goal2_table block1) (on r_gripper block2)";
  rai::String plan = "(pick block1 block1_ref r_gripper)\n(place block1 "
                     "r_gripper goal2_table)\n(pick block2 block2_ref "
                     "r_gripper)";

  bool visualize = false;
  rai::Configuration C;
  C.addFile(confFile);
  if (visualize)
    C.watch(true);

  rai::LGP_Tree lgp(C, folFile);
  lgp.fol.addTerminalRule(goal);

  auto node = lgp.walkToNode(plan);
  auto state_sequence = get_state_sequence(node);

  std::cout << "printing is subset " << std::endl;
  for (auto &s : state_sequence) {
    std::cout << is_subset(s, {"(on r_gripper block2)"}) << std::endl;
  }

  std::cout << "printing is subset " << std::endl;
  for (auto &s : state_sequence) {
    std::cout << is_subset(s, {"(on r_gripper block1)"}) << std::endl;
  }

  BOOST_CHECK(is_subsequence(state_sequence, {{"(on r_gripper block1)"}}));
  BOOST_CHECK(is_subsequence(state_sequence, {{"(on r_gripper block2)"}}));

  BOOST_CHECK(!is_subsequence(state_sequence, {{"(on l_gripper block2)"}}));
  BOOST_CHECK(!is_subsequence(state_sequence, {{"(on l_gripper block1)"}}));

  BOOST_CHECK(is_subsequence(
      state_sequence, {{"(on block1_ref block1)"}, {"(on r_gripper block1)"}}));

  BOOST_CHECK(!is_subsequence(
      state_sequence, {{"(on block1_ref block1)"}, {"(on r_gripper block2)"}}));

  BOOST_CHECK(is_subsequence(
      state_sequence, {{"(on block2_ref block2)"}, {"(on r_gripper block2)"}}));

  BOOST_CHECK(is_subsequence(
      state_sequence, {{"(on block1_ref block1)", "(on block2_ref block2)"},
                       {"(on r_gripper block1)"}}));

  BOOST_CHECK(!is_subsequence(
      state_sequence, {{"(on block2_ref block2)", "(on block1_ref block1)"},
                       {"(on r_gripper block2)"}}));

  {
    std::string _conflict =
        "Atom on(block2_ref, block2)\n\nAtom on(l_gripper, block2)\n";
    std::ofstream fout("conflict.tmp");
    fout << _conflict;
    fout.close();

    auto komo_style = komo_state_sequence_from_PDDL("conflict.tmp");

    BOOST_CHECK(!is_subsequence(state_sequence, komo_style));
  }
  {
    std::string _conflict =
        "Atom on(block2_ref, block2)\n\nAtom on(r_gripper, block2)\n";
    std::ofstream fout("conflict.tmp");
    fout << _conflict;
    fout.close();

    auto komo_style = komo_state_sequence_from_PDDL("conflict.tmp");

    BOOST_CHECK(is_subsequence(state_sequence, komo_style));
  }
}

BOOST_AUTO_TEST_CASE(extract_and_match_subgraphs) {
  // working: I can extract "pick" and "place".
  // Lets use this in the tree to try to reuse.

  rai::String folFile = "fol_lab_bench_two_one_object.g";
  rai::String confFile = "./lab_setting_two_two_objects.g";

  rai::String goal = "(on goal2_table block1) (on r_gripper block2)";
  rai::String plan = "(pick block1 block1_ref r_gripper)\n(place block1 "
                     "r_gripper goal2_table)\n(pick block2 block2_ref "
                     "r_gripper)";

  bool visualize = false;
  rai::Configuration C;
  C.addFile(confFile);
  if (visualize)
    C.watch(true);

  rai::LGP_Tree lgp(C, folFile);
  lgp.fol.addTerminalRule(goal);

  auto node = lgp.walkToNode(plan);

  auto bound = rai::BD_seq;
  bool verbose = false;

  StringAA collisions_list = {};
  StringAA collisions_pairwise = {};
  node->prepareProblem(bound, verbose, collisions_list, collisions_pairwise);

  auto factored_nlp = node->problem(bound).komo->mp_FineStructured();

  auto new_variableIndex =
      break_chain(factored_nlp->__variableIndex, separate_by_name,
                  [](auto &s) { return true; });

  factored_nlp->__variableIndex = new_variableIndex;
  factored_nlp->recompute_frame2ids_and_features();
  factored_nlp->komo.run_prepare(0);
  factored_nlp->report(cout, 3);

  arr default_x = factored_nlp->komo.x.copy();

  BGraph graph = create_boost_graph(*factored_nlp, false);

  {
    std::cout << "writing graph " << std::endl;
    std::ofstream file("tmp/g.dot");
    write_graphviz_easy(file, graph, make_label_writer3(graph));
    file.close();
    system_s("dot -Tpdf tmp/g.dot -o tmp/g.pdf");
    std::cin.get();
  }

  auto out = extract_interesting_subgraphs(graph);

  // lets write the interesting graphs into a folder.

  std::string path = "tmp_graphs/";
  std::string command = "mkdir -p " + path;
  std::string command2 = "mkdir -p " + path;
  system_s(command.c_str());
  system_s(command2.c_str());

  int id = 0;
  for (auto &g : out) {
    ofstream file_graphviz(path + "/" + std::to_string(id) + ".dot");
    write_graphviz_easy(file_graphviz, g, make_label_writer3(g));
    ofstream file_graphboost(path + "/" + std::to_string(id) + ".dat");
    file_graphboost << boost::write(g);
    id++;
  }
}

BOOST_AUTO_TEST_CASE(match_with_check) {

  // What does it mean to check?
  // Let's say I do check against variables and last set state.
  // Now I have to check that it fits. How?

  // For pick, the pose of the objects should be the same.
  // For fixed variables, I have to check that the numeric value is the same.

  // Issue: what happens with the constant variables? For example, I if set
  // on the table or I
}

// copy a var manually

void propagate_a_var(BGraph &graph, const std::string &name, int from) {
  // continue here!!

  auto [beg_ref, end_ref] = boost::vertices(graph);
  bool found = false;
  for (auto it = beg_ref; it != end_ref && !found; it++) {
    auto &v = graph[*it];
    if (__startsWith(v.name, name) && v.time == from && v.x.size()) {
      for (auto &it2 = beg_ref; it2 != end_ref && !found; it2++) {
        auto &v2 = graph[*it2];
        if (__startsWith(v2.name, name) && v2.time == from + 1 &&
            !v2.x.size()) {
          found = true;
          v2.x = v.x;
          std::cout << "propagation successful" << std::endl;
        }
      }
    }
  }
  if (!found) {
    std::cout << "warning: could not propagate" << std::endl;
  }
  //
}

BOOST_AUTO_TEST_CASE(propagate) {

  // Fix variables that are set
  rai::String folFile = "fol_lab_bench_two_one_object.g";
  rai::String confFile = "./lab_setting_two_two_objects.g";

  rai::String goal = "(on goal2_table block1) (on r_gripper block2)";
  // rai::String plan =
  //     "(pick block1 block1_ref r_gripper)\n(place block1 r_gripper
  //     block2)";
  rai::String plan = "(pick block1 block1_ref r_gripper)\n(pick block2 "
                     "block2_ref l_gripper)";

  bool visualize = false;
  rai::Configuration C;
  C.addFile(confFile);
  if (visualize)
    C.watch(true);

  rai::LGP_Tree lgp(C, folFile);
  lgp.fol.addTerminalRule(goal);

  auto node = lgp.walkToNode(plan);

  auto bound = rai::BD_seq;
  bool verbose = false;

  StringAA collisions_list = {};
  StringAA collisions_pairwise = {};
  node->prepareProblem(bound, verbose, collisions_list, collisions_pairwise);

  auto factored_nlp = node->problem(bound).komo->mp_FineStructured();

  auto new_variableIndex =
      break_chain(factored_nlp->__variableIndex, separate_by_name,
                  [](auto &s) { return true; });

  factored_nlp->__variableIndex = new_variableIndex;
  factored_nlp->recompute_frame2ids_and_features();
  factored_nlp->komo.run_prepare(0);
  factored_nlp->report(cout, 3);

  arr default_x = factored_nlp->komo.x.copy();

  BGraph graph = create_boost_graph(*factored_nlp, false);

  // lets solve the first state

  Solve_nlp_with_counter solver;
  BGraph graph_slice;
  int i = 0;
  boost::copy_graph(
      make_filtered_graph(graph, boost::keep_all(),
                          Keep_FF<BGraph>(
                              &graph,
                              [&](const GBNode &s) { return s.time == 0; },
                              [&](auto &s) { return true; }, true)),
      graph_slice);

  BOOST_ASSERT(
      solver.solve_nlp(factored_nlp, default_x, graph_slice, visualize, true));

  // fix the variables in the original one

  auto [beg, end] = boost::vertices(graph_slice);
  auto [beg_ref, end_ref] = boost::vertices(graph);
  for (auto it = beg; it != end; it++) {
    // look for the same vertex in the parent graph
    auto &v = graph_slice[*it];
    if (!v.is_var)
      continue;
    auto it_in_ref = std::find_if(beg_ref, end_ref, [&](const auto &s) {
      return graph[s].is_var && graph[s].id == v.id;
    });
    Qassert(it_in_ref != end_ref);

    arr x = factored_nlp->__variableIndex(v.id).getValue(true);
    graph[*it_in_ref].x = std::vector(x.begin(), x.end());
  }

  propagate_fixed_vars(graph);

  {
    std::cout << "writing graph " << std::endl;
    std::ofstream file("tmp/g.dot");
    write_graphviz_easy(file, graph, make_label_writer3(graph));
    file.close();
    system_s("dot -Tpdf tmp/g.dot -o tmp/g.pdf");
    std::cin.get();
  }

  {
    std::cout << "writing graph " << std::endl;
    std::ofstream file("tmp/g.dot");
    write_graphviz_easy(file, graph, make_label_writer3(graph));
    file.close();
    system_s("dot -Tpdf tmp/g.dot -o tmp/g.pdf");
    std::cin.get();
  }
}

BOOST_AUTO_TEST_CASE(reuse_samples) {

  // I want to do matches that fulfil the complete time steps.
  //
  // Let's Just Say: - number time steps of the time slice has
  // to be equal to the total time slice.

  rai::String folFile = "fol_lab_bench_two_one_object.g";
  rai::String confFile = "./lab_setting_two_two_objects.g";
  rai::String plan1 = "(pick block1 block1_ref r_gripper)";
  rai::String plan2 = "(pick block1 block1_ref r_gripper)";
  rai::String plan3 = "(pick block2 block2_ref r_gripper)"
                      "(place block2 r_gripper goal2_table)"
                      "(pick block1 block1_ref r_gripper)";

  bool visualize = false;
  rai::Configuration C;
  C.addFile(confFile);

  if (visualize)
    C.watch(true);

  rai::LGP_Tree lgp(C, folFile);
  auto node = lgp.walkToNode(plan1);
  auto [graph, factored_nlp] = get_nlps_and_graph(node);

  arr default_x = factored_nlp->komo.x.copy();
  // lets solve the first state

  Solve_nlp_with_counter solver;
  BOOST_CHECK(
      solver.solve_nlp(factored_nlp, default_x, graph, visualize, true));

  get_values_from_nlp(*factored_nlp, graph);

  auto interesting_subgraphs = extract_interesting_subgraphs(graph);

  // save graphs
  std::string path = "tmp_graphs/";
  std::string command = "rm -rf " + path;
  std::string command2 = "mkdir -p " + path;
  system_s(command.c_str());
  system_s(command2.c_str());

  int id = 0;
  for (auto &g : interesting_subgraphs) {
    ofstream file_graphviz(path + "/" + std::to_string(id) + ".dot");
    write_graphviz_easy(file_graphviz, g, make_label_writer3(g));
    ofstream file_graphboost(path + "/" + std::to_string(id) + ".dat");
    file_graphboost << boost::write(g);
    id++;
  }

  // extract

#if 0
  {
    std::cout << "now lets take plan 2" << std::endl;
    auto node2 = lgp.walkToNode(plan2);
    auto [graph2, factored_nlp2] = get_nlps_and_graph(node2);

    const arr default_x = factored_nlp2->komo.x;

    Solve_nlp_with_counter solver2;
    BGraph graph_slice;
    boost::copy_graph(
        make_filtered_graph(graph2, boost::keep_all(),
                            Keep_FF<BGraph>(
                                &graph2,
                                [&](const GBNode &s) { return s.time == 0; },
                                [&](auto &s) { return true; }, true)),
        graph_slice);
    BOOST_ASSERT(solver.solve_nlp(factored_nlp2, default_x, graph_slice,
                                  visualize, true));


    get_values_from_nlp(*factored_nlp2, graph_slice);
    fix_subgraph(graph2, graph_slice);

    // check which samples are applicable
    // I could also associate each sample to a transtion.
    std::vector<GraphStruct> graphs = load_graphs(path.c_str());

    //
    

    std::vector<std::vector<std::pair<int, int>>> matches;
    BGraph chosen;
    bool match_vars, match_fix_vars;

    for (auto &g_candidate : graphs) {
      auto &graph_small = g_candidate.g;
      auto &graph_big = graph2;
      matches.clear();
      subgraph_with_distance(graph_small, graph_big, matches, match_vars,
                             match_fix_vars);
      BOOST_CHECK(match_vars);
      BOOST_CHECK(match_fix_vars);
      if (match_vars && match_fix_vars) {
        chosen = graph_small;
      }
      break;
    }

    std::cout << "copying values " << std::endl;
    auto &graph_big = graph2;
    auto &graph_small = chosen;
    auto &match = matches.front();

    for (auto &[a, b] : match) {
      // a : small, b : big
      if (graph_small[a].is_var && !graph_big[b].x.size() ){
          graph_big[b].x = graph_small[a].x;
          if (false)
          {
            std::cout << "copy values of" << std::endl;
            std::cout << "***\n";
            std::cout << graph_small[a] << std::endl;
            std::cout << "--\n";
            std::cout << graph_big[b] << std::endl;
            std::cout << "***\n";
          }
        }
    }


    factored_nlp2->reset_subselect(default_x);

    std::cout << "Complet problem" << std::endl;
    factored_nlp2->report(std::cout, 5);

    std::cout << "propagating variables" << std::endl;
    propagate_fixed_vars(graph_big);
    propagate_a_var(graph_big, "l_panda", 0);

    bool fix_vars=true;
    warmstart_nlp_from_graph(graph_big, factored_nlp2,fix_vars);

    if (!factored_nlp2->subVars.N) {
      std::cout << "no free variables, only checking" << std::endl;
      factored_nlp2->subSelect({}, {});
      factored_nlp2->report(std::cout, 5);
      double tol = 1e-1;
      bool feasible = t_check_feasible(*factored_nlp2, tol);
      BOOST_CHECK(feasible);
      std::cout << "feasible is: " << feasible << std::endl;
    } else {
      std::cout << "still free variables, fixing" << std::endl;
      auto out = MP_Solver().setProblem(factored_nlp2).solve();
    }

    std::cout << "Done" << std::endl;
    factored_nlp2->report(std::cout, 5);
  }
#endif

#if 1
  {
    Solve_nlp_with_counter solver3;
    BGraph graph_slice3;

    std::cout << "now lets take plan 2" << std::endl;
    auto node3 = lgp.walkToNode(plan3);
    auto [graph3, factored_nlp3] = get_nlps_and_graph(node3);
    const arr default_x = factored_nlp3->komo.x;

    boost::copy_graph(
        make_filtered_graph(graph3, boost::keep_all(),
                            Keep_FF<BGraph>(
                                &graph3,
                                [&](const GBNode &s) { return s.time <= 2; },
                                [&](auto &s) { return true; }, true)),
        graph_slice3);
    BOOST_ASSERT(solver.solve_nlp(factored_nlp3, default_x, graph_slice3,
                                  visualize, true));

    get_values_from_nlp(*factored_nlp3, graph_slice3);
    fix_subgraph(graph3, graph_slice3);

    std::cout << "writing clean " << std::endl;
    std::ofstream file("tmp/g.dot");
    write_graphviz_easy(file, graph3, make_label_writer3(graph3));
    file.close();
    system_s("dot -Tpdf tmp/g.dot -o tmp/g.pdf");
    // system_s("zathura tmp/g.pdf & ");
    std::cin.get();

    std::vector<GraphStruct> graphs = load_graphs(path.c_str());
    std::vector<std::vector<std::pair<int, int>>> matches;
    BGraph chosen;
    bool match_vars, match_fix_vars;
    for (auto &g_candidate : graphs) {
      auto &graph_small = g_candidate.g;
      auto &graph_big = graph3;

      bool verbose = true;
      subgraph_with_distance(graph_small, graph_big, matches, match_vars,
                             match_fix_vars, verbose);
      BOOST_CHECK(match_vars);
      BOOST_CHECK(match_fix_vars);
      if (match_fix_vars && match_vars) {
        chosen = graph_small;
        break;
      }
    }

    auto &graph_big = graph3;
    auto &graph_small = chosen;
    auto &match = matches.front();

    {
      std::ofstream file("tmp/g2_.dot");
      write_graphviz_easy(file, graph3, make_label_writer3(graph3));
    }

    for (auto &[a, b] : match) {
      if (graph_small[a].is_var && !graph_big[b].x.size())
        graph_big[b].x = graph_small[a].x;
    }

    factored_nlp3->reset_subselect(default_x);

    std::cout << "Complet problem" << std::endl;
    factored_nlp3->report(std::cout, 5);

    std::cout << "graphviz-1" << std::endl;
    write_graphviz_easy(std::cout, graph3, make_label_writer3(graph3));

    {
      std::ofstream file("tmp/g2.dot");
      write_graphviz_easy(file, graph3, make_label_writer3(graph3));
      file.close();
      system_s("dot -Tpdf tmp/g2.dot -o tmp/g2.pdf");
    }

    std::cout << "propagating variables manually" << std::endl;
    propagate_a_var(graph_big, "l_panda", 2);
    propagate_fixed_vars(graph_big);

    std::cout << "graphviz-2" << std::endl;
    write_graphviz_easy(std::cout, graph3, make_label_writer3(graph3));

    bool fix_vars = true;
    warmstart_nlp_from_graph(graph_big, factored_nlp3, fix_vars);

    if (!factored_nlp3->subVars.N) {
      std::cout << "no free variables, only checking" << std::endl;
      factored_nlp3->subSelect({}, {});
      factored_nlp3->report(std::cout, 5);
      double tol = 1e-1;
      bool feasible = t_check_feasible(*factored_nlp3, tol);
      BOOST_CHECK(feasible);
      std::cout << "feasible is: " << feasible << std::endl;
    } else {
      std::cout << "still free variables, fixing" << std::endl;
      auto out = MP_Solver().setProblem(factored_nlp3).solve();
    }

    std::cout << "Done" << std::endl;
    factored_nlp3->subSelect({}, {});
    factored_nlp3->report(std::cout, 5);

    // lets do a warmstart
  }

#endif
}

BOOST_AUTO_TEST_CASE(grasp_models) {
  // lets test the quim grasp model...

  rai::String folFile = "fol_lab_bench_two_one_object.g";
  rai::String confFile = "./lab_setting_two_two_objects.g";
  rai::String plan1 = "(pick block1 block1_ref r_gripper)";
  rai::String plan2 = "(pick block1 block1_ref r_gripper)";
  rai::String plan3 = "(pick block2 block2_ref r_gripper)"
                      "(place block2 r_gripper goal2_table)"
                      "(pick block1 block1_ref r_gripper)";

  bool visualize = false;
  rai::Configuration C;
  C.addFile(confFile);

  if (visualize)
    C.watch(true);

  rai::LGP_Tree lgp(C, folFile);
  auto node = lgp.walkToNode(plan1);
  auto [graph, factored_nlp] = get_nlps_and_graph(node);

  arr default_x = factored_nlp->komo.x.copy();

  Solve_nlp_with_counter solver;
  BOOST_CHECK(
      solver.solve_nlp(factored_nlp, default_x, graph, visualize, true));

  // visualize

  factored_nlp->komo.view(true);
  factored_nlp->komo.view_play(true);

  std::cout << "writing clean " << std::endl;
  std::ofstream file("tmp/g.dot");
  write_graphviz_easy(file, graph, make_label_writer3(graph));
  file.close();
  system_s("dot -Tpdf tmp/g.dot -o tmp/g.pdf");
  // system_s("zathura tmp/g.pdf & ");
  std::cin.get();
}

BOOST_AUTO_TEST_CASE(small_vs_big_table) {

  // lets test the quim grasp model...

  // i have to put N objects in a table.
  // See if I find a diference based on the size of the table. 
  // use only one robot

}


BOOST_AUTO_TEST_CASE(two_robot_move) {

  // check if the intermediate position is constrained or not, based on the 
  // the table size.
  // How to warmstart each variable efficiently?

}

BOOST_AUTO_TEST_CASE(new_warmstart) {

  // check the warmstart of MARC.

}





