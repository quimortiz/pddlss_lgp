

#include "KOMO/local_utils.hpp"
#include <boost/program_options.hpp>

#include "Core/util.h"
#include "Core/util.ipp"
#include "KOMO/conv_komo_graph.h"
#include "KOMO/komo.h"
#include "LGP/LGP_symconflict.h"
#include "LGP/LGP_tree.h"
#include "PathAlgos/ConfigurationProblem.h"
#include "PathAlgos/RRT_PathFinder.h"
#include "magic_enum.hpp"
#include <LGP/LGP_graph.h>
#include <boost/config.hpp> // put this first to suppress some VC++ warnings
#include <filesystem>
#include <queue>

std::ofstream out_file;
int num_threads;

// we love good old global variables
bool VISUALIZE_KOMO = false;
std::string CONFLICT_FOLDER = "__conflicts";
bool INTERACTIVE_INTERESTING_SUBGRAPHS = false;
std::string PATH_FEASIBLE_SAMPLES = "tmp_feasible_graphs";
std::string PATH_TO_DOWNWARD_FOLDER = "/home/quim/stg/lgp-pddl/downward/";

std::vector<double> PANDA_Q0 = {0., -.5, 0., -2, -0., 2., -.5};

using Match = std::vector<std::pair<int, int>>;

rai::Transformation random_placement_on_table(const arr &table_size) {

  std::vector<double> lb(3, 0.);
  std::vector<double> ub(3, 0.);

  lb.at(0) = -table_size(0) / 2.;
  lb.at(1) = -table_size(1) / 2.;
  ub.at(0) = +table_size(0) / 2.;
  ub.at(1) = +table_size(1) / 2.;
  lb.at(2) = .1;
  ub.at(2) = .1;

  assert(lb.size() == 3);
  assert(ub.size() == 3);
  std::vector<double> out(4);

  for (size_t i = 0; i < 3; i++) {
    assert(lb.at(i) <= ub.at(i));
    double r = (double(rand()) / RAND_MAX);
    out.at(i) = lb.at(i) + r * (ub.at(i) - lb.at(i));
  }
  double yaw = -180 + (double(rand()) / RAND_MAX) * 360;

  out.at(3) = yaw;
  rai::Vector v = rai::Vector(out.at(0), out.at(1), out.at(2));
  rai::Quaternion q;
  q.setRadZ(yaw / 180. * M_PI);
  return rai::Transformation(v, q);
}

void inline create_dir_if_necessary(const char *file) {
  const std::filesystem::path path = std::filesystem::path(file).parent_path();
  if (!path.empty()) {
    std::filesystem::create_directories(path);
  }
}

inline std::string gen_random_id(const int len = 6) {
  static const char alphanum[] = "0123456789"
                                 "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                 "abcdefghijklmnopqrstuvwxyz";
  std::string tmp_s;
  tmp_s.reserve(len);

  for (int i = 0; i < len; ++i) {
    tmp_s += alphanum[rand() % (sizeof(alphanum) - 1)];
  }

  return tmp_s;
}

inline std::string get_time_stamp() {

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);

  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y--%H-%M-%S");
  auto str = oss.str();

  return str;
}

template <typename T> bool in_vector(const T &e, const std::vector<T> &v) {
  return std::find(v.begin(), v.end(), e) != v.end();
}

std::string clean_name(const std::string &name) {

  std::size_t found = name.find_first_not_of(".");
  if (found == std::string::npos) {
    std::cout << "name not working" << std::endl;
    throw std::runtime_error("BUU");
  }
  return name.substr(found);
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
    std::cout << "conditional vars are" << conditional_vars << std::endl;
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
                              auto &name = s.name;
                              auto name_clean = clean_name(name);

                              std::vector<std::string> del = {
                                  "F_qQuaternionNorms",
                                  "F_PositionRel/0-l_gripper-l_panda_base",
                                  "F_PositionRel/0-r_gripper-r_panda_base",
                                  "F_Pose/1-block", "F_Pose/2-block"};
                              return std::find_if(
                                         del.begin(), del.end(), [&](auto &j) {
                                           return _startsWith(name_clean, j);
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
      if (!v.is_var && in_vector(clean_name(v.name), name_ref)) {
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
    boost::copy_graph(make_filtered_graph(
                          graph, boost::keep_all(),
                          Keep_FF<BGraph>(
                              &graph, [](auto &s) { return true; },
                              [&](auto &s) {
                                return in_vector(clean_name(s.name), name_ref);
                              },
                              false, true)),
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

auto get_nlps_and_graph(rai::LGP_Node *node,
                        const StringAA &collisions_list = {},
                        const StringAA &collisions_pairwise = {}) {

  auto bound = rai::BD_seq;
  bool verbose = false;

  node->prepareProblem(bound, verbose, collisions_list, collisions_pairwise);
  auto factored_nlp = node->problem(bound).komo->mp_FineStructured();

  std::cout << "original factored" << std::endl;
  factored_nlp->subSelect({}, {});
  factored_nlp->report(std::cout, 3);

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

  if (!conflicts.size()) {
    std::cout << "ERROR: "
              << "conflicts.size() " << conflicts.size() << std::endl;
  }
  if (conflicts.size()) {
    auto &cconflict = conflicts.front();
    uintA vars_sorted, cons_sorted;
    boost::tie(vars_sorted, cons_sorted) = vars_cons_from_graph(cconflict);

    std::cout << "found infeasible relaxation pddlss_gnlpp" << std::endl;
    Relaxation infeas_relaxation = {
        0, 0, {false, vars_sorted, cons_sorted, {}}};

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
}

struct Compu_node {
  // lets try to order: minimize the time to reach the goal.
  // with values: heuristic, never added again
  // with variables: heuristic + expands, added again. I only do the
  // symbolic expansion once.
  // root is of type: with values

  Compu_node() { id = ++global_node_id; }

  int limit_infeas_childs = rai::getParameter<int>("ct/limit_infeas_childs", 5);
  int limit_num_free_vars_steps =
      rai::getParameter<int>("ct/limit_num_free_vars_steps", 10);
  bool only_optimize_goals =
      rai::getParameter<bool>("ct/only_optimize_goals", false);
  bool check_conflicts = rai::getParameter<int>("ct/check_conflicts", 1);
  bool fix_vars = true;
  bool fixed_non_used_robot = false;
  const int BIG_NUMBER_COST_TO_GO = 1000;
  const double FEASIBLE_THRESHOLD = .1;

  bool use_database_samples = true;

  bool feasible = true;
  int id;
  int num_infeas_childs = 0;
  rai::LGP_Node *node; // a Compu Node is a LGP Node with some
                       // assigned values!
  rai::Configuration *C;
  rai::LGP_Tree *tree;                   // or fol?
  std::vector<std::string> logic_states; // other representation?
  BGraph graph;
  std::map<std::string, arr> values;
  std::map<std::string, arr> new_values;
  TYPE type = TYPE::with_variables;
  Compu_node *parent = nullptr;
  std::vector<Compu_node *> childs;
  bool symbolic_expansion_done = false; // it can only be done once
  int depth =
      0; // takes into account the num of expands operations in the parent
  int tree_depth = 0; // does not consider this
  int heuristic = 0;
  int path_from_root = 0;
  bool terminal = false;
  int max_time_fixed = 0;
  std::vector<int> matches_in_child_assignment = {};
  int subgraph_matched = -1; // indicates which subgraph I have used for fixing
                             // variables. -1 means none.

  StringAA collisions_list = {};
  StringAA collisions_pairwise = {};

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
        path_to_node_has_conflict(node, CONFLICT_FOLDER)) {

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
        compute_pddl_heuristic(node, PATH_TO_DOWNWARD_FOLDER, CONFLICT_FOLDER);

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

  std::vector<std::map<std::string, arr>> get_compu_history() {
    std::vector<std::map<std::string, arr>> out;
    std::vector<Compu_node *> path = path_to_root();
    // only consider the ones with values

    std::vector<Compu_node *> path_only_value;

    std::copy_if(path.begin(), path.end(), std::back_inserter(path_only_value),
                 [](const auto &c) { return c->type == TYPE::with_values; });

    for (size_t i = 0; i < path_only_value.size(); i++) {
      auto node = path_only_value.at(path_only_value.size() - 1 - i);
      out.push_back(node->new_values);
    }
    out.push_back(new_values);
    return out;
  }

  void compute_values(bool *recompute_heuristic) {
    // recompute_heuristic is a flag to signal that I have
    // to update the PDDL model

    if (type == TYPE::with_values) {
      throw std::runtime_error("values are already computed");
    }
    type = TYPE::with_values;

    // collision list

    // StringAA collision_list = {{"block1"}, {"block2"}, {"block3"},
    //                            {"block4"}, {"block5"}, {"block6"}};

    std::cout << "check that this names are fine!! " << std::endl;
    StringAA collision_list = {{"f_block1_col"}, {"f_block2_col"},
                               {"f_block3_col"}, {"f_block4_col"},
                               {"f_block5_col"}, {"f_block6_col"}};

    // check
    StringAA collision_list_here = {};

    for (auto &c : collision_list) {
      auto &C = node->tree.kin;
      if (C.getFrame(c(0)))
        collision_list_here.append(c);
    }

    auto [graph, factored_nlp] = get_nlps_and_graph(node, collision_list_here);
    factored_nlp->report(cout, 3);

    arr default_x = factored_nlp->komo.x.copy();
    auto &vi = factored_nlp->__variableIndex;

    // example: id 0 dim 7 dofs [block1.202 ]  name block1 nameID block1.202
    // time 0
    std::cout << "variables of the problems are " << std::endl;
    for (auto &v : vi)
      std::cout << v << std::endl;

    auto [beg, end] = boost::vertices(graph);

    std::cout << "fixed vars are " << std::endl;
    for (auto &[k, v] : values) {
      std::cout << k << ":" << v << std::endl;
      auto it = std::find_if(vi.begin(), vi.end(), [&](auto &s) {
        return strcmp(s.nameID.p, k.c_str()) == 0;
      });
      Qassert(it != vi.end());
      it->setValue(v);

      auto it_g = std::find_if(beg, end, [&](const auto &s) {
        return graph[s].is_var && graph[s].name == k;
      });
      Qassert(it_g != end);
      graph[*it_g].x = std::vector(v.begin(), v.end());
    }

    auto it = std::max_element(beg, end, [&](const auto &a, const auto &b) {
      return graph[a].time < graph[b].time;
    });

    int max_time = graph[*it].time;
    Qassert(std::distance(beg, end));

    int max_time_conditional;
    if (values.size()) {
      auto it2 = std::max_element(beg, end, [&](const auto &a, const auto &b) {
        return bool(graph[a].x.size()) * graph[a].time <
               bool(graph[b].x.size()) * graph[b].time;
      });

      max_time_conditional = graph[*it2].time;
    } else {
      max_time_conditional = -1;
    }

    std::cout << "max time all and conditional " << max_time << " "
              << max_time_conditional << std::endl;

    std::cout << "propagating fixed vars " << std::endl;
    propagate_fixed_vars(graph);

    if (max_time_conditional >= 0) {

      if (fixed_non_used_robot) {

        std::cout << "fixing only last decision" << std::endl;
        std::stringstream ss;
        ss << *node->decision;
        std::string ss_str = ss.str();
        rai::String plan = node->getTreePathString('\n');
        std::cout << "plan here is " << plan << std::endl;

        std::vector<std::string> actions = string_to_tokens(plan.p);

        Qassert(actions.size() == max_time);

        for (size_t t = max_time_conditional; t < max_time; t++) {

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

      if (use_database_samples) {

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
            load_graphs(PATH_FEASIBLE_SAMPLES.c_str());

        // we have to take only the graph of: [ last fixed, max_time ]

        // i want the same max-min interval
        int it_counter = 0;

        std::vector<std::vector<std::pair<int, int>>> matches;
        BGraph chosen;
        bool match_vars, match_fix_vars;

        int match_index = 0;
        std::cout << "database of graphs " << graphs.size() << std::endl;
        for (auto &g_candidate : graphs) {

          if (in_vector(match_index, matches_in_child_assignment)) {
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
            subgraph_matched = match_index;
            break;
          } else {
            // I clear the vector because maybe there was a symbolic
            // match only.
            matches.clear();
          }
          match_index++;
        }

        if (matches.size()) {
          std::cout << "there is a match" << std::endl;
          for (auto &[a, b] : matches.front()) {
            if (chosen[a].is_var && !graph_slice[b].x.size())
              graph_slice[b].x = chosen[a].x;
          }
          fix_subgraph(graph, graph_slice);
        }
      }
    }

    // HERE I have to put the good warmstart!

    warmstart_nlp_from_graph(graph, factored_nlp, fix_vars);
    // check if there are still free variables

    if (factored_nlp->subFeats.N) {
      shared_ptr<SolverReturn> out;

      // update position of random frames -- random regularization for sampling!

      // get table size

      arr table_size =
          factored_nlp->komo.world.getFrame("goal_table")->getShape().size;

      // std::cout << "table size is " << table_size << std::endl;

      // update the position of objects!

      for (auto &f : factored_nlp->komo.pathConfig.frames) {
        // std::cout << f->name << " ID " << f->ID << std::endl;

        if (f->name == "f1_rand" || f->name == "f2_rand") {
          rai::Transformation T = random_placement_on_table(table_size);
          std::cout << "chosen random placement is " << T << std::endl;
          f->set_Q() = T;
          f->ensure_X();
        }
      }

      std::cout << "DONE" << std::endl;

      out = timed_solve(*factored_nlp).first;
      // out = NLP_Solver().setProblem(factored_nlp).solve();
      // out->write(std::cout);
      std::cout << std::endl;

      feasible = fabs(out->ineq) + fabs(out->eq) < FEASIBLE_THRESHOLD;
      CHECK_EQ(feasible, t_check_feasible(*factored_nlp, FEASIBLE_THRESHOLD),
               "");
      std::cout << "result of optimization:" << feasible << std::endl;
      factored_nlp->subSelect({}, {});
      // There could be small error on previously fixed variables

      std::cout << "check ALL:" << feasible << std::endl;
      if (feasible) {
        // TODO: why is this here? I don't understand the code anymore :/
        CHECK_EQ(feasible,
                 t_check_feasible(*factored_nlp, 2 * FEASIBLE_THRESHOLD), "");
      }

      if (VISUALIZE_KOMO) {
        factored_nlp->komo.view_play(true);
      }

    } else {
      // check all the constraints.
      factored_nlp->subSelect({}, {});
      feasible = t_check_feasible(*factored_nlp, FEASIBLE_THRESHOLD);
      std::cout << "no free variables after warmstart" << std::endl;
      std::cout << "assignment_feasible:" << feasible << std::endl;
    }

    if (!feasible) {
      // expand_symbolic = false;
      std::cout << "CAUTION, optimizaton is infeasible" << std::endl;
      std::cout << "setting expand symbolic to false" << std::endl;
      if (check_conflicts && max_time_conditional == 0) {
        std::cout << "Warning: check conflicts is only working  for no "
                     "condition variables"
                  << std::endl;
        find_conflict(node, factored_nlp, default_x, CONFLICT_FOLDER);
        if (recompute_heuristic)
          *recompute_heuristic = true;
      }

    } else {

      // I fix all the graph
      get_values_from_nlp(*factored_nlp, graph);

      auto interesting_subgraphs = extract_interesting_subgraphs(graph);

      // but I should only consider the newly assigned variables!
      // save to folder

      std::string command2 = "mkdir -p " + PATH_FEASIBLE_SAMPLES;
      system_s(command2.c_str());

      std::cout << "writing intersesting subgraphs" << std::endl;
      for (auto &g : interesting_subgraphs) {
        std::cout << global_id_interesting_subgraphs << std::endl;
        ofstream file_graphviz(PATH_FEASIBLE_SAMPLES + "/" +
                               int_to_string(global_id_interesting_subgraphs) +
                               ".dot");
        write_graphviz_easy(file_graphviz, g, make_label_writer3(g));
        ofstream file_graphboost(
            PATH_FEASIBLE_SAMPLES + "/" +
            int_to_string(global_id_interesting_subgraphs) + ".dat");
        file_graphboost << boost::write(g);

        if (INTERACTIVE_INTERESTING_SUBGRAPHS)
          std::cin.get();
        global_id_interesting_subgraphs++;
      }
    }
    // I have to fix the variables for the child!
    std::map<std::string, arr> all_values;
    for (auto &v : factored_nlp->__variableIndex) {
      std::cout << v << std::endl;
      all_values[v.nameID.p] = v.getValue(true);
    }
    std::cout << "warning: values is overwritten by new_values" << std::endl;

    // TODO: continue here, only with keys
    // CONTINUE HERE
    new_values.clear();

    // The difference of two sets is formed by the elements that are present in
    // the first set, but not in the second one. The elements copied by the
    // function come always from the first range, in the same order. The
    // elements in the both the ranges shall already be ordered.
    //
    //

    std::vector<std::string> all_values_keys(all_values.size());
    std::transform(all_values.begin(), all_values.end(),
                   all_values_keys.begin(), [](auto &kv) { return kv.first; });

    std::vector<std::string> values_keys(values.size());
    std::transform(values.begin(), values.end(), values_keys.begin(),
                   [](auto &kv) { return kv.first; });

    std::vector<std::string> new_keys;

    set_difference(all_values_keys.begin(), all_values_keys.end(),
                   values_keys.begin(), values_keys.end(),
                   std::back_inserter(new_keys));

    values = all_values;

    std::for_each(new_keys.begin(), new_keys.end(), [&](const auto &s) {
      std::cout << s << std::endl;
      auto k = values.at(s);
      new_values.insert({s, k});
    });

    std::cout << "after opti, fixed vars are " << std::endl;
    for (auto &[k, v] : values) {
      std::cout << k << ":" << v << std::endl;
    }

    max_time_fixed = max_time;
  }

  std::vector<Compu_node *> path_to_root() {
    std::vector<Compu_node *> path;
    Compu_node *t_parent = parent;
    while (t_parent) {
      path.push_back(t_parent);
      t_parent = t_parent->parent;
    }
    return path;
  }

  std::vector<Compu_node *> expand(bool *recompute_heuristic) {
    std::vector<Compu_node *> new_nodes{};
    num_expands++;
    bool expand_symbolic = true;
    ; // I only expand symbolic
      // if I have been able to
      // compute variables
    if (num_expands > 1 && type == TYPE::with_values) {
      throw std::runtime_error("node with values can only be expanded once");
    }
    if (type == TYPE::with_variables &&
        (!only_optimize_goals || node->isTerminal)) {

      Compu_node *cc = new Compu_node();
      cc->node = node; // it points to the same node
      cc->tree_depth = tree_depth;
      cc->values = values;
      cc->depth = depth + num_expands - 1; //
                                           // first time it will be = depth,
                                           // then depth+1,...
      cc->tree = tree;
      cc->heuristic = heuristic; // it has the same heuristic
      cc->terminal = terminal;
      // I have to block some subgraphs if they have been solved before.
      cc->matches_in_child_assignment = matches_in_child_assignment;
      cc->compute_values(recompute_heuristic);
      cc->parent = this;
      expand_symbolic = cc->feasible;

      // if the child used a subgrap, i update the matches in child assignment
      if (cc->subgraph_matched != -1)
        matches_in_child_assignment.push_back(cc->subgraph_matched);

      if (!cc->feasible)
        num_infeas_childs += 1;

      childs.push_back(cc);
      new_nodes.push_back(cc); // is it a goal
    }

    if (num_infeas_childs >= limit_infeas_childs) {
      std::cout << "limit infeas childs reached" << std::endl;
      feasible = false;
    }

    if (!symbolic_expansion_done && expand_symbolic && feasible) {

      Qassert(tree_depth >= max_time_fixed);
      if (tree_depth - max_time_fixed <= limit_num_free_vars_steps) {
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
          cc->parent = this;
          childs.push_back(cc);
          new_nodes.push_back(cc);
        };
        symbolic_expansion_done = true;
      } else {
        std::cout << "expansion not allowed! " << std::endl;
        std::cout << "depth:" << depth << " max_time_fixed:" << max_time_fixed
                  << "limit_num_free_vars_steps:" << limit_num_free_vars_steps
                  << std::endl;
      }
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
    os << "tree depth:" << tree_depth << std::endl;
    os << "depth:" << depth << std::endl;
    os << "heuristic:" << heuristic << std::endl;
    os << "path_from_root:" << path_from_root << std::endl;
    os << "id:" << id << std::endl;
    os << "type:" << int(type) << std::endl;
    os << "node in LGP:" << std::endl;
    node->write(os);
    os << "\nfixed values: " << std::endl;
    for (auto &[k, v] : values) {
      os << k << ":" << v << std::endl;
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
        Qassert(node->decision);
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

    if (subgraph_matched != -1 && type == TYPE::with_values)
      color = "darkgrey";

    if (terminal)
      color = type == TYPE::with_values ? "lightgreen" : "lightblue";

    if (highlight)
      color = "yellow";

    if (!feasible)
      color = "red";

    if (heuristic == BIG_NUMBER_COST_TO_GO && type == TYPE::with_values)
      color = "lightpink3";

    if (heuristic == BIG_NUMBER_COST_TO_GO && type == TYPE::with_variables)
      color = "lightpink";

    std::string out = "[shape=" + shape + " style=filled fillcolor=" + color +
                      " label=\"" + label + "\" ]";
    return out;
  };
};

void tree_to_graphviz_only_computed(std::ostream &os, Compu_node *root) {

  // MISSING: solve this!!
  // I need the full path as label.
  // label = "a:(I B1 B1_ref R)\n(    )\n(   ) h:3,d:1,td:1,ne:1,Tf:0,exp:1 "
  // digraph D {
  // 0 [shape=box style=filled fillcolor=lightgray label="
  // h:4,d:0,td:0,ne:1,Tf:0,exp:0 " ] 1 [shape=box style=filled fillcolor=white
  // label="a:(I B1 B1_ref R) h:3,d:1,td:1,ne:1,Tf:0,exp:1 " ] 2 [shape=box
  // style=filled fillcolor=white label="a:(I B2 B2_ref R)
  // h:3,d:1,td:1,ne:1,Tf:0,exp:2 " ] 3 [shape=box style=filled fillcolor=white
  // label="a:(I B3 B3_ref R) h:5,d:1,td:1,ne:0,Tf:0,exp:" ] 4 [shape=box
  // style=filled fillcolor=white label="a:(I B4 B4_ref R)
  // h:5,d:1,td:1,ne:0,Tf:0,exp:" ] 5 [shape=box style=filled fillcolor=white
  // label="a:(I B5 B5_ref R) h:5,d:1,td:1,ne:0,Tf:0,exp:" ] 6 [shape=box
  // style=filled fillcolor=white label="a:(I B6 B6_ref R)
  // h:5,d:1,td:1,ne:0,Tf:0,exp:" ] 7 [shape=box style=filled
  // fillcolor=lightgray label=" h:3,d:1,td:1,ne:1,Tf:1,exp:3 " ] 8 [shape=box
  // style=filled fillcolor=white label="a:(AonB B1 R block2)
  // h:4,d:2,td:2,ne:0,Tf:0,exp:" ] 9 [shape=box style=filled fillcolor=white
  // label="a:(AonB B1 R block3) h:4,d:2,td:2,ne:0,Tf:0,exp:" ] 10 [shape=box
  // style=filled fillcolor=white label="a:(AonB B1 R block4)
  // h:4,d:2,td:2,ne:0,Tf:0,exp:" ] 11 [shape=box style=filled fillcolor=white
  // label="a:(AonB B1 R block5) h:4,d:2,td:2,ne:0,Tf:0,exp:" ] 12 [shape=box
  // style=filled fillcolor=white label="a:(AonB B1 R block6)
  // h:4,d:2,td:2,ne:0,Tf:0,exp:" ] 13 [shape=box style=filled fillcolor=white
  // label="a:(AonT B1 R g_table) h:2,d:2,td:2,ne:1,Tf:0,exp:4 " ] 14 [shape=box
  // style=filled fillcolor=lightgray label=" h:3,d:1,td:1,ne:1,Tf:1,exp:5 " ]
  // 15 [shape=box style=filled fillcolor=white label="a:(AonB B2 R block1)
  // h:4,d:2,td:2,ne:0,Tf:0,exp:" ] 16 [shape=box style=filled fillcolor=white
  // label="a:(AonB B2 R block3) h:4,d:2,td:2,ne:0,Tf:0,exp:" ] 17 [shape=box
  // style=filled fillcolor=white label="a:(AonB B2 R block4)
  // h:4,d:2,td:2,ne:0,Tf:0,exp:" ] 18 [shape=box style=filled fillcolor=white
  // label="a:(AonB B2 R block5) h:4,d:2,td:2,ne:0,Tf:0,exp:" ] 19 [shape=box
  // style=filled fillcolor=white label="a:(AonB B2 R block6)
  // h:4,d:2,td:2,ne:0,Tf:0,exp:" ] 20 [shape=box style=filled fillcolor=yellow
  // label="a:(AonT B2 R g_table) h:2,d:2,td:2,ne:0,Tf:0,exp:" ] 21 [shape=box
  // style=filled fillcolor=white label="a:(AonB B1 R block2)
  // h:4,d:2,td:2,ne:0,Tf:1,exp:" ] 22 [shape=box style=filled fillcolor=white
  // label="a:(AonB B1 R block3) h:4,d:2,td:2,ne:0,Tf:1,exp:" ] 23 [shape=box
  // style=filled fillcolor=white label="a:(AonB B1 R block4)
  // h:4,d:2,td:2,ne:0,Tf:1,exp:" ] 24 [shape=box style=filled fillcolor=white
  // label="a:(AonB B1 R block5) h:4,d:2,td:2,ne:0,Tf:1,exp:" ] 25 [shape=box
  // style=filled fillcolor=white label="a:(AonB B1 R block6)
  // h:4,d:2,td:2,ne:0,Tf:1,exp:" ] 26 [shape=box style=filled fillcolor=white
  // label="a:(AonT B1 R g_table) h:2,d:2,td:2,ne:0,Tf:1,exp:" ] 27 [shape=box
  // style=filled fillcolor=lightgray label=" h:2,d:2,td:2,ne:0,Tf:2,exp:" ] 28
  // [shape=box style=filled fillcolor=white label="a:(I B1 g_T R)
  // h:3,d:3,td:3,ne:0,Tf:0,exp:" ] 29 [shape=box style=filled fillcolor=white
  // label="a:(I B2 B2_ref R) h:1,d:3,td:3,ne:0,Tf:0,exp:" ] 30 [shape=box
  // style=filled fillcolor=white label="a:(I B3 B3_ref R)
  // h:3,d:3,td:3,ne:0,Tf:0,exp:" ] 31 [shape=box style=filled fillcolor=white
  // label="a:(I B4 B4_ref R) h:3,d:3,td:3,ne:0,Tf:0,exp:" ] 32 [shape=box
  // style=filled fillcolor=white label="a:(I B5 B5_ref R)
  // h:3,d:3,td:3,ne:0,Tf:0,exp:" ] 33 [shape=box style=filled fillcolor=white
  // label="a:(I B6 B6_ref R) h:3,d:3,td:3,ne:0,Tf:0,exp:" ] 34 [shape=box
  // style=filled fillcolor=white label="a:(AonB B2 R block1)
  // h:4,d:2,td:2,ne:0,Tf:1,exp:" ] 35 [shape=box style=filled fillcolor=white
  // label="a:(AonB B2 R block3) h:4,d:2,td:2,ne:0,Tf:1,exp:" ] 36 [shape=box
  // style=filled fillcolor=white label="a:(AonB B2 R block4)
  // h:4,d:2,td:2,ne:0,Tf:1,exp:" ] 37 [shape=box style=filled fillcolor=white
  // label="a:(AonB B2 R block5) h:4,d:2,td:2,ne:0,Tf:1,exp:" ] 38 [shape=box
  // style=filled fillcolor=white label="a:(AonB B2 R block6)
  // h:4,d:2,td:2,ne:0,Tf:1,exp:" ] 39 [shape=box style=filled fillcolor=white
  // label="a:(AonT B2 R g_table) h:2,d:2,td:2,ne:0,Tf:1,exp:" ] 0 -> 1 0 -> 2
  // 0 -> 3
  // 0 -> 4
  // 0 -> 5
  // 0 -> 6
  // 1 -> 7
  // 1 -> 8
  // 1 -> 9
  // 1 -> 10
  // 1 -> 11
  // 1 -> 12
  // 1 -> 13
  // 2 -> 14
  // 2 -> 15
  // 2 -> 16
  // 2 -> 17
  // 2 -> 18
  // 2 -> 19
  // 2 -> 20
  // 7 -> 21
  // 7 -> 22
  // 7 -> 23
  // 7 -> 24
  // 7 -> 25
  // 7 -> 26
  // 13 -> 27
  // 13 -> 28
  // 13 -> 29
  // 13 -> 30
  // 13 -> 31
  // 13 -> 32
  // 13 -> 33
  // 14 -> 34
  // 14 -> 35
  // 14 -> 36
  // 14 -> 37
  // 14 -> 38
  // 14 -> 39
  // }
}

void tree_to_graphviz(std::ostream &os, Compu_node *root,
                      bool only_expanded = false) {

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
      if (!only_expanded || c->num_expands || c->type == TYPE::with_values) {
        global_id++;
        queue.push(std::make_pair(global_id, c));
        childs.push_back(global_id);
      }
    }
    childs_map.insert({id, childs});
  }

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

enum class HEURISTIC {
  BEST_PATH_WITHOUT_LEVEL,
  BEST_PATH,
  COMPU_COST,
  BEST_PATH_COMPU_TIE,
};

enum class Exit_search {
  max_it,
  max_sol,
  max_time,
};

struct Compu_tree {

  bool visualize_solution = false;
  int max_depth_discrete = 10;
  int max_it = 1000;
  int max_num_goals = 10;
  bool interactive = false;
  bool open_viewer_first_time = true;
  bool block_found_plans = false;
  std::vector<std::unique_ptr<KOMO>> komo_problems_solution;
  HEURISTIC heuristic = HEURISTIC::BEST_PATH;

  Compu_node *root;
  std::vector<Compu_node *> open_nodes; // queue for search
  std::vector<Compu_node *> goals;

  Compu_node *pop_and_erase_best() {

    std::cout << "printing nodes " << std::endl;
    for (auto &n : open_nodes) {
      n->write(std::cout);
    }

    auto score_with_expands = [](const auto &a) {
      return a->heuristic + a->depth + a->num_expands;
    };

    std::function<bool(Compu_node *, Compu_node *)> funx;

    switch (heuristic) {
    case HEURISTIC::BEST_PATH_WITHOUT_LEVEL:
      funx = [](const auto &a, const auto &b) {
        return a->heuristic + a->depth < b->heuristic + b->depth;
      };
      break;
    case HEURISTIC::BEST_PATH:
      funx = [&](const auto &a, const auto &b) {
        return score_with_expands(a) < score_with_expands(b);
      };
      break;
    case HEURISTIC::COMPU_COST:
      funx = [&](const auto &a, const auto &b) {
        return score_with_expands(a) - a->max_time_fixed <
               score_with_expands(b) - b->max_time_fixed;
      };
      break;
    case HEURISTIC::BEST_PATH_COMPU_TIE:
      funx = [&](const auto &a, const auto &b) {
        auto sa = score_with_expands(a);
        auto sb = score_with_expands(b);
        if (sa == sb) {
          return -a->max_time_fixed < -b->max_time_fixed;
        } else {
          return sa < sb;
        }
      };
      break;
    default:
      throw std::runtime_error("use valid heuristic");
    }

    Qassert(open_nodes.size());

    // fifo best first queue
    auto it = std::min_element(open_nodes.begin(), open_nodes.end(), funx);

    auto min = *it;
    open_nodes.erase(it);
    return min;
  }

  Exit_search search() {

    int it = 0;
    int expand_counter = 0;

    Exit_search out;

    while (true) {

      if (it++ > max_it) {
        out = Exit_search::max_it;
        break;
      }
      if (goals.size() >= max_num_goals) {
        out = Exit_search::max_sol;
        break;
      }
      // TODO: comput time

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

        node->highlight = true;
        std::ofstream file_("g_.dot");
        tree_to_graphviz(file_, root, true);
        file_.close();
        system_s("dot -Tpdf g_.dot -o g_.pdf");
        if (open_viewer_first_time && it == 1)
          system_s("zathura g_.pdf & ");
        std::cin.get();
        node->highlight = false;
      }

      std::cout << "Chosen node is" << std::endl;
      node->write(std::cout);
      node->expands.push_back(expand_counter++);

      bool recompute_heuristic = false;

      int symbolic_depth = node->node->getTreePath().N - 1;

      if (symbolic_depth > max_depth_discrete) {
        std::cout
            << "skipping a node expansion because symbolic depth is too big"
            << std::endl;
        std::cout << "symbolic_depth " << symbolic_depth
                  << "max_depth_discrete " << max_depth_discrete << std::endl;
        continue;
      }

      auto new_nodes = node->expand(&recompute_heuristic);

      if (node->type == TYPE::with_variables && node->feasible) {
        std::cout << "node has variables, so I add it again" << std::endl;
        open_nodes.push_back(node);
      }

      for (auto &n : new_nodes) {
        if (n->terminal && n->type == TYPE::with_values && n->feasible) {
          goals.push_back(n);
          if (visualize_solution) {
            std::cout << "displaying a solution! " << std::endl;
            n->node->problem(rai::BD_seq).komo->view_play(true);
          }
          komo_problems_solution.emplace_back(std::make_unique<KOMO>());
          komo_problems_solution.back()->clone(
              *n->node->problem(rai::BD_seq).komo, false);

          // KOMO komo_out;
          // std::unique_ptr<KOMO> komo_solution = std::make_unique<KOMO>();
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
        } else {
          open_nodes.push_back(n);
        }
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

        node->highlight = true;
        std::ofstream file_("g_.dot");
        tree_to_graphviz(file_, root, true);
        file_.close();
        system_s("dot -Tpdf g_.dot -o g_.pdf");
        std::cin.get();
        node->highlight = false;
      }
    }

    std::cout << "Termination Criteria " << std::endl;
    std::cout << magic_enum::enum_name(out) << std::endl;
    return out;
  }
};
