#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
#define BOOST_TEST_DYN_LINK

#include "pddlss_gnlpp.hpp"

BOOST_AUTO_TEST_CASE(hello) { std::cout << "hello world " << std::endl; }

BOOST_AUTO_TEST_CASE(store_graph) {
  // BOOST_CHECK(false);
}

BOOST_AUTO_TEST_CASE(easy_minimal_conflict) {}

BOOST_AUTO_TEST_CASE(plan_solve_join) {

  rai::setParameter<bool>("KOMO/mimicStable", false);

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
    C.view(true);

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
    auto out = NLP_Solver().setProblem(factored_nlp).solve();

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
    C.view(true);

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

  {
    std::ofstream file("g.dot");
    tree_to_graphviz(file, &croot);
  }

  {
    std::ofstream file("g_.dot");
    bool only_expanded = true;
    tree_to_graphviz(file, &croot, only_expanded);
  }
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
    C.view(true);

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

  system_s("rm -rf " + CONFLICT_FOLDER + " && mkdir " + CONFLICT_FOLDER);
  system_s("rm -rf " + PATH_FEASIBLE_SAMPLES + " && mkdir " +
           PATH_FEASIBLE_SAMPLES);

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
  auto out = NLP_Solver().setProblem(factored_nlp).solve();

  std::map<std::string, arr> new_values;
  for (auto &v : factored_nlp->__variableIndex) {
    std::cout << v << std::endl;
    new_values[v.nameID.p] = v.getValue(true);
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
    C.view(true);

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
    C.view(true);

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

  rai::String folFile = "fol_lab_bench_two_one_object.g";
  rai::String confFile = "./lab_setting_two_two_objects.g";

  rai::String goal = "(on goal2_table block1) (on r_gripper block2)";
  rai::String plan = "(pick block1 block1_ref r_gripper)\n(pick block2 "
                     "block2_ref l_gripper)";

  bool visualize = false;
  rai::Configuration C;
  C.addFile(confFile);
  if (visualize)
    C.view(true);

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
    C.view(true);

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
      auto out = NLP_Solver().setProblem(factored_nlp2).solve();
    }

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
      auto out = NLP_Solver().setProblem(factored_nlp3).solve();
    }

    factored_nlp3->subSelect({}, {});
    factored_nlp3->report(std::cout, 5);

    // lets do a warmstart
  }

#endif
}

BOOST_AUTO_TEST_CASE(grasp_models) {
  // lets test the quim grasp model...

  rai::setParameter<bool>("KOMO/mimicStable", false);

  rai::String folFile = "fol_lab_bench_two_one_object.g";
  rai::String confFile = "./lab_setting_two_two_objects.g";
  rai::String plan1 = "(pick block1 block1_ref r_gripper)\n"
                      "(place block1 r_gripper goal2_table)";

  // rai::String plan2 = "(pick block1 block1_ref r_gripper)";
  // rai::String plan3 = "(pick block2 block2_ref r_gripper)"
  //                     "(place block2 r_gripper goal2_table)"
  //                     "(pick block1 block1_ref r_gripper)";

  bool visualize = false;
  rai::Configuration C;
  C.addFile(confFile);

  if (visualize)
    C.view(true);

  rai::LGP_Tree lgp(C, folFile);
  auto node = lgp.walkToNode(plan1);
  auto [graph, factored_nlp] = get_nlps_and_graph(node);

  {
    std::cout << "writing clean " << std::endl;
    std::ofstream file("tmp/g.dot");
    write_graphviz_easy(file, graph, make_label_writer3(graph));
    file.close();
    system_s("dot -Tpdf tmp/g.dot -o tmp/g.pdf");
  }

  arr default_x = factored_nlp->komo.x.copy();
  auto &k = factored_nlp->komo;

  for (auto &d : k.pathConfig.activeDofs) {
    if (in_vector(d->joint()->type,
                  {rai::JT_hingeX, rai::JT_hingeY, rai::JT_hingeZ})) {
      d->sampleSdv = .1;
    }
  }

  k.initRandom(3);
  k.view(true);
  k.view_play(true);

  // double sampleSdv=.01; //sdv of gaussian around default

  Solve_nlp_with_counter solver;
  bool feas = solver.solve_nlp(factored_nlp, default_x, graph, visualize, true);
  Qassert(feas);

  // visualize

  std::cout << "reporting" << std::endl;
  factored_nlp->komo.view(true);
  factored_nlp->komo.view_play(true);

  std::cout << "writing clean " << std::endl;
  std::ofstream file("tmp/g.dot");
  write_graphviz_easy(file, graph, make_label_writer3(graph));
  file.close();
  system_s("dot -Tpdf tmp/g.dot -o tmp/g.pdf");
}

BOOST_AUTO_TEST_CASE(small_vs_big_table) {

  std::string timesstamp = get_time_stamp();
  int argc = boost::unit_test::framework::master_test_suite().argc;
  char **argv = boost::unit_test::framework::master_test_suite().argv;

  rai::initCmdLine(argc, argv);
  rai::setParameter<bool>("KOMO/mimicStable", false);

  bool interactive = rai::getParameter<bool>("interactive", false);
  const bool check_plan = rai::getParameter<bool>("check_plan", false);
  const bool solve_gnlp = rai::getParameter<bool>("solve_gnlp", false);
  const bool solve_ssgnlp = rai::getParameter<bool>("solve_ssgnlp", true);

  int visualize = rai::getParameter<int>("vis", 0);
  int seed = rai::getParameter<int>("seed", 0);

  double table_size = rai::getParameter<double>("table_size", -1.);

  rai::String base_path = rai::getParameter<rai::String>(
      "base_path", "small_vs_big_table/" + rai::String(timesstamp));
  bool visualize_solutions = rai::getParameter<bool>("vs", false);
  rai::String folFile =
      rai::getParameter<rai::String>("folFile", "./fol_lab_bench_auto.g");
  rai::String confFile =
      rai::getParameter<rai::String>("confFile", "./models/table_experiment.g");
  rai::String goalFile = rai::getParameter<rai::String>("goalFile", "NONE");
  rai::String goal = rai::getParameter<rai::String>("goal", "NONE");
  rai::String plan = rai::getParameter<rai::String>("plan", "NONE");
  rai::String plan_file = rai::getParameter<rai::String>("plan_file", "NONE");

  if (base_path.getLastN(1) != "/")
    base_path = base_path + "/";

  rai::setParameter<rai::String>("out_file_benchmark",
                                 base_path + "__benchmark.txt");

  std::filesystem::create_directories(base_path.p);

  CONFLICT_FOLDER = base_path.p + CONFLICT_FOLDER;
  PATH_FEASIBLE_SAMPLES = base_path.p + PATH_FEASIBLE_SAMPLES;

  std::cout << "time " << timesstamp << std::endl;
  std::cout << "** parsed parameters:\n" << rai::getParameters() << '\n';

  if (seed >= 0) {
    srand(seed);
    rnd.seed(seed);
  } else {
    srand(time(NULL));
    rnd.seed(time(NULL));
  }

  if (plan == "NONE") {
    if (plan != "NONE") {
      std::ifstream file(plan_file);
      std::string plan_str;
      plan_str = std::string((std::istreambuf_iterator<char>(file)),
                             std::istreambuf_iterator<char>());
      plan = plan_str.c_str();
    }
  }
  std::cout << "parsed plan is " << std::endl << plan << std::endl;

  rai::Configuration C;
  std::cout << "loading config file " << confFile << std::endl;

  // check that files exists

  if (!std::filesystem::exists(folFile.p)) {
    std::cout << "file:" << folFile << " does not exist" << std::endl;
    exit(0);
  }

  if (!std::filesystem::exists(confFile.p)) {
    std::cout << "file:" << confFile << " does not exist" << std::endl;
    exit(0);
  }

  // throw -1;

  C.addFile(confFile);

  // TODO: refactor!!!

  if (table_size > 0) {
    // in the first problem, i can change the shape of the table number 0
    if (auto f = C.getFrame("goal_table", false); f)
      f->setShape(rai::ST_ssBox, {table_size, table_size, .04, .01});

    // in the second problem, i can change the shape of the table number 1
    if (auto f = C.getFrame("goal1_table", false); f)
      f->setShape(rai::ST_ssBox, {table_size, table_size, .04, .01});
  }

  if (visualize)
    C.view(true);

  bool generate_auto = false;
  if (folFile.contains("auto"))
    generate_auto = true;

  rai::LGP_Tree lgp(C, folFile, generate_auto);

  if (goal != "NONE") {
    lgp.fol.addTerminalRule(goal);
  } else if (goalFile != "NONE") {
    std::ifstream file(goalFile);
    std::string goal_str = std::string((std::istreambuf_iterator<char>(file)),
                                       std::istreambuf_iterator<char>());
    std::cout << "GOAL IS " << goal_str << std::endl;
    lgp.fol.addTerminalRule(goal_str.c_str());
  }

  if (check_plan) {
    CHECK(plan != "NONE", "no plan specified");
    std::cout << "plan is " << plan << std::endl;
    auto node = lgp.walkToNode(plan);

    StringAA collision_list = {
        {"f_block1_col"}, {"f_block2_col"},  {"f_block3_col"},
        {"f_block4_col"}, {"f_block5_col"},  {"f_block6_col"},
        {"obs_1"},        {"col_l_gripper"}, {"col_r_gripper"}};
    StringAA collision_list_here = {};

    for (auto &c : collision_list) {
      auto &C = node->tree.kin;
      if (C.getFrame(c(0)))
        collision_list_here.append(c);
    }

    std::cout << "collision list is " << collision_list_here << std::endl;

    auto [graph, factored_nlp] = get_nlps_and_graph(node, collision_list_here);

    {
      std::cout << "writing clean " << std::endl;
      std::string rand_id = gen_random_id();
      std::string filename = "/tmp/pddlss_gnlpp/g" + rand_id + ".dot";
      create_dir_if_necessary(filename.c_str());
      std::cout << "writing file " << filename << std::endl;
      std::ofstream file(filename);
      write_graphviz_easy(file, graph, make_label_writer3(graph));
      file.close();
      system_s(string_format("dot -Tpdf {} -o {}.pdf", filename, filename));
    }

    factored_nlp->subSelect({}, {});
    auto out = NLP_Solver().setProblem(factored_nlp).solve();

    factored_nlp->komo.view_play(true);

    std::cout << *out << std::endl;
    std::cout << "***" << std::endl;

    // factored_nlp->komo.view(true);
    // factored_nlp->komo.view_play(true);

    return;
  }
  if (solve_gnlp) {
    // solve using Graph NLP

    OptConflict opt_conf{
        .nlp_heu = true,
        .graph = true,
        .relax = true,
        .filter = false,
        .linear_filtering = true, // I can use Linear or QuickXplain
        .mode = RELAX_MODE::def,
        .use_gnn = false,
    };

    OptGraphLGP opt{
        .lgp = &lgp,
        .obsFile = rai::getParameter<rai::String>("obsFile", "none"),
        .obsFilePairs = "none",
        .verbosity = 0,
        .visualize = static_cast<bool>(visualize),
        .path_to_downward_folder = PATH_TO_DOWNWARD_FOLDER,
        .path_to_utils = "../rai/utils",
        .id = "qq",
        .opt = opt_conf,
        .compute_full_traj = false,
        .imax = 100,
    };

    Datasd data;
    Datass datas;

    std::pair<bool, rai::String> out = graph_lgp(opt, &data, &datas);

    std::cout << out.first << " " << out.second << std::endl;
    BOOST_CHECK(out.first);
    return;
  }
  if (solve_ssgnlp) {
    // solve using ssgnlp

    auto root = lgp.root;
    std::cout << "lgp node " << std::endl;
    root->write(std::cout);

    VISUALIZE_KOMO = visualize > 1;
    Compu_node croot;
    croot.node = root;
    croot.tree = &lgp;

    system_s("rm -rf " + CONFLICT_FOLDER + " && mkdir " + CONFLICT_FOLDER);
    system_s("rm -rf " + PATH_FEASIBLE_SAMPLES + " && mkdir " +
             PATH_FEASIBLE_SAMPLES);

    croot.compute_heuristic();
    bool recompute_heuristic;
    croot.compute_values(&recompute_heuristic);
    Qassert(!recompute_heuristic); // the initial state should be feasible!

    Compu_tree compu_tree;

    // CONTINUE HERE: GO back to check the conditional sampling operations!!

    compu_tree.visualize_solution = visualize > 0;
    compu_tree.max_depth_discrete =
        rai::getParameter<int>("ct/max_depth_discrete", 10);
    compu_tree.max_it = rai::getParameter<int>("ct/max_it", 100);

    rai::String heu_str =
        rai::getParameter<rai::String>("ct/heuristic", "BEST_PATH_COMPU_TIE");

    compu_tree.heuristic = magic_enum::enum_cast<HEURISTIC>(heu_str.p).value();

    compu_tree.max_num_goals = rai::getParameter<int>("ct/max_num_goals", 1);
    compu_tree.root = &croot;
    compu_tree.open_nodes = {&croot};
    compu_tree.interactive = interactive;
    compu_tree.block_found_plans =
        rai::getParameter<bool>("ct/block_found_plans", false);
    Exit_search status;
    try {
      status = compu_tree.search();
    } catch (std::exception &e) {
      std::cout << "ERROR " << std::endl;
      std::cout << e.what() << std::endl;
      status = Exit_search::error;
    }

    {
      rai::String file_report = base_path + "__report.txt";
      std::cout << "file_report: " << file_report << std::endl;
      std::ofstream report_out(base_path + "__report.txt");
      if (compu_tree.goals.size()) {
        std::cout << "problem solved! " << std::endl;
        report_out << "solved: " << 1 << std::endl;
        using namespace magic_enum::ostream_operators;
        report_out << "exit_status: " << status << std::endl;
      } else {
        report_out << "solved: " << 0 << std::endl;
        using namespace magic_enum::ostream_operators;
        report_out << "exit_status: " << status << std::endl;
      }
    }

    {
      std::ofstream file(base_path + "g.dot");
      tree_to_graphviz(file, &croot);
      file.close();

      std::ofstream file2(base_path + "g_.dot");
      tree_to_graphviz(file2, &croot, true);
      file2.close();

      system_s(
          ("dot -Tpdf " + base_path + "g.dot -o " + base_path + "g.pdf").p);

      system_s(
          ("dot -Tpdf " + base_path + "g_.dot -o " + base_path + "g_.pdf").p);
    }

    BOOST_TEST(compu_tree.goals.size());

    // print the solutions

    auto print_solutions = [&](std::ostream &out) {
      for (auto &g : compu_tree.goals) {
        out << "***\nGOAL NODE***\n" << std::endl;
        g->write(out);
        out << std::endl;
        auto compu = g->get_compu_history();
        out << "**compu history**" << std::endl;
        for (auto &vv : compu) {
          out << "compu decision" << std::endl;
          for (auto &[k, v] : vv) {
            out << k << ": " << v << std::endl;
          }
        }
      }
    };

    print_solutions(std::cout);
    ofstream out_solutions(base_path + "solutions.txt");
    print_solutions(out_solutions);

    if (visualize_solutions) {
      int counter = 0;
      for (auto &komo : compu_tree.komo_problems_solution) {
        std::cout << "diplaying solution " << counter++ << std::endl;
        std::string video_path = std::string(base_path.p) + "solutions/s" +
                                 std::to_string(counter) + "/";

        std::filesystem::create_directories(video_path);
        komo->view_play(false, .05, video_path.c_str());
        komo->view_close();
        int out_code = system_s("cd " + video_path + "&&" +
                                "ffmpeg -framerate 2 -i \"%04d.png\" out.mp4");
        BOOST_CHECK(out_code == 0);
      }

      // for (auto &g : compu_tree.goals) {
      //   std::cout << "diplaying solution " << counter++ << std::endl;
      //   g->node->problem(rai::BoundType::BD_seq).komo->view(true);
      //   g->node->problem(rai::BoundType::BD_seq).komo->view_play(true);
      // }
    }
  }
}

BOOST_AUTO_TEST_CASE(repeat_ops) {
  // TODO: make sure that I don't use exactly the same
  // Graph to warmstart another expansion. Then I should store a kind of list.

  rai::setParameter<bool>("KOMO/mimicStable", false);
  rai::String folFile = "./fol_lab_bench_auto.g";
  rai::String confFile = "./models/table_experiment_minimal.g";

  rai::String goal = "(on goal_table block1)";
  rai::String plan = "(pick block1 block1_ref r_gripper)\n(placeontable block1 "
                     "r_gripper goal_table)\n";

  bool visualize = true;
  rai::Configuration C;
  C.addFile(confFile);
  auto f = C.getFrame("goal_table");
  f->setShape(rai::ST_ssBox, {.2, .2, .04, .01});

  // TODO: I need option to create SOS cost, on the init guess
  // of the placement of objects. (and grasps?) -- I could do the
  // same for top grasp with one parameter :)
  // f->getShape().createMeshes();

  if (visualize)
    C.view(true);

  bool generate_auto = true;
  rai::LGP_Tree lgp(C, folFile, generate_auto);
  lgp.fol.addTerminalRule(goal);

  bool check_plan = false;
  bool solve_gnlp = false;
  bool solve_ssgnlp = true;

  if (check_plan) {
    auto node = lgp.walkToNode(plan);

    // StringAA collision_list = {{"block1"}, {"block2"}, {"block3"},
    //                            {"block4"}, {"block5"}, {"block6"}};

    StringAA collision_list = {};

    auto [graph, factored_nlp] = get_nlps_and_graph(node, collision_list);

    {
      std::cout << "writing clean " << std::endl;
      std::ofstream file("tmp/g.dot");
      write_graphviz_easy(file, graph, make_label_writer3(graph));
      file.close();
      system_s("dot -Tpdf tmp/g.dot -o tmp/g.pdf");
    }

    factored_nlp->subSelect({}, {});
    auto out = NLP_Solver().setProblem(factored_nlp).solve();
    BOOST_CHECK(out->feasible);

    std::cout << *out << std::endl;
    std::cout << "***" << std::endl;

    factored_nlp->komo.view(true);
    factored_nlp->komo.view_play(true);
  }
  if (solve_gnlp) {
    // solve using Graph NLP

    OptConflict opt_conf{
        .nlp_heu = true,
        .graph = true,
        .relax = true,
        .filter = false,
        .linear_filtering = true, // I can use Linear or QuickXplain
        .mode = RELAX_MODE::def,
        .use_gnn = false,
    };

    OptGraphLGP opt{
        .lgp = &lgp,
        .obsFile = "none",
        .obsFilePairs = "none",
        .verbosity = 0,
        .visualize = false,
        .path_to_downward_folder = PATH_TO_DOWNWARD_FOLDER,
        .path_to_utils = "../rai/utils",
        .id = "qq",
        .opt = opt_conf,
        .compute_full_traj = false,
        .imax = 100,
    };

    Datasd data;
    Datass datas;

    std::pair<bool, rai::String> out = graph_lgp(opt, &data, &datas);

    std::cout << out.first << " " << out.second << std::endl;
    BOOST_CHECK(out.first);
  }
  if (solve_ssgnlp) {
    // solve using ssgnlp

    auto root = lgp.root;
    std::cout << "lgp node " << std::endl;
    root->write(std::cout);

    VISUALIZE_KOMO = false;
    Compu_node croot;
    croot.node = root;
    croot.tree = &lgp;

    system_s("rm -rf " + CONFLICT_FOLDER + " && mkdir " + CONFLICT_FOLDER);
    system_s("rm -rf " + PATH_FEASIBLE_SAMPLES + " && mkdir " +
             PATH_FEASIBLE_SAMPLES);

    croot.compute_heuristic();
    bool recompute_heuristic;
    croot.compute_values(&recompute_heuristic);
    Qassert(!recompute_heuristic); // the initial state should be feasible!

    Compu_tree compu_tree;
    compu_tree.max_it = 200;
    compu_tree.max_num_goals = 5;
    compu_tree.root = &croot;
    compu_tree.open_nodes = {&croot};
    compu_tree.interactive = true;
    compu_tree.block_found_plans = false;
    compu_tree.search();
    BOOST_ASSERT(compu_tree.goals.size());

    std::ofstream file("g.dot");
    tree_to_graphviz(file, &croot);

    std::ofstream file2("g_.dot");
    tree_to_graphviz(file2, &croot, true);

    // print the solutions
    for (auto &g : compu_tree.goals) {
      std::cout << "***\nGOAL NODE***\n" << std::endl;
      g->write(std::cout);
      std::cout << std::endl;
      auto compu = g->get_compu_history();
      std::cout << "**compu history**" << std::endl;
      for (auto &vv : compu) {
        std::cout << "compu decision" << std::endl;
        for (auto &[k, v] : vv) {
          std::cout << k << ": " << v << std::endl;
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(rrt) {

  auto filename = "./models/check_rrt.g";
  rai::Configuration C;
  C.addFile(filename);

  ConfigurationProblem P(C);
  const arr starts = C.getJointState();

  const arr goals = {1.5, -.5, -.0, -2, -0, 2, -.5};
  double _stepsize = .2;
  int verbose = 3;
  RRT_PathFinder pathfinder(P, starts, goals, _stepsize, verbose);

  pathfinder.planConnect();

  auto out = make_shared<PathResult>(pathfinder.path);
  std::cout << out->path << std::endl;

  std::cout << *out << std::endl;
}

BOOST_AUTO_TEST_CASE(two_robot_move) {

  // check if the intermediate position is constrained or not, based on the
  // the table size.
  // How to warmstart each variable efficiently?
}

BOOST_AUTO_TEST_CASE(new_warmstart) {

  // check the warmstart of MARC.
}
#if 0
BOOST_AUTO_TEST_CASE(two_robots) {


  int argc = boost::unit_test::framework::master_test_suite().argc;
  char **argv = boost::unit_test::framework::master_test_suite().argv;

  std::string timesstamp = get_time_stamp();
  int visualize = 0;
  bool interactive = false;
  double table_size = .2;
  rai::String base_path = "exp_two_robots/" + rai::String(timesstamp) + "/";

  rai::initCmdLine(argc, argv);
  rai::setParameter<bool>("KOMO/mimicStable", false);

  visualize = rai::getParameter<int>("vis", visualize);
  interactive = rai::getParameter<bool>("interactive", interactive);
  table_size = rai::getParameter<double>("table_size", table_size);
  int seed = rai::getParameter<int>("seed", 0);
  base_path = rai::getParameter<rai::String>(
      "base_path", "small_vs_big_table/" + rai::String(timesstamp));


  rai::String folFile = "fol_lab_bench_two_one_object.g";
  rai::String confFile = "./lab_setting_two_two_objects.g";

  if (seed >= 0)
    srand(seed);
  else {
    srand(time(NULL));
  }

  rai::Configuration C;

  rai::String goal = "(on goal1_table block1) (on block1 block2)";
  rai::String plan = "(pick block1 block1_ref r_gripper)\n(place block1 "
                     "r_gripper goal1_table)\n(pick block2 block2_ref "
                     "r_gripper)\n(place block2 r_gripper block1)";


  // try to solve with nonlinear optimization 
  C.addFile(confFile);
  if (visualize)
    C.view(true);

  rai::LGP_Tree lgp(C, folFile);
  lgp.fol.addTerminalRule(goal);








}

#endif
