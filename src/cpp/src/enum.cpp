/*
 * enum.cpp
 *
 *  Created on: 29-sep-2015
 *      Author: M. El-Kebir
 */

#include <lemon/arg_parser.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "utils.h"
#include "rootedcladisticancestrygraph.h"
#include "rootedcladisticnoisyancestrygraph.h"
#include "rootedcladisticenumeration.h"
#include "rootedcladisticnoisyenumeration.h"
#include "character.h"
#include "charactermatrix.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

typedef RootedCladisticAncestryGraph::StateTreeVector StateTreeVector;

std::map<std::string, StlDoubleVector> read_proportions(std::istream& inFileProportions)
{
  std::map<std::string, StlDoubleVector> props;
  std::string line;
  
  // skip first line
  getline(inFileProportions, line);
  
  StringVector s;
  
  while (inFileProportions.good())
  {
    getline(inFileProportions, line);
    
    boost::split(s, line, boost::is_any_of(","));
    if (s.size() == 0 || line.empty())
    {
      break;
    }
    
    props[s[0]] = StlDoubleVector(s.size() - 1);
    
    for (size_t i = 1; i < s.size() ; ++i)
    {
      props[s[0]][i - 1] = boost::lexical_cast<double>(s[i]);
    }
  }
  
  return props;
}

std::map<std::string, std::string> read_tree(std::istream& inFileTree)
{
  std::map<std::string, std::string> parentMap;
  
  std::string line;
  StringVector s;
  
  while (inFileTree.good())
  {
    getline(inFileTree, line);
    
    boost::split(s, line, boost::is_any_of(","));
    if (s.size() == 0 || line.empty())
    {
      break;
    }
    
    parentMap[s[1]] = s[0];
  }
  
  return parentMap;
}

std::string identify_root(const std::map<std::string, std::string>& parentMap)
{
  std::set<std::string> sourceSet;
  std::set<std::string> targetSet;
  
  for (const auto& kv : parentMap)
  {
    targetSet.insert(kv.first);
    sourceSet.insert(kv.second);
  }
  
  std::set<std::string> res;
  
  std::set_difference(sourceSet.begin(), sourceSet.end(),
                      targetSet.begin(), targetSet.end(),
                      std::inserter(res, res.begin()));
  
  assert(res.size() == 1);
  
  return *res.begin();
}

std::map<std::string, std::set<std::string>> get_parent_to_children(const std::map<std::string, std::string>& parentMap)
{
  std::map<std::string, std::set<std::string>> parentToChildrenMap;
  
  std::set<std::string> nodeSet;
  for (const auto& kv : parentMap)
  {
    nodeSet.insert(kv.first);
    nodeSet.insert(kv.second);
  }
  
  for (const std::string& parent : nodeSet)
  {
    parentToChildrenMap[parent];
    for (const auto& kv : parentMap)
    {
      if (kv.second == parent)
      {
        parentToChildrenMap[parent].insert(kv.first);
      }
    }
  }
  
  return parentToChildrenMap;
}

void get_cumulative_proportions(const std::string& node,
                                const std::map<std::string, std::set<std::string>>& parentToChildrenMap,
                                const std::map<std::string, StlDoubleVector>& props,
                                std::map<std::string, StlDoubleVector>& cumProps)
{
  const int m =  props.begin()->second.size();
  
  // initialize cumulative proportions to 0
  cumProps[node] = props.at(node);
  for (const std::string& child : parentToChildrenMap.at(node))
  {
    get_cumulative_proportions(child, parentToChildrenMap, props, cumProps);

    assert(props.at(node).size() == props.at(child).size());
    for (size_t p = 0; p < m; ++p)
    {
      cumProps[node][p] += cumProps.at(child)[p];
    }
  }
}

std::map<std::string, StlDoubleVector> get_cumulative_proportions(const std::map<std::string, std::string>& parentMap,
                                                                  const std::map<std::string, StlDoubleVector>& props)
{
  assert(props.size());

  const int m =  props.begin()->second.size();

  // initialize cumulative proportions to 0
  std::map<std::string, StlDoubleVector> cumProps;
  for (const auto& kv : parentMap)
  {
    cumProps[kv.first] = StlDoubleVector(m, 0);
    cumProps[kv.second] = StlDoubleVector(m, 0);
  }
  
  const std::string root = identify_root(parentMap);
  const std::map<std::string, std::set<std::string>> parentToChildrenMap = get_parent_to_children(parentMap);
  
  get_cumulative_proportions(root, parentToChildrenMap, props, cumProps);
  
  return cumProps;
}

void read_paction_input(const std::string& fileCnaTree,
                        const std::string& fileCnaProportions,
                        const std::string& fileSnvTree,
                        const std::string& fileSnvProportions,
                        std::stringstream& ss)
{
  std::ifstream inFileCnaTree(fileCnaTree.c_str());
  if (!inFileCnaTree.good())
  {
    std::cerr << "Error opening '" << fileCnaTree << "' for reading" << std::endl;
    exit(1);
  }
  
  std::ifstream inFileCnaProportions(fileCnaProportions.c_str());
  if (!inFileCnaProportions.good())
  {
    std::cerr << "Error opening '" << fileCnaProportions << "' for reading" << std::endl;
    exit(1);
  }
  
  std::ifstream inFileSnvTree(fileSnvTree.c_str());
  if (!inFileSnvTree.good())
  {
    std::cerr << "Error opening '" << fileSnvTree << "' for reading" << std::endl;
    exit(1);
  }
  
  std::ifstream inFileSnvProportions(fileSnvProportions.c_str());
  if (!inFileSnvProportions.good())
  {
    std::cerr << "Error opening '" << fileSnvProportions << "' for reading" << std::endl;
    exit(1);
  }
  
  auto cnaTree = read_tree(inFileCnaTree);
  auto cnaProps = read_proportions(inFileCnaProportions);
  auto cnaRoot = identify_root(cnaTree);
//  auto cnaCumProps = get_cumulative_proportions(cnaTree, cnaProps);
  
  // mapping of index to cna node labels and vice versa
  StringVector cnaLabels;
  cnaLabels.push_back(cnaRoot);
  std::map<std::string, int> cnaInvLabels;
  cnaInvLabels[cnaRoot] = 0;
  for (const auto& kv : cnaProps)
  {
    if (kv.first != cnaRoot)
    {
      cnaInvLabels[kv.first] = cnaLabels.size();
      cnaLabels.push_back(kv.first);
    }
  }
  
  auto snvTree = read_tree(inFileSnvTree);
  auto snvProps = read_proportions(inFileSnvProportions);
  auto snvRoot = identify_root(snvTree);
//  auto snvCumProps = get_cumulative_proportions(snvTree, snvProps);
  
  // mapping of index to snv node labels and vice versa
  StringVector snvLabels;
  snvLabels.push_back(snvRoot);
  std::map<std::string, int> snvInvLabels;
  snvInvLabels[snvRoot] = 0;
  for (const auto& kv : snvProps)
  {
    if (kv.first != snvRoot)
    {
      snvInvLabels[kv.first] = snvLabels.size();
      snvLabels.push_back(kv.first);
    }
  }
  
  assert(cnaProps.begin()->second.size() == snvProps.begin()->second.size());
  
  const int k = std::max(cnaTree.size(), snvTree.size()) + 1;
  const int m = cnaProps.begin()->second.size();
  const int n = 2;
  
  ss << k << " #k" << std::endl;
  ss << m << " #m" << std::endl;
  ss << 2 << " #n" << std::endl;
  
  for (int i = 0; i < k; ++i)
  {
    assert(!(i < cnaLabels.size()) || (cnaProps.count(cnaLabels[i]) == 1));
    assert(!(i < snvLabels.size()) || (snvProps.count(snvLabels[i]) == 1));
    
    for (int p = 0; p < m; ++p)
    {
      bool first = true;
      for (int c = 0; c < n; ++c)
      {
        if (first)
          first = false;
        else
          ss << " ";
        
        if (c == 1)
        {
          if (i < cnaLabels.size())
          {
            ss << cnaProps[cnaLabels[i]][p];
          }
          else
          {
            ss << 0;
          }
        }
        else
        {
          assert(c == 0);
          if (i < snvLabels.size())
          {
            ss << snvProps[snvLabels[i]][p];
          }
          else
          {
            ss << 0;
          }
        }
      }
      ss << std::endl;
    }
  }
  
  ss << std::endl;
  for (int p = 0; p < m; ++p)
  {
    if (p > 0)
      ss << " ";
    ss << "sample_" << p;
  }
  ss << std::endl;
  
  ss << "SNV CNA" << std::endl;
  
  ss << std::endl;
  
  // write state trees
  int count = 1;
  ss << -1;
  for (const std::string& snv_label : snvLabels)
  {
    if (snv_label != snvRoot)
    {
      ++count;
      const std::string& parent = snvTree.at(snv_label);
      ss << " " << snvInvLabels.at(parent);
    }
  }
  for (; count < k; ++count)
  {
    ss << " -2";
  }
  ss << std::endl;
  ss << snvRoot;
  for (const std::string& snv_label : snvLabels)
  {
    if (snv_label != snvRoot)
    {
      ss << " " << snv_label;
    }
  }
  ss << std::endl;
  
  count = 1;
  ss << -1;
  for (const std::string& cna_label : cnaLabels)
  {
    if (cna_label != cnaRoot)
    {
      ++count;
      const std::string& parent = cnaTree.at(cna_label);
      ss << " " << cnaInvLabels.at(parent);
    }
  }
  for (; count < k; ++count)
  {
    ss << " -2";
  }
  ss << std::endl;
  ss << cnaRoot;
  for (const std::string& cna_label : cnaLabels)
  {
    if (cna_label != cnaRoot)
    {
      ss << " " << cna_label;
    }
  }
  ss << std::endl;
  
//  std::cout << ss.str();
}

void read_cladistic_tensor(std::istream& inFile,
                           RealTensor& F,
                           StateTreeVector& S)
{
  g_lineNumber = 0;
  
  inFile >> F;
  F.setLabels(inFile);
  std::string line;
  getline(inFile, line);
  
  S = StateTreeVector(F.n(), StateTree(F.k()));
  
  for (int c = 0; c < F.n(); ++c)
  {
    try
    {
      inFile >> S[c];
    }
    catch (std::runtime_error& e)
    {
      std::cerr << "State tree file. " << getLineNumber() << e.what() << std::endl;
      exit(1);
    }
  }
}

void read_cladistic_tensor_noisy(std::ifstream& inFile,
                                 RealTensor& F,
                                 StateTreeVector& S,
                                 RealTensor& F_lb,
                                 RealTensor& F_ub)
{
  read_cladistic_tensor(inFile, F, S);
  
  std::string line;
  getline(inFile, line);
  inFile >> F_lb;
  getline(inFile, line);
  inFile >> F_ub;
}

void enumerate(int limit,
               int timeLimit,
               int threads,
               int state_tree_limit,
               std::istream& inFile,
               const std::string& intervalFile,
               int lowerbound,
               bool monoclonal,
               const std::string& purityString,
               bool writeCliqueFile,
               bool readCliqueFile,
               const std::string& cliqueFile,
               int offset,
               double delta,
               const IntSet& whiteList,
               SolutionSet& sols)
{
  RealTensor F;
  StateTreeVector S;
  
  try
  {
    read_cladistic_tensor(inFile, F, S);
  }
  catch (std::runtime_error& e)
  {
    std::cerr << "Input file. " << e.what() << std::endl;
    exit(1);
  }
  
  if (delta == 0)
  {
    RootedCladisticAncestryGraph G(F, S);
    G.init();
    G.setLabels(F);
//    G.writeDOT(std::cout);

    RootedCladisticEnumeration alg(G, limit, timeLimit, threads, 0, false, false, whiteList);

    alg.run();

    alg.populateSolutionSet(sols);
  }
  else
  {
    RealTensor F_lb(F);
    RealTensor F_ub(F);

    const int k = F.k();
    const int m = F.m();
    const int n = F.n();
    
    for (int i = 0; i < k; ++i)
    {
      for (int p = 0; p < m; ++p)
      {
        for (int c = 0; c < n; ++c)
        {
          F_lb.set(i, p, c, std::max(0., F(i, p, c) - delta));
          F_ub.set(i, p, c, std::min(1., F(i, p, c) + delta));
        }
      }
    }
    
    RootedCladisticNoisyAncestryGraph G(F, S, F_lb, F_ub);
    G.init();
    G.setLabels(F);

    RootedCladisticNoisyEnumeration alg(G, limit, timeLimit, threads, 0, false, false, whiteList);

    alg.run();

    alg.populateSolutionSet(sols);
  }
    
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cerr << "Generated " << sols.solutionCount() << " solutions" << std::endl;
  }
}

void writeSolTree(const Solution& sol, std::ostream& out)
{
  const int n = sol.A().n();
  PerfectPhyloTree tree(sol.A(), sol.S());
  for (Digraph::ArcIt a(tree.T()); a != lemon::INVALID; ++a)
  {
    const IntPair& ci = tree.nodeToCharState(tree.T().source(a));
    const IntPair& dj = tree.nodeToCharState(tree.T().target(a));
    StlIntVector charStateVector_ci = tree.charStateVector(tree.T().source(a));
    StlIntVector charStateVector_dj = tree.charStateVector(tree.T().target(a));
  
    for (int c = 0; c < n; ++c)
    {
      if (c == 0)
      {
        out << "(";
      }
      else
      {
        out << ", ";
      }
  
      out << sol.S()[c].label(charStateVector_ci[c]);
    }
    out << ")";
    out << "\t";
    for (int c = 0; c < n; ++c)
    {
      if (c == 0)
      {
        out << "(";
      }
      else
      {
        out << ", ";
      }
  
      out << sol.S()[c].label(charStateVector_dj[c]);
    }
    out << ")";
    out << std::endl;
  }
}

void writeSolProps(const Solution& sol, std::ostream& out)
{
  const auto& A = sol.A();
  const auto& U = sol.U();
  PerfectPhyloTree tree(A, sol.S());
//  std::cout << A << std::endl;
//  std::cout << U << std::endl;

  const int n = sol.observedF().n();
  const int m = sol.observedF().m();

  out << "clone\tsnv\tcna";
  for (int p = 0; p < m; ++p)
  {
    out << "\t" << sol.observedF().getRowLabel(p);
  }
  out << std::endl;
  
  int uIdx = 0;
  for (int i = 0; i < A.k(); ++i)
  {
    for (int c = 0; c < A.n(); ++c)
    {
      if (A.defined(c, i))
      {
        Digraph::Node v_ci = tree.charStateToNode(c, i);
        StlIntVector charStateVector = tree.charStateVector(v_ci);

        for (int cc = 0; cc < n; ++cc)
        {
          if (cc == 0)
          {
            out << "(";
          }
          else
          {
            out << ", ";
          }

          out << sol.S()[cc].label(charStateVector[cc]);
        }
        out << ")";
        out << "\t" << sol.S()[0].label(charStateVector[0]);
        out << "\t" << sol.S()[1].label(charStateVector[1]);
        for (int p = 0; p < m; ++p)
        {
          out << "\t" << U(p, uIdx);
        }
        
        out << std::endl;
        uIdx++;
      }
      else if (i != 0)
      {
        uIdx++;
      }
    }
  }
}

int main(int argc, char** argv)
{
  double delta = 0;
  std::stringstream ss;
  
  if (argc >= 5)
  {
    read_paction_input(argv[1], argv[2], argv[3], argv[4], ss);
//    std::cout << ss.str();
//    return 0;
  }
  else
  {
    std::cerr << argv[0] << " <cna_tree> <cna_props> <snv_tree> <snv_props> (<delta> <output_prefix>)" << std::endl;
    return 0;
  }
  
  if (argc > 5)
  {
    delta = boost::lexical_cast<double>(argv[5]);
  }
    
  SolutionSet sols;
  enumerate(-1, -1, 1, -1, ss, "", 0, false, "", false, false, "", 0, delta, IntSet(), sols);
  
  for (int solIdx = 0; solIdx < sols.solutionCount(); ++solIdx)
  {
    const Solution& sol = sols.solution(solIdx);
    
    if (argc <= 6)
    {
      writeSolProps(sol, std::cout);
      writeSolTree(sol, std::cout);
      std::cout << std::endl;
    }
    else
    {
      std::ofstream outProps((std::string(argv[6]) + "_props_" + boost::lexical_cast<std::string>(solIdx) + ".tsv").c_str());
      writeSolProps(sol, outProps);
      outProps.close();
      
      std::ofstream outTree((std::string(argv[6]) + "_tree_" + boost::lexical_cast<std::string>(solIdx) + ".tsv").c_str());
      writeSolTree(sol, outTree);
      outTree.close();
    }
  }
  
  return 0;
}
