#include <iostream>
#include <vector>
#include <climits>
#include <chrono>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <random>
using namespace std;

struct Class;
struct Station;
struct Parent;

using VI = vector<int>;
using VC = vector<Class>;
using VS = vector<Station>;
using MI = vector<VI>;
using MP = vector<Parent>;
using is = ifstream;
using os = ofstream;

int UNDEF = -1;

struct Class
{
  int id, ncars;
  VI improvements;
};

struct Station
{
  int requirements;
  VI line;
};

struct Parent
{
  VI solution;
  int penalization;
  double fitness;
};

double now()
{
  return clock() / double(CLOCKS_PER_SEC);
}

os open(const string &out)
{
  os f;
  f.open(out, os::out | os::trunc);
  f.setf(ios::fixed);
  f.precision(1);
  return f;
}

void write_solution(const VI &cs, int cp, const double &elapsed_time, const string &out)
{
  os f = open(out);
  f << cp << " " << elapsed_time << endl;

  f << cs[0];
  for (int i = 1; i < cs.size(); ++i)
    f << " " << cs[i];

  f << endl;

  f.close();
}

/* Returns a vector containing for each index the number of cars that are left for that class */
VI count_cleft(const VC &classes)
{
  int K = classes.size();
  VI cleft(K);
  for (int i = 0; i < K; ++i)
    cleft[i] = classes[i].ncars;
  return cleft;
}

int count_penalization(const VI &line, int ce, int ne)
{

  int total_penalizations, requirements;
  total_penalizations = requirements = 0;

  // count inicial windows
  int uw = 0;
  for (int i = ne; i > 0; --i)
  {
    requirements += line[uw];
    total_penalizations += max(requirements - ce, 0);
    ++uw;
  }

  // count mid penalizations
  for (; uw < line.size(); ++uw)
  {
    requirements += line[uw];
    requirements -= line[uw - ne];
    total_penalizations += max(requirements - ce, 0);
  }

  // count end windows
  while (requirements > 0)
  {
    requirements -= line[uw - ne];
    total_penalizations += max(requirements - ce, 0);
    ++uw;
  }

  return total_penalizations;
}

MI createSM(const VI &P, int C, int M, const VC &classes)
{

  MI SM(M, VI(C));
  for (int i = 0; i < C; ++i)
    for (int j = 0; j < M; ++j)
      SM[j][i] = classes[P[i]].improvements[j];
  return SM;
}

int penalization(const VI &P, const VI &ce, const VI &ne, const VC &classes)
{
  int C = P.size();
  int M = ce.size();
  MI stations_matrix = createSM(P, C, M, classes);
  int total_pen = 0;

  for (int i = 0; i < M; ++i)
    total_pen += count_penalization(stations_matrix[i], ce[i], ne[i]);
  return total_pen;
}

void generate_permutation(VI &parent)
{
  random_shuffle(parent.begin(), parent.end());
}

VI create_parent(int C, const VI &cleft)
{
  VI parent;
  parent.reserve(C);
  for (int i = 0; i < cleft.size(); ++i)
  {
    int cl = cleft[i];
    for (; cl > 0; --cl)
      parent.push_back(i);
  }
  return parent;
}

MP generate_parents(int C, const VI &cleft, const VI &ce, const VI &ne, const VC &classes)
{
  VI P = create_parent(C, cleft);

  // select the population size (it will remain constant during the iterations)
  int psize = 600;

  MP parents(psize);

  for (int i = 0; i < psize; ++i)
  {
    generate_permutation(P);
    parents[i] = Parent{P, penalization(P, ce, ne, classes), 0};
  }
  return parents;
}

bool improved(Parent &best_individual, const Parent &P)
{
  if (P.fitness > best_individual.fitness)
  {
    best_individual = P;
    return true;
  }
  return false;
}

bool best_fit(const Parent &P1, const Parent &P2)
{
  return P1.fitness > P2.fitness;
}

void sort_parents(MP &parents)
{
  sort(parents.begin(), parents.end(), best_fit);
}

VI generate_random_cut_points(int C)
{
  random_device rd;                            // obtain a random number from hardware
  mt19937 gen(rd());                           // seed the generator
  uniform_int_distribution<> distr1(0, C - 2); // define the range
  int cut1 = distr1(gen);
  uniform_int_distribution<> distr2(cut1 + 1, C - 1);
  int cut2 = distr2(gen);
  return {cut1, cut2};
}

int most_cleft_index(const VI &cleft)
{
  auto it = max_element(cleft.begin(), cleft.end());
  return distance(cleft.begin(), it);
}

void adapt_solution(VI &child, VI &cleft, int left, int right)
{
  int C = child.size();
  int j = 0;

  while (j < C)
  {
    if (j == left)
    {
      j = right + 1;
      if (j == C)
        break;
    }
    if (int a = child[j]; cleft[a] > 0)
      --cleft[a];
    else
    {
      child[j] = most_cleft_index(cleft);
      --cleft[child[j]];
    }
    ++j;
  }
}

void adapt_solution2(VI &child, VI &cleft, int left, int right)
{
  int C = child.size();

  for (int j = left; j <= right; ++j)
  {
    if (int a = child[j]; cleft[a] > 0)
      --cleft[a];
    else
    {
      child[j] = most_cleft_index(cleft);
      --cleft[child[j]];
    }
  }
}

VI uniform(const VI &P1, const VI &P2, VI cleft)
{
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> distr(1, 2);

  VI child(P1.size());

  for (int i = 0; i < P1.size(); ++i)
  {
    if (distr(gen) == 0)
    {
      int a = P1[i];
      if (cleft[a] > 0)
      {
        child[i] = a;
        --cleft[a];
      }
      else if (cleft[P2[i]] > 0)
      {
        child[i] = P2[i];
        --cleft[P2[i]];
      }
      else
      {
        int b = most_cleft_index(cleft);
        child[i] = b;
        --cleft[b];
      }
    }
    else
    {
      int a = P2[i];
      if (cleft[a] > 0)
      {
        child[i] = a;
        --cleft[a];
      }
      else if (cleft[P1[i]] > 0)
      {
        child[i] = P1[i];
        --cleft[P1[i]];
      }
      else
      {
        int b = most_cleft_index(cleft);
        child[i] = b;
        --cleft[b];
      }
    }
  }
  return child;
}

MI recombinate(const VI &P1, const VI &P2, int left, int right, const VI &cleft)
{
  int C = P1.size();

  VI children1, children2, children3, children4;
  children1 = children3 = P1;
  children2 = children4 = P2;
  VI cleft1, cleft2, cleft3, cleft4;
  cleft1 = cleft2 = cleft3 = cleft4 = cleft;

  // fisrt pair of children
  for (int i = left; i <= right; ++i)
  {
    children1[i] = P2[i];
    --cleft1[P2[i]];
    children2[i] = P1[i];
    --cleft2[P1[i]];
  }

  // second pair of children
  for (int i = 0; i < left; ++i)
  {
    children3[i] = P2[i];
    --cleft3[P2[i]];
    children4[i] = P1[i];
    --cleft4[P1[i]];
  }
  for (int i = right + 1; i < C; ++i)
  {
    children3[i] = P2[i];
    --cleft3[P2[i]];
    children4[i] = P1[i];
    --cleft4[P1[i]];
  }

  adapt_solution(children1, cleft1, left, right);
  adapt_solution(children2, cleft2, left, right);
  adapt_solution2(children3, cleft3, left, right);
  adapt_solution2(children4, cleft4, left, right);

  VI children5 = uniform(P1, P2, cleft);

  return {children1, children2, children3, children4, children5};
}

MP recombination(const VI &P1, const VI &P2, const VI &cleft, const VI &ce, const VI &ne, const VC &classes)
{
  int C = P1.size();

  VI cuts = generate_random_cut_points(C);
  MI recombination = recombinate(P1, P2, cuts[0], cuts[1], cleft);
  MP children;
  children.reserve(5);

  for (int i = 0; i < 5; ++i)
    children.push_back(Parent{recombination[i], penalization(recombination[i], ce, ne, classes), 0});

  return children;
}

void mutation(MP &children)
{
  // add a mutation to the children with propability 0.01

  if (rand() / (double)RAND_MAX < 0.01)
  {

    int C = children[0].solution.size();
    int j = rand() % C;
    int k = rand() % C;

    for (int i = 0; i < children.size(); ++i)
      swap(children[i].solution[j], children[i].solution[k]);
  }
}

void update_fitness(MP &parents)
{
  int n = parents.size();
  double total_sum = 0;
  for (const Parent &parent : parents)
    total_sum += parent.penalization;

  for (Parent &parent : parents)
    parent.fitness = (double)(total_sum - parent.penalization) / (total_sum * (n - 1));
}

MP roulette_wheel_selection(const MP &parents)
{
  double offset = 0.0;
  VI picks;
  picks.reserve(2);

  while (picks.size() != 2)
  {
    double rndNumber = rand() / (double)RAND_MAX;
    for (int i = 0; i < parents.size(); i++)
    {
      offset += parents[i].fitness;
      if (rndNumber < offset)
      {
        if (picks.size() == 0 or picks[0] != i)
          picks.push_back(i);
        else
        {
          int j = (i - 1 >= 0 ? i - 1 : i + 1);
          picks.push_back(j);
        }
        break;
      }
    }
    offset = 0.0;
  }

  return {parents[picks[0]], parents[picks[1]]};
}

void genetic(const string &out, int C, const VI &ce, const VI &ne, const VC &classes, const VI &cleft)
{
  double start = now();

  MP parents = generate_parents(C, count_cleft(classes), ce, ne, classes);
  int termination_conditions = 100000;

  update_fitness(parents);
  sort_parents(parents);

  Parent best_individual = parents[0];

  // termination conditions: we stop if in the last -termination_condition- generations
  // it hasn't been an improvement in the solution
  int tc = termination_conditions;

  while (tc > 0 and (now() - start) < 60)
  {

    MP selected_parents = roulette_wheel_selection(parents);
    VI SP1 = selected_parents[0].solution;
    VI SP2 = selected_parents[1].solution;

    MP children = recombination(SP1, SP2, cleft, ce, ne, classes);

    mutation(children);

    for (int i = 0; i < 5; ++i)
      parents.push_back(children[i]);

    update_fitness(parents);
    sort_parents(parents);

    for (int i = 0; i < 5; ++i)
      parents.pop_back();

    if (improved(best_individual, parents[0]))
      tc = termination_conditions;
    else
      --tc;
    }

  write_solution(best_individual.solution, best_individual.penalization, now() - start, out);
}

void read_input(is &in, int &C, int &M, int &K, VI &ce, VI &ne, VC &classes)
{
  in >> C >> M >> K;

  ce.resize(M);
  ne.resize(M);
  classes.resize(K);

  for (auto &capacity : ce)
    in >> capacity;
  for (auto &window : ne)
    in >> window;

  for (int i = 0; i < K; ++i)
  {
    in >> classes[i].id >> classes[i].ncars;
    VI improvements(M);
    for (auto &bit : improvements)
      in >> bit;
    classes[i].improvements = improvements;
  }
}

int main(int argc, char *argv[])
{

  if (argc != 3)
  {
    cout << "Syntax: " << argv[0] << " input_file output_file" << endl;
    exit(1);
  }

  ifstream in(argv[1]);

  int C, M, K;
  VI ce, ne;
  VC classes;

  read_input(in, C, M, K, ce, ne, classes);

  genetic(argv[2], C, ce, ne, classes, count_cleft(classes));
}