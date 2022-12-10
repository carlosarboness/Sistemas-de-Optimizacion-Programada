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
  int fitness;
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

VI generate_permutation(VI parent)
{
  random_shuffle(parent.begin(), parent.end());
  return parent;
}

/* Updates each individual line and retuns the penalitzations of it */
int add(int bit, Station &st, int ce, int ne, bool end)
{

  // update the window
  int uw = st.line.size() - ne;
  if (uw >= 0)
    st.requirements -= st.line[uw];

  st.line.push_back(bit);
  st.requirements += bit;

  int penalitzation = max(st.requirements - ce, 0);
  // we add to the penalitzarions the final windows
  if (end)
  {
    ++uw;
    for (; st.requirements > 0; ++uw)
    {
      st.requirements -= st.line[uw];
      penalitzation += max(st.requirements - ce, 0);
    }
  }
  return penalitzation;
}

int update(const VI &imp, VS &stations, const VI &ce, const VI &ne, bool end)
{
  int M = ce.size();
  int total_penalization = 0;
  for (int i = 0; i < M; ++i)
    total_penalization += add(imp[i], stations[i], ce[i], ne[i], end);
  return total_penalization;
}

/* Returns a vector with the inicialized stations of the algorithm. Firstly the requiremnts are
equal to zero and the line is empty */
VS inicialize_stations(int C, int M)
{
  VS stations(M);

  int requirements = 0;
  VI inicial_line;
  inicial_line.reserve(C);

  for (Station &st : stations)
    st = Station{requirements, inicial_line};

  return stations;
}

int count(const VI &P, const VI &ce, const VI &ne, const VC &classes)
{
  int C = P.size();
  int K = classes.size();
  int M = ne.size();
  VS stations = inicialize_stations(C, M);

  int total_penalization = 0;
  for (int i = 0; i < C; ++i)
    total_penalization += update(classes[P[i]].improvements, stations, ce, ne, i == C - 1);
  return total_penalization;
}

int fitness(const VI &P, const VI &ce, const VI &ne, const VC &classes)
{
  int penalizations = count(P, ce, ne, classes);
  return penalizations;
}

MP generate_parents(int C, const VI &cleft, const VI &ce, const VI &ne, const VC &classes)
{
  VI P = create_parent(C, cleft);

  // select the population size (it will remain constant during the iterations)
  int psize = 400;

  MP parents(psize);

  for (int i = 0; i < psize; ++i)
  {
    P = generate_permutation(P);
    parents[i] = Parent{P, fitness(P, ce, ne, classes)};
  }
  return parents;
}

Parent find_best_individual(const MP &parents)
{
  int n = parents.size();
  Parent bi = Parent{{}, INT_MAX};
  for (int i = 0; i < n; ++i)
  {
    if (parents[i].fitness < bi.fitness)
      bi = parents[i];
  }
  return bi;
}

bool best_fit(const Parent &P1, const Parent &P2)
{
  return P1.fitness < P2.fitness;
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
    children.push_back(Parent{recombination[i], count(recombination[i], ce, ne, classes)});

  return children;
}

void swap(Parent &P, int i, int j)
{

  int z = P.solution[i];
  P.solution[i] = P.solution[j];
  P.solution[j] = z;
}

void mutate(MP &children)
{

  int C = children[0].solution.size();
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> distr(0, C - 1);

  for (int i = 0; i < 5; ++i)
  {
    swap(children[i], distr(gen), distr(gen));
  }
}

void mutation(MP &children)
{
  // add a mutation to the children with propability 0.01

  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> distr(0, 9);
  if (distr(gen) == distr(gen))
    mutate(children);
}

void genetic(const string &out, int C, const VI &ce, const VI &ne, const VC &classes, const VI &cleft)
{
  double start = now();
  MP parents = generate_parents(C, count_cleft(classes), ce, ne, classes);
  Parent best_individual = find_best_individual(parents);
  int termination_conditions = 10000;

  sort_parents(parents);

  // termination conditions: we stop if in the last -termination_condition- generations
  // it hasn't been an improvement in the solution
  int tc = termination_conditions;

  while (tc > 0 and (now() - start) < 60)
  {

    Parent P1 = parents[0];
    Parent P2 = parents[1];

    MP children = recombination(P1.solution, P2.solution, cleft, ce, ne, classes);

    mutation(children);

    for (int i = 0; i < 5; ++i)
      parents.push_back(children[i]);

    sort_parents(parents);

    for (int i = 0; i < 5; ++i)
      parents.pop_back();

    Parent bi = find_best_individual(parents);

    if (bi.fitness < best_individual.fitness)
    {
      best_individual = bi;
      tc = termination_conditions;
    }
    else
      --tc;
  }
  write_solution(best_individual.solution, best_individual.fitness, now() - start, out);
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