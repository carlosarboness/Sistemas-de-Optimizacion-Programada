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
using VD = vector<double>;
using VC = vector<Class>;
using VS = vector<Station>;
using MI = vector<VI>;
using MP = vector<Parent>;
using MD = vector<VD>;
using VB = vector<bool>;
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
  assert(parent.size() == C);
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
  VI parent = create_parent(C, cleft);

  // select the population size (it will remain constant during the iterations)
  int psize = 50;

  MP parents(psize);

  for (int i = 0; i < psize; ++i)
  {
    VI P = generate_permutation(parent);
    parents[i] = Parent{P, fitness(P, ce, ne, classes)};
  }

  return parents;
}

Parent find_best_individual(const MP &parents)
{
  Parent bi = Parent{{}, INT_MAX};
  for (int i = 0; i < parents.size(); ++i)
    if (parents[i].fitness < bi.fitness)
      bi = parents[i];
  return bi;
}

bool best_fit(const Parent &P1, const Parent &P2)
{
  return P1.fitness <= P2.fitness;
}

MP select_parents(MP &parents)
{
  sort(parents.begin(), parents.end(), best_fit);
  return parents;
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

void adapt_solution(const VI &P, VI &child, VI &cleft, int left, int right)
{
  int C = P.size();
  int j = 0;

  while (j < C)
  {
    if (j == left)
    {
      j = right + 1;
      if (j == C)
        break;
    }

    int a = child[j];
    if (cleft[a] > 0)
      --cleft[a];
    else
    {
      child[j] = most_cleft_index(cleft);
      --cleft[child[j]];
    }
    ++j;
  }
}

MI recombinate(const VI &P1, const VI &P2, int left, int right, const VI &cleft)
{
  int C = P1.size();

  VI children1 = P1;
  VI cleft1 = cleft;
  VI children2 = P2;
  VI cleft2 = cleft;

  for (int i = left; i <= right; ++i)
  {
    children1[i] = P2[i];
    --cleft1[P2[i]];
    children2[i] = P1[i];
    --cleft2[P1[i]];
  }

  adapt_solution(P1, children1, cleft1, left, right);
  adapt_solution(P2, children2, cleft2, left, right);

  return {children1, children2};
}

MP recombination(const VI &P1, const VI &P2, const VI &cleft, const VI &ce, const VI &ne, const VC &classes)
{
  int C = P1.size();

  VI cuts = generate_random_cut_points(C);

  MI recombination = recombinate(P1, P2, cuts[0], cuts[1], cleft);

  Parent children1 = Parent{recombination[0], count(recombination[0], ce, ne, classes)};
  Parent children2 = Parent{recombination[1], count(recombination[1], ce, ne, classes)};

  return {children1, children2};
}

void swap(const VI &ch, VI &c, int i, int j)
{
  c[i] = ch[j];
  c[j] = ch[i];
}

void mutate(MP &children)
{
  int C = children[0].solution.size();
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> distr(0, C - 1);

  VI c1 = children[0].solution;
  VI c2 = children[1].solution;

  int i = distr(gen);
  int j = distr(gen);

  swap(c1, c1, i, j);
  swap(c2, c2, i, j);

  children[0].solution = c1;
  children[1].solution = c2;
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
  int termination_conditions = 100000;

  // termination conditions: we stop if in the last -termination_condition- generations
  // it hasn't been an improvement in the solution
  int tc = termination_conditions;

  while (tc > 0)
  {
    parents = select_parents(parents);

    Parent P1 = parents[0];
    Parent P2 = parents[1];

    MP children = recombination(P1.solution, P2.solution, cleft, ce, ne, classes);

    mutation(children);

    parents[0] = children[0];
    parents[1] = children[1];

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