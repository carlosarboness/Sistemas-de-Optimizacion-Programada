#include <iostream>
#include <vector>
#include <climits>
#include <chrono>
#include <cassert>
#include <fstream>
#include <algorithm>
using namespace std;

struct Class;
struct Station;

using VI = vector<int>;
using VD = vector<double>;
using VC = vector<Class>;
using VS = vector<Station>;
using MD = vector<VD>;
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

struct fit
{
  int id_car = UNDEF, cl = UNDEF;
  double cost = UNDEF;
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

/* Updates each individual line and retuns the penalitzations of it */
int add_car(int bit, Station &st, int ce, int ne, bool end)
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

int update_station_cs(const VI &imp, VS &stations, const VI &ce, const VI &ne, bool end)
{
  int M = ce.size();
  int total_penalization = 0;
  for (int i = 0; i < M; ++i)
    total_penalization += add_car(imp[i], stations[i], ce[i], ne[i], end);
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

int count_penalization_cs(const VI &sol, int C, const VI &ce, const VI &ne, const VC &classes)
{
  int M = ne.size();
  VS stations = inicialize_stations(C, M);

  int total_penalization = 0;
  for (int i = 0; i < C; ++i)
    total_penalization += update_station_cs(classes[sol[i]].improvements, stations, ce, ne, i == C - 1);
  return total_penalization;
}

// returns the argmax of V
int argmax(const VI &vi, const VD &vd)
{
  if (vi.size() != 0)
  {
    auto it = max_element(vi.begin(), vi.end());
    return distance(vi.begin(), it);
  }
  else
  {
    auto it = max_element(vd.begin(), vd.end());
    return distance(vd.begin(), it);
  }
}

void update_fitness(int id_car, MD &fitness, VI &cleft)
{
  int K = cleft.size();
  for (int i = 0; i < K; ++i)
  {
    double c = cleft[id_car] / fitness[i][id_car];
    fitness[i][id_car] = (cleft[id_car] - 1) / c;
  }
  --cleft[id_car];
}

int best(const VD &f, VI &cleft)
{
  int id = UNDEF;
  double fit = UNDEF;
  for (int i = 0; i < cleft.size(); ++i)
  {
    if (f[i] > fit or (f[i] == fit and cleft[i] > cleft[id]))
    {
      id = i;
      fit = f[i];
    }
  }
  return id;
}

int best_fit(VI &cleft, MD &fitness, int id_last_car)
{

  int id_car = best(fitness[id_last_car], cleft);
  update_fitness(id_car, fitness, cleft);
  return id_car;
}

double calculate_cost(const VI &ce, const VI &ne, const VI &imp_i, const VI &imp_j)
{
  double cost = 0;

  for (int k = 0; k < ce.size(); ++k)
    if (imp_i[k] == 1 and imp_j[k] == 1)
      cost += ne[k] / ce[k];

  return cost;
}

MD generate_fitness_matrix(int K, const VI &cleft, const VI &ce, const VI &ne, const VC &classes)
{
  MD fitness(K, VD(K));
  for (int i = 0; i < K; ++i)
    for (int j = 0; j < K; ++j)
    {
      double c = calculate_cost(ce, ne, classes[i].improvements, classes[j].improvements);
      fitness[i][j] = (double)cleft[j] / (c == 0 ? 1 : c);
    }
  return fitness;
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

VI base_case(int C, VI &cleft, MD &fitness)
{
  VI solution(C);
  int base = argmax(cleft, {});
  solution[0] = base;
  update_fitness(base, fitness, cleft);
  return solution;
}

VI gen_sol(int C, const VI &ce, const VI &ne, const VC &classes)
{

  VI cleft = count_cleft(classes);
  MD fitness = generate_fitness_matrix(classes.size(), cleft, ce, ne, classes);
  // Calculate solution
  VI solution = base_case(C, cleft, fitness);

  for (int k = 1; k < C; ++k)
    solution[k] = best_fit(cleft, fitness, solution[k - 1]);

  return solution;
}

void greedy(const string &s, int C, const VI &ce, const VI &ne, const vector<Class> &classes)
{
  double start = now();
  VI solution = gen_sol(C, ce, ne, classes);
  double end = now();
  int penalization = count_penalization_cs(solution, C, ce, ne, classes);
  write_solution(solution, penalization, end - start, s);
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

  greedy(argv[2], C, ce, ne, classes);
}

  

  
   
