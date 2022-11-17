#include <iostream>
#include <vector>
#include <climits>
#include <chrono>
#include <queue>
#include <cassert>
#include <fstream>
using namespace std::chrono;
using namespace std;

using VI = vector<int>;
using VD = vector<double>;
using MD = vector<VD>;

double now()
{
  return clock() / double(CLOCKS_PER_SEC);
}

ofstream open(const string &s)
{
  ofstream f;
  f.open(s, ofstream::out | ofstream::trunc);
  f.setf(ios::fixed);
  f.precision(1);
  return f;
}

struct Class
{
  int id, n; // id and # of cars of the Class
  VI imp;    // improvements
};

struct Pen
{
  int sum;      // suma total de 1's que hi ha dins la cua
  queue<int> q; // cua que conté els últims ne elements
};

struct lakjdfaklj
{
  int id, left;
  int costs;
}

void
write_solution(const VI &current_sol, int pen, const double &time, const string &s)
{
  ofstream f = open(s);
  f << pen << " " << time << endl;
  bool primer = true;
  for (auto &b : current_sol)
  {
    if (primer)
      primer = false;
    else
      f << " ";
    f << b;
  }
  f << endl;
  f.close();
}

double calculate_cost(const VI &ce, const VI &ne, const VI &imp_i, const VI &imp_j)
{
  double cost = 0;
  int M = imp_i.size();
  for (int k = 0; k < M; ++k)
    if (imp_i[k] == imp_j[k])
      if (ne[k] == 2)
        cost += 1;
      else
        cost += ce[k] / ne[k];
  return cost;
}

VI gen_sol(int C, int &pen, const VI &ce, const VI &ne, const vector<Class> &classes, const string &s)
{
  // Generate data structues
  int K = classes.size();
  int M = ne.size();
  VI solution(C);
  VI cars_left(K);
  for (int i = 0; i < K; ++i)
    cars_left[i] = classes[i].n;
  vector<Pen> pens(M);
  for (int i = 0; i < M; ++i)
  {
    pens[i].sum = 0;
    for (int j = 0; j < ne[i]; ++j)
      pens[i].q.push(0);
  }
  MD costs(C, VD(C));
  for (int i = 0; i < C; ++i)
    for (int j = i; j < C; ++j)
      costs[i][j] = costs[j][i] = calculate_cost(ce, ne, classes[i].imp, classes[j].imp);

  // Calculate solution

  /* Algorisme:
  1r:
  int k = 0
  id = car = 0
  for car_it in cars_left:
    if car_it > car:
      car = car_it
      id = i
  solution[k] = id;
  ++k;
  2n:
  //calculate ce_max
  while(k <= ce_max){
    solution[k] = priorize_id_more_cars_left
    if equal then priorize_id_less_cost
    if equal then smaller_id
    ++k;
    actualize_pen;
  }
  3r:
  while(k < C){
    solution[k] = priorize_id_less_cost
    if equal then priorize_id_more_cars_left
    if equal then smaller_id
    ++k;
    actualize_pen;
  }
  */
  return solution;
}

void greedy(const string &s, int C, const VI &ce, const VI &ne, const vector<Class> &classes)
{
  int pen = 0;
  double start = now();
  VI solution = gen_sol(C, pen, ce, ne, classes, start, s);
  double end = now();
  double elapsed_time = end - start;
  write_solution(solution, pen, elapsed_time, s);
}

void read_input(ifstream &f, int &C, int &M, int &K, VI &ce, VI &ne, vector<Class> &classes)
{
  f >> C >> M >> K;
  ce.resize(M);
  ne.resize(M);
  classes.resize(K);
  for (auto &capacity : ce)
    f >> capacity;
  for (auto &num_cars : ne)
    f >> num_cars;
  for (int i = 0; i < K; ++i)
  {
    f >> classes[i].id >> classes[i].n;
    VI imp(M);
    for (auto &imprv : imp)
      f >> imprv;
    classes[i].imp = imp;
  }
}

int main(int argc, char *argv[])
{

  assert(argc == 3);

  ifstream f(argv[1]);

  int C, M, K;
  VI ce, ne;
  vector<Class> classes;

  read_input(f, C, M, K, ce, ne, classes);

  greedy(argv[2], C, ce, ne, classes);
}