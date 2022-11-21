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

/*
Writes the solution in corrent_sol with penalitation pen that took time seconts to be
found in the ouput file s with the appropiate format
*/
void write_solution(const VI &current_sol, int pen, const double &time, const string &s)
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

/*
Returns the aproximate cost of going from a car with the improvements of imp_i to a
car with the improvemnts of imp_j. The aproximate cost is the sum of ne/ce ot the
improvments that both cars have.
*/
double calculate_cost(const VI &ce, const VI &ne, const VI &imp_i, const VI &imp_j)
{
  double cost = 0;
  int M = imp_i.size();
  for (int k = 0; k < M; ++k)
    if (imp_i[k] == imp_j[k])
      cost += ne[k] / ce[k];
  return cost;
}

/*
Returns the argmax of the vector v
*/
int argmax_VI(const VI &v)
{
  int n = v.size();
  int k = -1;
  int max_elem = 0;
  for (int i = 0; i < n; ++i)
    if (v[i] > max_elem)
    {
      max_elem = v[i];
      k = i;
    }
  return k;
}

/*
Returns the class of the car that is written next in the solution following the greedy algorithm
*/
int conditions(int last_car, VI &cars_left, const VD &costs)
{
  int K = cars_left.size();
  int cl = 0;
  double cost = -1;
  int id_car = -1;
  for (int i = 0; i < K; ++i)
  {

    if (cars_left[i] > cl or (cars_left[i] == cl and costs[i] < cost))
    {
      id_car = i;
      cl = cars_left[i];
      cost = costs[i];
    }
  }
  --cars_left[id_car];
  return id_car;
}

/*
Returns the vector with the solution following the greedy algorithm
*/
VI gen_sol(int C, const VI &ce, const VI &ne, const vector<Class> &classes)
{
  // Generate data structues
  int K = classes.size();
  VI cars_left(K);
  for (int i = 0; i < K; ++i)
    cars_left[i] = classes[i].n;
  MD costs(K, VD(K)); // matrix where in the position [i][j] stores the approximate cost of going from i to j
  for (int i = 0; i < K; ++i)
    for (int j = i; j < K; ++j)
      costs[i][j] = costs[j][i] = calculate_cost(ce, ne, classes[i].imp, classes[j].imp);
  VI solution(C);

  // Calculate solution
  solution[0] = argmax_VI(cars_left);
  --cars_left[solution[0]];
  for (int k = 1; k < C; ++k)
    solution[k] = conditions(solution[k - 1], cars_left, costs[solution[k - 1]]);
  return solution;
}

/*
Returns the penalization of adding the improvment improv
and updates the Pen pena
*/
int count_pen_millora(int improv, Pen &pena, int ce_i, bool final)
{
  pena.sum -= pena.q.front();
  pena.q.pop();
  pena.q.push(improv);
  pena.sum += improv;
  int pen = max(pena.sum - ce_i, 0);
  if (final and pen > 0) // if is the final car of the solution, the queue is emptied
  {
    while (not pena.q.empty())
    {
      pena.sum -= pena.q.front();
      pena.q.pop();
      pen += max(pena.sum - ce_i, 0);
    }
  }
  return pen;
}

/*
Returns the penalization of adding in the solution the class
with improvment imp and updates the vector of Pen pens
*/
int count_pen_tot(const VI &imp, vector<Pen> &pens, const VI ce, bool final)
{
  int M = pens.size();
  int pen_tot = 0;
  for (int j = 0; j < M; ++j)
    pen_tot += count_pen_millora(imp[j], pens[j], ce[j], final);
  return pen_tot;
}

/*
Returns the penalization of the soltion sol
*/
int count_pen(const vector<int> &sol, int C, const VI &ce, const VI &ne, const vector<Class> &classes)
{
  int M = ne.size();
  int total_pen = 0;
  vector<Pen> pens(M); // Data structure to calculate the penalizations
  for (int i = 0; i < M; ++i)
  {
    pens[i].sum = 0;
    for (int j = 0; j < ne[i]; ++j)
      pens[i].q.push(0);
  }
  for (int i = 0; i < C; ++i)
    total_pen += count_pen_tot(classes[sol[i]].imp, pens, ce, i == C - 1);
  return total_pen;
}

/*
Given an imput in the data structures, writes the solution following the greedy algorithm with
its cost and the time it took to calculate it

The greedy algorithm is the following:
  - A car of the class with more cars left is added to the solution. If there is more than one maximum
    the car added is the one that has the less aproximate cost given the last car added to the solution.
    If there is still a draw, the car with the smallest id is added.
*/
void greedy(const string &s, int C, const VI &ce, const VI &ne, const vector<Class> &classes)
{
  double start = now();
  VI solution = gen_sol(C, ce, ne, classes);
  double end = now();
  int pen = count_pen(solution, C, ce, ne, classes);
  double elapsed_time = end - start;
  write_solution(solution, pen, elapsed_time, s);
}

/*
Reads the input from the input file f and stores the data in the corresponding data structure
*/
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