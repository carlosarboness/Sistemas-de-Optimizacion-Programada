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
using MI = vector<VI>;

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

int count_pen_millora(int improv, Pen &pena, int ce_i)
{
  pena.sum -= pena.q.front();
  pena.q.pop();
  pena.q.push(improv);
  pena.sum += improv;
  return max(pena.sum - ce_i, 0);
}

int count_pen_tot(const VI &imp, vector<Pen> &pens, const VI ce)
{
  int M = pens.size();
  int pen_tot = 0;
  for (int i = 0; i < M; ++i)
  {
    int x = count_pen_millora(imp[i], pens[i], ce[i]);
    pen_tot += x;
  }
  return pen_tot;
}

int min_VI(const VI &v)
{
  int n = v.size();
  int min_elem = INT_MAX;
  for (int i = 0; i < n; ++i)
    if (v[i] < min_elem)
      min_elem = v[i];
  return min_elem;
}

int lower_bound(int i, int current_pen, const VI &current_sol, const VI &cars_left, const MI &costs)
{
  if (i == 0)
    return -1;
  int lb = current_pen;
  int K = cars_left.size();
  for (int j = 0; j < K; ++j)
    if (cars_left[j] > 0)
    {
      int min = INT_MAX;
      for (int k = 0; k < K; ++k)
        if (cars_left[k] > 0 and costs[j][k] < min)
          min = costs[j][k];
      lb += min * cars_left[j];
    }
  return lb;
}

void exh_rec(int i, int current_pen, int &min_pen, VI &current_sol, VI &cars_left, vector<Pen> &pens,
             const VI ce, const VI ne, const vector<Class> &classes, const MI &costs, const double &start, const string &s)
{
  int C = current_sol.size();
  int K = classes.size();
  if (current_pen >= min_pen)
    return;
  if (i == C)
  {
    min_pen = current_pen;

    double end = now();
    double elapsed_time = end - start;
    write_solution(current_sol, current_pen, elapsed_time, s);
  }
  else if (lower_bound(i, current_pen, current_sol, cars_left, costs) < min_pen)
  {
    for (int j = 0; j < K; ++j)
    {
      if (cars_left[j] > 0)
      {
        --cars_left[j];
        current_sol[i] = j;
        vector<Pen> pens_rec = pens;
        exh_rec(i + 1, current_pen + count_pen_tot(classes[j].imp, pens, ce),
                min_pen, current_sol, cars_left, pens, ce, ne, classes, costs, start, s);
        pens = pens_rec;
        ++cars_left[j];
      }
    }
  }
}

int calculate_cost(const VI &ne, const VI &imp_i, const VI &imp_j)
{
  int cost = 0;
  int M = imp_i.size();
  for (int k = 0; k < M; ++k)
    if (imp_i[k] == imp_j[k] and ne[k] == 2)
      ++cost;
  return cost;
}

void exh(const string &s, int C, const VI &ce, const VI &ne, const vector<Class> &classes)
{
  int K = classes.size();
  int M = ne.size();
  VI current_sol(C);
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
  MI costs(K, VI(K));
  for (int i = 0; i < K; ++i)
    for (int j = i; j < K; ++j)
      costs[i][j] = costs[j][i] = calculate_cost(ne, classes[i].imp, classes[j].imp);

  double start = now();
  int min_pen = INT_MAX;
  exh_rec(0, 0, min_pen, current_sol, cars_left, pens, ce, ne, classes, costs, start, s);
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

  exh(argv[2], C, ce, ne, classes);
}
