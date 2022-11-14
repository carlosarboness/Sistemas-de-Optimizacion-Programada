#include <iostream>
#include <vector>
#include <climits>
#include <chrono>
#include <iomanip>
#include <queue>
#include <cassert>
#include <fstream>
using namespace std::chrono;
using namespace std;

using VI = vector<int>;

double now()
{
  return clock() / double(CLOCKS_PER_SEC);
}

struct Class
{
  int id, n;
  VI imp; // improvements
};

struct Pen
{
  int sum;      // suma total de 1's que hi ha dins la cua
  queue<int> q; // cua que conté els últims ne elements
};

void write_solution(const VI &best_sol, int pen, const double &time, ofstream &s)
{
  s << pen << " " << setprecision(1) << time << endl;
  bool primer = true;
  for (auto &b : best_sol)
  {
    if (primer)
      primer = false;
    else
      s << " ";
    s << b;
  }
  s << endl;
}

int count_pen_millora(int improv, Pen &pena, int ce_i)
{
  pena.sum += improv;
  pena.sum -= pena.q.front();
  pena.q.pop();
  return max(pena.sum - ce_i, 0);
}

// falta sumar la part del final que no te mida n!!!
int count_pen_tot(const VI &imp, vector<Pen> &pens, const VI ce)
{
  int M = pens.size();
  int pen_tot = 0;
  for (int i = 0; i < M; ++i)
    pen_tot += count_pen_millora(imp[i], pens[i], ce[i]);
  return pen_tot;
}

void exh_rec(int i, int current_pen, int min_pen, VI &best_sol, VI &cars_left, vector<Pen> &pens,
             const VI ce, const VI ne, const vector<Class> &classes, const double &start, ofstream &s)
{
  int C = best_sol.size();
  int K = classes.size();
  if (current_pen >= min_pen)
    return;
  if (i == C)
  {
    double end = now();
    double elapsed_time = end - start;
    write_solution(best_sol, current_pen, elapsed_time, s);
  }
  else
  {
    for (int j = 0; j < K; ++j)
    {
      if (cars_left[j] > 0)
      {
        --cars_left[j];
        best_sol[i] = j;
        exh_rec(i + 1, current_pen + count_pen_tot(classes[j].imp, pens, ce),
                min_pen, best_sol, cars_left, pens, ce, ne, classes, start, s);
        ++cars_left[j];
      }
    }
  }
}

void exh(ofstream &s, int C, const VI &ce, const VI &ne, const vector<Class> &classes)
{
  int K = classes.size();
  int M = ne.size();
  VI best_sol(C);
  VI cars_left(K);
  for (int i = 0; i < K; ++i)
    cars_left[i] = classes[i].n;
  vector<Pen> pens(M);
  for (int i = 0; i < M; ++i)
  {
    pens[i].sum = 0;
    for (int j = 0; j < ne[i]; ++j)
      pens[i].q.push(0);
    // assert(pens[i].q.size() == ne);
  }
  double start = now();
  exh_rec(0, 0, INT_MAX, best_sol, cars_left, pens, ce, ne, classes, start, s);
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
  string fe = argv[1], fs = argv[2];

  ifstream f(fe);
  ofstream s(fs);

  int C, M, K;
  VI ce, ne;
  vector<Class> classes;

  read_input(f, C, M, K, ce, ne, classes);

  exh(s, C, ce, ne, classes);
}
