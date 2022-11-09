#include <iostream>
#include <vector>
#include <climits>
#include <chrono>
#include <iomanip>
#include <queue>
using namespace std::chrono;
using namespace std;

using VI = vector<int>;

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

void write_solution(const VI &best_sol, int pen, double time)
{
  cout << pen << " " << setprecision(1) << time << endl;
  bool primer = true;
  for (auto &s : best_sol)
  {
    if (primer)
      primer = false;
    else
      cout << " ";
    cout << s;
  }
  cout << endl;
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
             const VI ce, const VI ne, const vector<Class> &classes, const auto &start)
{
  int C = best_sol.size();
  int K = classes.size();
  if (current_pen >= min_pen)
    return;
  if (i == C)
  {
    auto now = high_resolution_clock::now();
    write_solution(best_sol, current_pen, now - start);
  }
  else
  {
    for (int j = 0; j < K; ++j)
    {
      if (cars_left[j] > 0)
      {
        --cars_left[j];
        best_sol[i] = j;
        exh_rec(i + 1, current_pen + count_pen_tot(classes[j].imp, pens, ce), min_pen, best_sol, cars_left, ce, ne, classes);
        ++cars_left[j];
      }
    }
  }
}

void exh(int C, const VI &ce, const VI &ne, const vector<Class> &classes)
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
  auto start = high_resolution_clock::now();
  exh_rec(0, 0, INT_MAX, best_sol, cars_left, pens, ce, ne, classes, start);
}

int main()
{
  int C, M, K;
  cin >> C >> M >> K;
  VI ce(M), ne(M);
  for (auto &capacity : ce)
    cin >> capacity;
  for (auto &num_cars : ne)
    cin >> num_cars;
  vector<Class> classes(K);
  for (int i = 0; i < K; ++i)
  {
    cin >> classes[i].id >> classes[i].n;
    VI imp(M);
    for (auto &imprv : imp)
      cin >> imprv;
    classes[i].imp = imp;
  }
  exh(C, ce, ne, classes);
}