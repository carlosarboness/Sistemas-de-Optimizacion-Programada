#include <iostream>
#include <vector>
#include <climits>
using namespace std;

using VI = vector<int>;

struct Class
{
  int id, n;
  VI imp; // improvements
};

void exh_rec(int i, int current_pen, int min_pen, VI &best_sol, VI &cars_left,
             const VI &ce, const VI &ne, const vector<Class> &classes)
{
  int C = best_sol.size();
  int K = classes.size();
  if (current_pen >= min_pen)
    return;
  if (i == C)
    write_solution(best_sol);
  else
  {
    for (int j = 0; j < K; ++j)
    {
      if (cars_left[j] > 0)
      {
        --cars_left[j];
        best_sol[i] = j;
        exh_rec(i + 1, current_pen + count_pen(), min_pen, best_sol, cars_left, ce, ne, classes);
        ++cars_left[j];
      }
    }
  }
}

void exh(int C, const VI &ce, const VI &ne, const vector<Class> &classes)
{
  int K = classes.size();
  VI best_sol(C);
  VI cars_left(K);
  for (int i = 0; i < K; ++i)
    cars_left[i] = classes[i].n;
  exh_rec(0, 0, INT_MAX, best_sol, cars_left, ce, ne, classes);
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