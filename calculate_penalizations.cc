#include <iostream>
#include <vector>
#include <climits>
#include <queue>
#include <fstream>
#include <sstream>
using namespace std;

using VI = vector<int>;

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

// falta sumar la part del final que no te mida n!!!
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

int count_pen(const vector<int> &sol, int C, const VI &ce, const VI &ne, const vector<Class> &classes)
{
  int K = classes.size();
  int M = ne.size();
  int total_pen = 0;
  vector<Pen> pens(M);
  for (int i = 0; i < M; ++i)
  {
    pens[i].sum = 0;
    for (int j = 0; j < ne[i]; ++j)
      pens[i].q.push(0);
  }
  for (int i = 0; i < C; ++i)
  {
    total_pen += count_pen_tot(classes[sol[i]].imp, pens, ce);
  }
  return total_pen;
}

void read_input(ifstream &f, int &C, int &M, int &K, VI &ce, VI &ne,
                vector<Class> &classes, vector<vector<int>> &solutions)
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

  string s;
  getline(f, s);
  while (getline(f, s))
  {
    vector<int> sol(C);
    istringstream iss(s);
    for (int &x : sol)
      f >> x;
    solutions.push_back(sol);
  }
}

int main(int argc, char *argv[])
{

  /*

  Input:

  10 5 3
  1 1 1 2 1
  2 2 2 3 2
  0 4 1 1 0 0 1
  1 3 0 1 0 1 0
  2 3 0 0 1 0 0
                        <- deixar una separació
  0 1 0 1 2 0 2 0 2 1
         .
         .              <- posar tantes solucions com es vulguin (que compleixin els requeriments)
         .
  0 1 4 1 2 0 0 0 2 1

  Output(terminal):

    cost primera solucio
            .
            .
            .
    cost última solució

  */

  ifstream f(argv[1]);

  int C, M, K;
  VI ce, ne;
  vector<Class> classes;

  vector<vector<int>> solutions;
  read_input(f, C, M, K, ce, ne, classes, solutions);

  for (vector<int> &sol : solutions)
    cout << count_pen(sol, C, ce, ne, classes) << endl;
}
