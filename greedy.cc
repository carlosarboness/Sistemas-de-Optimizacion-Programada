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

void greedy()
{
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

  greedy();
}
