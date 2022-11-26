#include <iostream>
#include <vector>
#include <climits>
#include <chrono>
#include <cassert>
#include <fstream>
using namespace std;

struct Class;
struct Station;

using VI = vector<int>;
using VC = vector<Class>;
using VS = vector<Station>;
using is = ifstream;
using os = ofstream;

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

void restore(VS &st, const VI &ne)
{
  int M = ne.size();
  for (int i = 0; i < M; ++i)
  {
    int n = st[i].line.size();
    st[i].requirements -= st[i].line[n - 1];
    if (n - ne[i] - 1 >= 0)
      st[i].requirements += st[i].line[n - ne[i] - 1];
    st[i].line.pop_back();
  }
}

int geom_sum(int n)
{
  return (n * (n + 1)) / 2;
}

/* Updates each individual line and retuns the penalitzations of it */
int update_station(int bit, Station &st, int ce, int ne, bool end)
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
    int rec = st.requirements;
    for (; rec > 0; ++uw)
    {
      rec -= st.line[uw];
      penalitzation += max(rec - ce, 0);
    }
  }
  return penalitzation;
}

/* Update Production Lines: updates the production lines of all the stations because of the
addition of a new car in the production line. The function returns the total sum of penalizations
of each line, updated with the new car */
int UPL(const VI &improvements, VS &stations, const VI &ce, const VI &ne, bool end)
{
  int M = stations.size();
  int total_penalitzations = 0;
  for (int i = 0; i < M; ++i)
    total_penalitzations += update_station(improvements[i], stations[i], ce[i], ne[i], end);
  return total_penalitzations;
}

int calculate_ones(int j, const VI &cleft, const VC &classes)
{
  int K = cleft.size();
  int ones = 0;
  for (int i = 0; i < K; ++i)
    ones += cleft[i] * classes[i].improvements[j];
  return ones;
}

int lb_station(int j, int i, const VI &cs, const VI &cleft, int ce, int ne, const VC &classes)
{
  int C = cs.size();
  int initial_zeros = 0;
  for (int k = i; k >= 0 and classes[cs[k]].improvements[j] == 0; --k)
    ++initial_zeros;

  initial_zeros = max(ne - ce, 0);
  int ones = calculate_ones(j, cleft, classes);
  int zeros = (C - i - 1) - ones - initial_zeros;

  if (zeros < 0)
    return (ones - ne + 1) * (ne - ce) + 2 * geom_sum(ne - ce - 1); //+ geom_sum(-zeros);

  ones -= ce;
  while (zeros > 0 and ones > 0)
  {
    zeros -= (ne - ce);
    ones -= ce;
  }

  if (zeros <= 0 and ones > 0)
    return (ones - ne + 1) * (ne - ce) + 2 * geom_sum(ne - ce - 1) + geom_sum(-zeros);

  return 0;
}

int lower_bound(int i, int cp, const VI &cs, const VI &cleft, const VI &ce, const VI &ne, const VC &classes)
{
  if (i < 0)
    return 0;
  int M = ne.size();
  int lb = cp;
  for (int j = 0; j < M; ++j)
    lb += lb_station(j, i, cs, cleft, ce[j], ne[j], classes);
  return lb;
}

void exhaustive_search_rec(int i, int cp, int &mp, VI &cs, VI &cleft, VS &stations, const VI &ce,
                           const VI &ne, const VC &classes, const double &start, const string &out)
{
  int C = cs.size();
  int K = classes.size();

  if (i == C)
  {
    mp = cp;
    write_solution(cs, cp, now() - start, out);
  }
  if (lower_bound(i - 1, cp, cs, cleft, ce, ne, classes) < mp)
  {
    for (int cl = 0; cl < K; ++cl)
    {
      if (cleft[cl] > 0)
      {
        --cleft[cl];
        cs[i] = cl;

        // VS stationsr = stations;

        if (int up = UPL(classes[cl].improvements, stations, ce, ne, i + 1 == C); up + cp < mp)
          exhaustive_search_rec(i + 1, cp + up, mp, cs, cleft, stations, ce, ne, classes, start, out);

        // stations = stationsr;
        restore(stations, ne);
        ++cleft[cl];
      }
    }
  }
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

/* Returns a vector containing for each index the number of cars that are left for that class */
VI count_cleft(const VC &classes)
{
  int K = classes.size();
  VI cleft(K);
  for (int i = 0; i < K; ++i)
    cleft[i] = classes[i].ncars;
  return cleft;
}

void exhaustive_search(int C, const VI &ce, const VI &ne, const VC &classes, const string &out)
{
  int M = ce.size();

  VI cs(C); // current solution
  VI cleft = count_cleft(classes);
  VS stations = inicialize_stations(C, M);
  int mp = INT_MAX; // minimum penalization

  double start = now(); // inicialize the counter
  exhaustive_search_rec(0, 0, mp, cs, cleft, stations, ce, ne, classes, start, out);
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

  exhaustive_search(C, ce, ne, classes, argv[2]);
}