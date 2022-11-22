#include <iostream>
#include <vector>
#include <climits>
#include <chrono>
#include <queue>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <random>
using namespace std::chrono;
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

void restore(VS &stations, const VI &ne)
{
  int M = ne.size();
  for (int i = 0; i < M; ++i)
  {
    int n = stations[i].line.size();
    stations[i].requirements -= stations[i].line[n - 1];
    if (n - ne[i] >= 0)
      stations[i].requirements += stations[i].line[n - ne[i]];
    stations[i].line.pop_back();
  }
}

/* Updates each individual line and retuns the penalitzations of it */
int update_station(int bit, Station &st, int ce, int ne, bool end)
{
  st.line.push_back(bit);
  st.requirements += bit;

  // update the window
  int uw = st.line.size() - ne - 1;
  if (uw >= 0)
    st.requirements -= st.line[uw];

  int penalitzation = max(st.requirements - ce, 0);
  // we add to the penalitzarions the final windows
  if (end)
  {
    ++uw;
    --ne;
    for (; st.requirements > 0 and ne > 0; ++uw, --ne)
    {
      st.requirements -= st.line[uw];
      penalitzation += max(st.requirements - ce, 0);
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
  else
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
  int K = classes.size();
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
