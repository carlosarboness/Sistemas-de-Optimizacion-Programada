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
  int ones = 0;
  for (int i = 0; i < cleft.size(); ++i)
    ones += cleft[i] * classes[i].improvements[j];
  return ones;
}

int count_penalization(const VI &line, int ce, int ne)
{

  int total_penalizations, requirements;
  total_penalizations = requirements = 0;

  // count inicial windows
  int uw = 0;
  for (int i = ne; i > 0; --i)
  {
    requirements += line[uw];
    total_penalizations += max(requirements - ce, 0);
    ++uw;
  }

  // count mid penalizations
  for (; uw < line.size(); ++uw)
  {
    requirements += line[uw];
    requirements -= line[uw - ne];
    total_penalizations += max(requirements - ce, 0);
  }

  // count end windows
  while (requirements > 0)
  {
    requirements -= line[uw - ne];
    total_penalizations += max(requirements - ce, 0);
    ++uw;
  }

  return total_penalizations;
}

int lb_station(int j, int i, const VI &cs, const VI &cleft, VI line, int ce, int ne, const VC &classes)
{
  int C = cs.size();
  int ones = calculate_ones(j, cleft, classes);
  int zeros = (C - i - 1) - ones;

  if (i < ne)
  {
    int first_ones = 0;
    for (int k = 0; k <= i; ++k)
      if (line[k])
        ++first_ones;

    for (int k = ce - first_ones; k > 0 and ones > 0; --k)
    {
      --ones;
      line.push_back(1);
      ++i;
    }
  }
  for (; i < ne; ++i)
  {
    --zeros;
    line.push_back(0);
  }

  while (zeros > 0 and ones > 0)
  {
    int bit = line[i - ne + 1];
    line.push_back(bit);
    (bit ? --ones : --zeros);
    ++i;
  }

  for (; zeros > 0; --zeros)
    line.push_back(0);

  for (; ones > 0; --ones)
    line.push_back(1);

  if (cs[0] == 2 and cs[1] == 1 and cs[2] == 3 and cs[3] == 1 and cs[4] == 2 and cs[5] == 0 and cs[6] == 2)
  {
    cout << j << ":   ";

    for (int s = 0; s <= i; ++s)
      cout << line[s] << " ";

    cout << " | ";

    for (int s = i + 1; s < line.size(); ++s)

      cout << line[s] << " ";
    cout << "  ce: " << ce << "  ne: " << ne << "   pen: " << count_penalization(line, ce, ne) << endl;
  }

  return count_penalization(line, ce, ne);
}

int lower_bound(int i, int cp, const VI &cs, const VI &cleft, const VS &st, const VI &ce, const VI &ne, const VC &cl)
{
  if (i == -1)
    return 0;

  int lb = 0;
  for (int j = 0; j < ne.size(); ++j)
    lb += lb_station(j, i, cs, cleft, st[j].line, ce[j], ne[j], cl);

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
  if (lower_bound(i - 1, cp, cs, cleft, stations, ce, ne, classes) < mp)
  {
    for (int cl = 0; cl < K; ++cl)
    {
      if (cleft[cl] > 0)
      {
        --cleft[cl];
        cs[i] = cl;

        if (int up = UPL(classes[cl].improvements, stations, ce, ne, i + 1 == C); up + cp < mp)
          exhaustive_search_rec(i + 1, cp + up, mp, cs, cleft, stations, ce, ne, classes, start, out);

        restore(stations, ne);
        ++cleft[cl];
      }
    }
  }
  /*
  else
  {
    for (int s = 0; s < i; ++s)
    {
      cout << cs[s] << " ";
    }
    cout << endl;
    cout << lower_bound(i - 1, cp, cs, cleft, stations, ce, ne, classes) << " > " << mp << endl
         << endl;
  }
  */
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