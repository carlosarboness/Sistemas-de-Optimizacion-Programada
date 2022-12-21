#include <iostream>
#include <vector>
#include <climits>
#include <chrono>
#include <cassert>
#include <fstream>
#include <numeric>
using namespace std;

struct Station;
using VI = vector<int>;
using MI = vector<VI>;
using VS = vector<Station>;
using VD = vector<double>;
using MD = vector<VD>;

/* ------- GLOBAL VARIABLES ------- */

double start;
string output_file;

/* -------------------------------- */

struct Window
{
  VI ce, ne; // Stores the ce (ne) of a particular station in the position of the vector
};

// Stores all the information given by the input
struct Data
{
  int C, M, K;     // Number of cars, improvements and classes respectively
  Window w;        // Stores the windows informations
  MI improvements; // Matrix of integers (1/0) that indicate if the i-th class requires j-th improvment
  VI cleft;        // Vector that saves the cars left of the class i
};

// Stores the number of needed requirements in a particular
struct Station
{
  int requirements; // number of requirements in window of size ne in a particular moment of the station
  VI line;          // stores the requirements in order that need to be satisfied by the station in a particular moment
};

/* Returns the value in seconds of the processor time consumed by the program. */
double now()
{
  return clock() / double(CLOCKS_PER_SEC);
}

/* Opens the file -output_file- removing all the data that has in it */
ofstream open_truncated_file(const string &output_file)
{
  ofstream out;
  out.open(output_file, ofstream::out | ofstream::trunc);
  out.setf(ios::fixed);
  out.precision(1);
  return out;
}

/* Writes the solution with its penalization and the time used to find that solution (elapsed_time)*/
void write_solution(const VI &solution, int penalization, const double &elapsed_time)
{
  ofstream out = open_truncated_file(output_file);
  out << penalization << " " << elapsed_time << endl;

  out << solution[0];
  for (int i = 1; i < solution.size(); ++i)
    out << " " << solution[i];
  out << endl;

  out.close();
}

// benet
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
// benet

/* Update Production Lines: updates the production lines of all the stations because of the
addition of a new car in the production line.
The function returns the total sum of penalizations of each line, updated with the new car */
int UPL(const VI &imp, VS &stations, const Window &w, bool end)
{
  const auto &[ce, ne] = w;

  int M = stations.size();
  int total_penalitzations = 0;
  for (int i = 0; i < M; ++i)
    total_penalitzations += update_station(imp[i], stations[i], ce[i], ne[i], end);
  return total_penalitzations;
}

/* Returns the penalization of the station line given a window with ce_i and ne_i */
int penalization(const VI &line, int ce_i, int ne_i)
{

  int total_penalizations, requirements;
  total_penalizations = requirements = 0;

  // count inicial windows
  int uw = 0;
  for (int i = ne_i; i > 0; --i)
  {
    requirements += line[uw];
    total_penalizations += max(requirements - ce_i, 0);
    ++uw;
  }

  // count mid penalizations
  for (; uw < line.size(); ++uw)
  {
    requirements += line[uw];
    requirements -= line[uw - ne_i];
    total_penalizations += max(requirements - ce_i, 0);
  }

  // count end windows
  while (requirements > 0)
  {
    requirements -= line[uw - ne_i];
    total_penalizations += max(requirements - ce_i, 0);
    ++uw;
  }

  return total_penalizations;
}

/* Returns the improvements that are still needed in the station of the column j of improvements*/
int ImprovementsLeft(int j, const VI &cleft, const MI &improvements)
{
  int ones = 0;
  for (int i = 0; i < cleft.size(); ++i)
    ones += cleft[i] * improvements[i][j];
  return ones;
}

/* Calculates the lower bownd for a particular station. The lower bound consist in counting how many improvements are needed in
that station and distribute them in the best possible way such as the penalization is minimum */
int lower_bound_station(int i, int j, const VI &cs, const VI &cleft, VI line, int ce_i, int ne_j, const MI &imp)
{
  int ones = ImprovementsLeft(j, cleft, imp); // number of improvements (ones) left in that station
  int zeros = (cs.size() - i - 1) - ones;     // number of "not improvements" (zeros) left in that station

  int requierements = 0; // Counts how many requirements of improvements are needed in the last window
  for (int z = i; z >= 0 and z > i - ne_j; --z)
    if (line[z]) // if is needed an improvement
      ++requierements;

  while (zeros > 0 and ones > 0)
  {
    int substract = i - ne_j + 1; // postion of the vector that gets outside the window when adding a new element
    if (substract >= 0)
      requierements -= line[substract];

    // if there aren't zeros left or adding an improvement doesn't add penalization
    if (ones > 0 and (zeros == 0 or requierements + 1 <= ce_i))
    {
      line.push_back(1);
      ++requierements;
      --ones;
    }
    else // if there are no ones left or adding an improvement adds penalization
    {
      line.push_back(0);
      --zeros;
    }
    ++i;
  }

  // Returns the penalization of the line created which is the smaller one given the cs.
  return penalization(line, ce_i, ne_j);
}

/* Returns a lower bound for the current solution cs, this is the sum of the minimum penalization reachable of each station
(without taking into account the others) if the classes (hence the improvements of each class) were distributed in a way
that minimizes the penalization. The cs (current solution) vector is filled up to the position i (i included), so the lower
bound will be created from position i. */
int lower_bound(int i, const VI &cs, const VI &cleft, const VS &st, const Window &w, const MI &imp)
{
  int lb = 0;
  if (i >= 0)
  {
    const auto &[ce, ne] = w;
    int M = ne.size();

    // calculate lower bound for each station (each station is independent)
    for (int j = 0; j < M; ++j)
      lb += lower_bound_station(i, j, cs, cleft, st[j].line, ce[j], ne[j], imp);
  }
  return lb;
}

/* Updates de elements of the vector generator in order to
optimize the order in wich the exhaustive search is done*/
VI update_generator(VI generator)
{
  // random_shuffle(generator.begin(), generator.end());
  int K = generator.size();
  for (int i = 0; i < K; ++i)
    generator[i] = (generator[i] + 1) % K;

  return generator;
}

/* Restores all the changes by adding the last class_id to the stations */
void restore_stations(VS &st, const VI &ne)
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

/* Exhaustive search parameters:

1. i - current solution index: indicates the recursive call, the position i of the current solution will be filled in this call
2. cs - current solution: current solution of that recursive call
3. cp - current penalization: penalization of the current solution
4. w - Window: contains vectors ce i ne, that represent the capacity and penalization from each station
5. cl - cars left: contains how many cars of each class are still waiting to be chosen
6. st - stations: vector of stations filled up to position i, with the cars already chosen
7. imp - improvements: imp[i] ( i = 0, ..., K ) contains the improvements that class -i- requires
8. gen - generator: order in which the exhaustive search will be done, contains all classes 0,...,K

*/
void exhaustive_search_rec(int i, VI &cs, int cp, int &mp, const Window &w, VI &cl, VS &st, const MI &imp, const VI &gen)
{
  int C = cs.size();

  if (i == C) // the current solution is better than the best one yet
  {
    mp = cp;
    write_solution(cs, cp, now() - start);
  }

  if (lower_bound(i - 1, cs, cl, st, w, imp) < mp)
  {
    for (int class_id : gen) // for each class_id in the generator vector
    {
      if (cl[class_id] > 0) // check if there are still cars of that class
      {
        --cl[class_id];
        cs[i] = class_id;

        // penalzation increse when adding class_id to the stations
        int updated_penalization = UPL(imp[class_id], st, w, i + 1 == C);

        if (cp + updated_penalization < mp)
          exhaustive_search_rec(++i, cs, cp + updated_penalization, mp, w, cl, st, imp, update_generator(gen));

        // Restore all changes done
        restore_stations(st, w.ne);
        ++cl[class_id];
      }
    }
  }
}

/* Returns a vector of stations of size M, where each station is inicialized with null requirements
and an empty line, that will be filled with C elements */
VS inicialize_stations(int C, int M)
{
  VS stations(M);

  int requirements = 0; // when the line is empty there are no requirements
  VI inicial_line;
  inicial_line.reserve(C);

  for (Station &station : stations)
    station = Station{requirements, inicial_line};

  return stations;
}

/* Returns a vector of size K where each position is the correspondent index.
Example: if K = 10 ----> generator = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9} */
VI create_generator(int K)
{
  VI generator(K);                             // vector with K ints
  iota(generator.begin(), generator.end(), 0); // Fill with 0, 1, ..., K
  return generator;
}

/* Given a data and the output file name writes in the output file the solution of the exhaustive search */
void exhaustive_search(const Data &data)
{
  auto [C, M, K, w, improvements, cleft] = data; // get the data

  VS stations = inicialize_stations(C, M); // vector of empty stations that will be filled during the search
  VI generator = create_generator(K);      // Vector that inticates search space order in the exhaustive search

  VI current_solution(C); // current solution of the problem, will vary during the executions
  int mp = INT_MAX;       // inicial minimum penalization, INT_MAX works as infinity

  start = now(); // inicialize the counter
  exhaustive_search_rec(0, current_solution, 0, mp, w, cleft, stations, improvements, generator);
}

void read_input(ifstream &input_file, Data &data)
{
  auto &[C, M, K, w, improvements, cleft] = data;
  auto &[ce, ne] = w;
  input_file >> C >> M >> K;

  ce.resize(M);
  ne.resize(M);
  cleft.resize(K);
  improvements.resize(K, VI(M));

  for (auto &capacity : ce)
    input_file >> capacity;
  for (auto &window : ne)
    input_file >> window;

  for (int i = 0; i < K; ++i)
  {
    int aux;
    input_file >> aux >> cleft[i];
    for (int j = 0; j < M; ++j)
      input_file >> improvements[i][j];
  }
}

int main(int argc, char *argv[])
{
  if (argc != 3)
  {
    cout << "Syntax: " << argv[0] << " input_file output_file" << endl;
    exit(1);
  }

  ifstream input_file(argv[1]);
  output_file = argv[2];

  Data data;

  read_input(input_file, data);

  exhaustive_search(data);
}