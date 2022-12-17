#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <random>
using namespace std;

using VI = vector<int>;
using MI = vector<vector<int>>;

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
  VI cleft;
};

double now()
{
  return clock() / double(CLOCKS_PER_SEC);
}

ofstream open(const string &out)
{
  ofstream f;
  f.open(out, ofstream::out | ofstream::trunc);
  f.setf(ios::fixed);
  f.precision(1);
  return f;
}

void write_solution(const VI &cs, int cp, const double &elapsed_time, const string &out)
{
  ofstream f = open(out);
  f << cp << " " << elapsed_time << endl;

  f << cs[0];
  for (int i = 1; i < cs.size(); ++i)
    f << " " << cs[i];

  f << endl;

  f.close();
}

/* Returns the penalization produced in the station -station- */
int count_penalization(const VI &station, int ce, int ne)
{

  int total_penalizations, requirements;
  total_penalizations = requirements = 0;

  // count inicial windows
  int uw = 0;
  for (int i = ne; i > 0; --i)
  {
    requirements += station[uw];
    total_penalizations += max(requirements - ce, 0);
    ++uw;
  }

  // count mid penalizations
  for (; uw < station.size(); ++uw)
  {
    requirements += station[uw];
    requirements -= station[uw - ne];
    total_penalizations += max(requirements - ce, 0);
  }

  // count end windows
  while (requirements > 0)
  {
    requirements -= station[uw - ne];
    total_penalizations += max(requirements - ce, 0);
    ++uw;
  }

  return total_penalizations;
}

/* Generates a stations representation matrix. Each row is a different station, with values
{0, 1} in each position, depending if the there is an improvement in that position or not.
Each column is a different class vector of improvements, depending on the -solution- vector */
MI generate_stations(const VI &solution, const MI &improvements)
{

  int C = solution.size();
  int M = improvements[0].size();
  MI stations(M, VI(C));
  for (int i = 0; i < C; ++i)
    for (int j = 0; j < M; ++j)
      stations[j][i] = improvements[solution[i]][j];
  return stations;
}

/* Counts the penalization of the solution in P */
int penalization(const VI &solution, const Data &data)
{
  MI stations = generate_stations(solution, data.improvements);

  // for each station count its penalziation
  int total_penalization = 0;
  for (int i = 0; i < data.M; ++i)
    total_penalization += count_penalization(stations[i], data.w.ce[i], data.w.ne[i]);

  return total_penalization;
}

/* Selects a random neighbour from the solution -current_solution-. The defined neighbourhood
is all the solutions that can be reachead by permutating 2 elements of -current_solution- */
VI select_neighbour(const VI &current_solution)
{
  int C = current_solution.size();
  VI neighbour = current_solution;
  swap(neighbour[rand() % C], neighbour[rand() % C]);
  return neighbour;
}

/* Simulated Annealing algorithm consists of following steps:

Step 1: Initialize – Start with a random initial solution and set initial parameters t0, a.

Step 2: Move – Select a random neighbour from the current solution.

Step 3: Calculate score – calculate the change in performance due to the move made.

Step 4: Choose – Depending on the change in performance, accept or reject the move. The prob of acceptance depending on the current “temperature”.

Step 5: Update and repeat – Update the temperature value by lowering it by a factor a.

Go back to Step 2.

The process is done until “Freezing Point” is reached.

Acceptance criteria:

 - If the performance of neighboring state is better than current state (delta >= 0) - Accept
 - If the performance of neighboring state is worse than current state (delta < 0) - Accept with probability e^{-delta/t} where, t is current temperature and delta is the performance difference between the current state and the neighboring state

 */

void run_algorithm(VI current_solution, int current_pen, double t, const double &a, const Data &data, const double &start, const string &output_file)
{

  VI best_solution = current_solution; // saves the better solution
  int best_penalization = current_pen; // saves the better solution penalization

  while (now() - start < 60)
  {
    //  select a neighbour from the current solution
    VI neighbour = select_neighbour(current_solution);
    int neighbour_pen = penalization(neighbour, data); // calculate the penalization of the neighbour
    int delta = current_pen - neighbour_pen;           // difference of performance

    if (delta >= 0) // the neighbour solution is equal or better than the current one
    {
      current_solution = neighbour;
      current_pen = neighbour_pen;

      if (current_pen < best_penalization)
      {
        best_solution = current_solution;
        best_penalization = current_pen;
        write_solution(best_solution, best_penalization, now() - start, output_file);
      }
    }
    else // the neighbour solution is worse than the current one
    {
      double rndNumber = rand() / (double)RAND_MAX; // random number in range [0, 1)
      if (rndNumber < exp(delta / t))               // acceptance of a worsening move probability
      {
        current_solution = neighbour;
        current_pen = neighbour_pen;
      }
    }

    t = t * a; // update (decrease) the temperature
  }
}

/* Generates an initial random solution for the problem. First, all the classes defined by the problem are put in order, with their respective multiplicity. For example : s = {0, 0, 1, 1, 1, 2, 2, 2, 2, 2}. Then, in order to have a better distributed initial solution a random permutation is done to vector.
For example : s = {2, 0, 1, 1, 0 ,2, 2, 2, 1, 2} */
VI initial_solution(const Data &data)
{
  VI s;
  s.reserve(data.C);

  for (int class_id = 0; class_id < data.K; ++class_id) // for each class
    for (int m = 0; m < data.cleft[class_id]; ++m)      // for each multiplicity
      s.push_back(class_id);                            // push_back of the class

  random_shuffle(s.begin(), s.end()); // random permutation of the elements in s
  return s;
}

/* Annealing parameters:

1. s - initial solution: An initial solution to the problem is created in order to start the algorithm.

2. p - cost, penalization, performance ... : Evaluation of the solution to compare it with others.

3. t0 - initial temperature: Due to the way the probability is calculated, when the temperature is higher, is it more likely that the algorithm accepts a worse solution. This promotes Exploration of the search space and allows the algorithm to more likely travel down a sub-optimal path to potentially find a global maximum.

4. a - geometric parameter: Factor by which temperature is scaled. Lower values of a restrict the search space at a faster rate than higher values.

*/
void simulated_annealing(const Data &data, const string &output_file)
{
  VI s = initial_solution(data); // create random initial solution
  int p = penalization(s, data); // calculate its penalization

  double t0 = 100; // set inicial temperature
  double a = 0.95; // temperature reduction will follow a geometric rule: t = t * a

  double start = now(); // inicialize counter

  run_algorithm(s, p, t0, a, data, start, output_file);
}

void read_input(ifstream &in, Data &data)
{
  auto &[C, M, K, w, improvements, cleft] = data;
  auto &[ce, ne] = w;
  in >> C >> M >> K;

  ce.resize(M);
  ne.resize(M);
  cleft.resize(K);
  improvements.resize(K, VI(M));

  for (auto &capacity : ce)
    in >> capacity;
  for (auto &window : ne)
    in >> window;

  int aux;
  for (int i = 0; i < K; ++i)
  {
    in >> aux >> cleft[i];
    for (int j = 0; j < M; ++j)
      in >> improvements[i][j];
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

  Data data;
  read_input(input_file, data);

  simulated_annealing(data, argv[2]);
}