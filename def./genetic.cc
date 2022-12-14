#include <iostream>
#include <vector>
#include <climits>
#include <chrono>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <map>
#include <random>
#include <set>
using namespace std;

struct Class;
struct Station;
struct Individual;
using offspring = vector<vector<int>>;

using VI = vector<int>;
using MVI = vector<VI>;
using MI = vector<Individual>;
using is = ifstream;
using os = ofstream;

int UNDEF = -1;
double UNDEFINED = -1;

struct Window
{
  VI ce, ne;
};

// Stores all the information given by the input
struct Data
{
  int C, M, K;
  Window w;
  MVI improvements;
  VI cleft; // vector with the total cars of each class (cars left)
};

struct Individual
{
  VI solution;
  int penalization;
  double fitness;
};

struct crossover
{
  map<int, int> P1map, P2map;
  VI c1, c2;
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

MVI createSM(const VI &P, int C, int M, const MVI &improvements)
{

  MVI SM(M, VI(C));
  for (int i = 0; i < C; ++i)
    for (int j = 0; j < M; ++j)
      SM[j][i] = improvements[P[i]][j];
  return SM;
}
/* Counts the penalization of the solution in P */
int penalization(const VI &P, const Data &data)
{
  // generate the matrix with the stations representation, to make easier the counting
  MVI stations_matrix = createSM(P, data.C, data.M, data.improvements);
  int total_pen = 0;

  // for each station count its penalziation
  for (int i = 0; i < data.M; ++i)
    total_pen += count_penalization(stations_matrix[i], data.w.ce[i], data.w.ne[i]);

  return total_pen;
}
/* permutates randomly the elements of the vector -individual- */
void generate_permutation(VI &individual)
{
  // shuffle(individual.begin(), individual.end(), default_random_engine(individual.size()));
  random_shuffle(individual.begin(), individual.end());
}
/* Returns a vector of size C in which the elements are the classes that the final solution will contain
the elements of the final solution in order */
VI create_generator_individual(int C, const VI &cleft)
{
  VI individual;
  individual.reserve(C);
  for (int i = 0; i < cleft.size(); ++i)
  {
    int cl = cleft[i];
    for (; cl > 0; --cl)
      individual.push_back(i);
  }
  return individual;
}

/* Returns true if the current_best_individual is better than the best_individual ever found,
else returns false. If returns true it also updates the best_individual*/
bool improved(Individual &best_individual, const Individual &current_best_individual)
{
  if (current_best_individual.fitness < best_individual.fitness)
  {
    best_individual = current_best_individual;
    return true;
  }
  return false;
}

/* Auxiliar comparating function for -sort_population()-. Individuals with less penalzaitions go first */
bool lower_penalization_first(const Individual &P1, const Individual &P2)
{
  return P1.penalization < P2.penalization;
}

/* Input is a matrix of indivuals. We sort this individuals by -lower_penalization_first- criteria defined above */
void sort_population(MI &population)
{
  sort(population.begin(), population.end(), lower_penalization_first);
}

/* Returns a pair of two integers corresponding to two random cuts in range [0, C).
The first cut is always lower that the second one (by construction of the distributions. */
pair<int, int> generate_random_cut_points(int C)
{
  random_device rng;  // obtain a random number from hardware
  mt19937 gen(rng()); // seed the generator

  uniform_int_distribution<int> distribution1(0, C - 2); // uniform distribution in range[0, C-2]
  int first_cut_point = distribution1(gen);              // get a number in that range with uniform probability

  // uniform distribution at the right part of the first cut
  uniform_int_distribution<int> distribution2(first_cut_point + 1, C - 1);
  int second_cut_point = distribution2(gen); // get a number in that range with uniform probability

  // create the pair of cuts
  pair<int, int> cuts(first_cut_point, second_cut_point);

  return cuts;
}

int most_cleft_index(const VI &cleft)
{
  auto it = max_element(cleft.begin(), cleft.end());
  return distance(cleft.begin(), it);
}

void adapt_solution(VI &child, VI cleft)
{
  int C = child.size();

  for (int i = 0; i < C; ++i)
  {
    if (cleft[child[i]] > 0)
      --cleft[child[i]];
    else
    {
      child[i] = most_cleft_index(cleft);
      --cleft[child[i]];
    }
  }
}

void mid_permutation_recombination(const VI &P1, const VI &P2, const pair<int, int> &cuts, const VI &cleft, offspring &children)
{
  int left = cuts.first;
  int right = cuts.second;
  VI children1 = P1, children2 = P2;

  for (int i = left; i <= right; ++i)
  {
    children1[i] = P2[i];
    children2[i] = P1[i];
  }

  adapt_solution(children1, cleft);
  adapt_solution(children2, cleft);

  children.push_back(children1);
  children.push_back(children2);
}

void edge_permutation_recombination(const VI &P1, const VI &P2, const pair<int, int> &cuts, const VI &cleft, offspring &children)
{
  int C = P1.size();
  int left = cuts.first;
  int right = cuts.second;
  VI children1 = P1, children2 = P2;

  for (int i = 0; i < left; ++i)
  {
    children1[i] = P2[i];
    children2[i] = P1[i];
  }
  for (int i = right + 1; i < C; ++i)
  {
    children1[i] = P2[i];
    children2[i] = P1[i];
  }

  adapt_solution(children1, cleft);
  adapt_solution(children2, cleft);

  children.push_back(children1);
  children.push_back(children2);
}

void uniform_recombination(const VI &P1, const VI &P2, const VI &cleft, offspring &children)
{
  int C = P1.size();
  VI children1(C), children2(C);

  for (int i = 0; i < C; ++i)
  {
    if (rand() / (double)RAND_MAX < 0.5)
    {
      children1[i] = P1[i];
      children2[i] = P2[i];
    }
    else
    {
      children1[i] = P2[i];
      children2[i] = P1[i];
    }
  }

  adapt_solution(children1, cleft);
  adapt_solution(children2, cleft);

  children.push_back(children1);
  children.push_back(children2);
}

crossover generate_inicial_data(const VI &P1, const VI &P2, int left, int right)
{
  int C = P1.size();
  map<int, int> P1map, P2map;
  vector<bool> used(C, false);
  VI c1(C, -1), c2(C, -1);

  for (int i = 0; i < C; ++i)
  {
    c1[i] = i;
    P1map[i] = P1[i];
  }

  for (int i = 0; i < C; ++i)
  {
    for (int j = 0; j < C; ++j)
    {
      if (not used[j] and P2[i] == P1[j])
      {
        used[j] = true;
        c2[i] = j;
        P2map[j] = P2[i];
        break;
      }
    }
  }

  return crossover{P1map, P2map, c1, c2};
}

void partially_mapped_recombination(const VI &P1, const VI &P2, const pair<int, int> &cuts, const VI &cleft, offspring &children)
{
  int C = P1.size();
  int left = cuts.first;
  int right = cuts.second;
  auto [P1map, P2map, c1, c2] = generate_inicial_data(P1, P2, left, right);
  map<int, int> permutations;

  for (int i = left; i <= right; ++i)
  {
    permutations[c1[i]] = c2[i];
    permutations[c2[i]] = c1[i];
    swap(c1[i], c2[i]);
  }

  for (int i = 0; i < left; ++i)
  {
    if (permutations.count(c1[i]))
      c1[i] = permutations[c1[i]];
    if (permutations.count(c2[i]))
      c2[i] = permutations[c2[i]];
  }

  for (int i = right + 1; i < C; ++i)
  {
    if (permutations.count(c1[i]))
      c1[i] = permutations[c1[i]];
    if (permutations.count(c2[i]))
      c2[i] = permutations[c2[i]];
  }
  VI children1(C), children2(C);
  for (int i = 0; i < C; ++i)
  {
    children1[i] = P1map[c1[i]];
    children2[i] = P2map[c2[i]];
  }

  adapt_solution(children1, cleft);
  adapt_solution(children2, cleft);

  children.push_back(children1);
  children.push_back(children2);
}

set<int> generate_set(int C)
{
  set<int> S;
  for (int i = 0; i < C; ++i)
    S.insert(i);
  return S;
}

void order_recombination(const VI &P1, const VI &P2, const pair<int, int> &cuts, offspring &children)
{
  int C = P1.size();
  int left = cuts.first;
  int right = cuts.second;
  auto [P1map, P2map, c1, c2] = generate_inicial_data(P1, P2, left, right);

  set<int> s1, s2;
  s1 = s2 = generate_set(C);
  VI order1, order2;
  order1.reserve(C - (right - left + 1));
  order2.reserve(C - (right - left + 1));

  for (int i = left; i <= right; ++i)
  {
    s1.erase(c1[i]);
    s2.erase(c2[i]);
  }

  for (int i = 0; i < C; ++i)
  {
    if (s1.count(c2[i]))
      order1.push_back(c2[i]);

    if (s2.count(c1[i]))
      order2.push_back(c1[i]);
  }

  int i = 0;
  for (int j = 0; j < order1.size(); ++j)
  {
    if (i == left)
      i = right + 1;
    c1[i] = order1[j];
    c2[i] = order2[j];
    ++i;
  }

  VI children1(C), children2(C);
  for (int i = 0; i < C; ++i)
  {
    children1[i] = P1map[c1[i]];
    children2[i] = P2map[c2[i]];
  }

  children.push_back(children1);
  children.push_back(children2);
}

offspring recombinate(const VI &P1, const VI &P2, const VI &cleft)
{
  int C = P1.size();
  int numChildren = 10;
  offspring children;
  children.reserve(numChildren);

  mid_permutation_recombination(P1, P2, generate_random_cut_points(C), cleft, children);
  edge_permutation_recombination(P1, P2, generate_random_cut_points(C), cleft, children);
  partially_mapped_recombination(P1, P2, generate_random_cut_points(C), cleft, children);
  order_recombination(P1, P2, generate_random_cut_points(C), children);
  uniform_recombination(P1, P2, cleft, children);

  return children;
}

MI recombination(const VI &P1, const VI &P2, const Data &data)
{

  MVI recombination = recombinate(P1, P2, data.cleft);
  MI children;
  children.reserve(recombination.size());

  for (const VI &offspring : recombination)
    children.push_back(Individual{offspring, penalization(offspring, data), UNDEFINED});

  return children;
}

/* Adds a mutation to the children. A mutation is a permutation of two elements of each children */
void mutation(MI &children, double probability, const Data &data)
{
  int C = data.C;
  double rndNumber = rand() / (double)RAND_MAX; // random number in range [0, 1)

  // add a mutation to the children with propability 0.01
  if (rndNumber < probability)
  {
    for (Individual &individual : children)
    {
      swap(individual.solution[rand() % C], individual.solution[rand() % C]);
      individual.penalization = penalization(individual.solution, data);
    }
  }
}

// updates the fitness of every parent
void update_fitness(MI &population, const int &psum)
{
  int n = population.size();

  for (Individual &individual : population)
    individual.fitness = (double)(psum - individual.penalization) / (psum * (n - 1));
}
/* Returns a the index of the parent chosen by the roulette wheel selection method.
Parents with higher fitness have more probability to be chosen */
int roulette_wheel_selection(const MI &parents)
{
  /*
                       | random number
                       |
                       X
    +--------------+-----+----+--+--++---+-----+
    |              |     |    |  |  ||   |     |
    +--------------+-----+----+--+--++---+-----+
   0.0            p0   p0+p1                  1.0 == sum(p[i])
  */

  double rndNumber = rand() / (double)RAND_MAX; // random number in range [0, 1)
  double offset = 0.0;                          // acumulated probability
  int pick = 0;                                 // roulette winner

  for (int i = 0; i < parents.size(); i++)
  {
    offset += parents[i].fitness;
    if (rndNumber < offset)
    {
      pick = i;
      break;
    }
  }
  return pick;
}

void update_population(MI &population, const MI &children, int &psum)
{
  int population_size = population.size();
  int numChildren = children.size();

  for (int i = 0; i < numChildren; ++i)
  {
    psum -= population.back().penalization;
    population.pop_back();
  }

  population.reserve(population_size);

  for (int i = 0; i < numChildren; ++i)
  {
    population.push_back(children[i]);
    psum += children[i].penalization;
  }
}

/* Returns a vector of individuals, that is the inicial random generated population, int he variable psum
we save the sum of all the individuals penalization, in order to save computations later */
MI generate_inicial_population(const Data &data, int &psum)
{
  // vector with all the cars classes in order (p.e [0, 0, 0, 1, 1, 1, 2, 2, 2, 2])
  VI Generator = create_generator_individual(data.C, data.cleft);

  // select the population size (it will remain constant during the iterations)
  int psize = 1000;

  MI inicial_population(psize);

  for (int i = 0; i < psize; ++i)
  {
    generate_permutation(Generator); // random permutation of the vector generator, that is a new individual
    inicial_population[i] = Individual{Generator, penalization(Generator, data), UNDEFINED};
    psum += inicial_population[i].penalization;
  }

  return inicial_population;
}

int check(const int &pick1, const int &pick2)
{
  if (pick1 != pick2)
    return pick2;
  return (pick2 - 1 >= 0 ? pick2 - 1 : pick2 + 1);
}

/* A percentatge of the current population is replaced by new random generated individuals in order to trying to
escape from a local minima. The sum of the total penalizations is updated with the new individuals penalization */
void refresh_population(MI &population, int &psum, const Data &data)
{
  int psize = population.size();
  // percentatge of population replaced, only 10% of individuals of the current
  // population will remain in the next genration, the other individuals will be replaced
  double PopulationReplacedSize = 0.9;

  for (int i = 0; i < psize * PopulationReplacedSize; ++i)
  {
    psum -= population.back().penalization;
    population.pop_back();
  }

  population.reserve(psize);

  VI Generator = create_generator_individual(data.C, data.cleft);

  while (population.size() < psize)
  {
    generate_permutation(Generator);
    int pen = penalization(Generator, data);
    population.push_back(Individual{Generator, pen, UNDEFINED});
    psum += pen;
  }
}

/* Returns the probability of accepting a worsening move following the Bolltzman distribution
Precondition: diff >= 0*/
double probability(const double &T, int diff)
{
  return 1 / exp(diff / T);
}

/* Metaheuristics applied to solution to try to find a better one. In each iteration we visit a neighbour
of our current solution. If this neighbour is a better one, we keep it, if not we also keep with some
probability that will decrease during the iterations. Each time a better solution is found is  written in
the -out- file. */
void simulated_annealing(VI current_solution, int current_pen, const Data &data, const double &start, const string &out)
{
  int C = data.C;

  VI best_solution = current_solution; // saves the better solution
  int best_penalization = current_pen; // saves the lowest penalization
  double T = 100;                      // inicial temperature

  while (now() - start < 59)
  {
    cout << T << endl;
    VI neighbour = current_solution;
    swap(neighbour[rand() % C], neighbour[rand() % C]); // select the neighbour
    int pen = penalization(neighbour, data);            // count its penalization

    if (pen < current_pen)
    {
      current_solution = neighbour;
      current_pen = pen;

      if (current_pen < best_penalization)
      {
        best_solution = current_solution;
        best_penalization = current_pen;
        write_solution(best_solution, best_penalization, now() - start, out);
      }
    }
    else
    {
      double rndNumber = rand() / (double)RAND_MAX;
      if (rndNumber < probability(T, pen - current_pen))
      {
        current_solution = neighbour;
        current_pen = pen;
      }
    }

    T = T * 0.95; // update (decrease) the temperature
  }
}

void genetic(const string &out, const Data &data)
{
  double start = now(); // inicialize the counter

  int psum = UNDEF;                                        // total sum of penalization of the enitre population
  MI population = generate_inicial_population(data, psum); // random generated vector of parents (inicial population)

  update_fitness(population, psum); // updates the fitness of every individual in the population
  sort_population(population);      // sorts the population by penalization, less penalizations solutions go first

  // after sorting the individual with less penalizations will be in the postion 0 of the vector
  int individual_with_less_penalization = 0;
  Individual best_individual = population[individual_with_less_penalization]; // saves the best individual ever generated

  double refreshing_rate = 0.2;    // indicates how fast the population is refreshed
  double mutation_probability = 0; // inicial mutation probability, will increase during the generations

  // termination conditions: we stop if in the last -termination_condition- generations hasn't been an
  // improvement in the solution. This parameter can be changed depending on the statement
  int generations_without_improvements = 100000;
  int tc = generations_without_improvements;

  // do genetic algorithm while termination conditions are not met or until the selected time has expired
  while (tc > 0 and (now() - start) < 0)
  {
    // execute the roulette wheel selection method to select 2 parents
    int pick1 = roulette_wheel_selection(population);
    int pick2 = roulette_wheel_selection(population);

    // selection of the two parents to execute the recombination
    Individual Parent1 = population[pick1];
    Individual Parent2 = population[check(pick1, pick2)]; // we want to make sure to select differt parents

    //  -numChildren- children are generated by recombination from the previous parents
    MI children = recombination(Parent1.solution, Parent2.solution, data);

    // mutation probability increases during the generations in order to avoid reaching a local minimum
    mutation_probability += 0.001;

    // in some cases, a mutation appears in the children in order to allow the appereance of new traits
    mutation(children, mutation_probability, data);

    // the children are introduced into the population
    update_population(population, children, psum);

    // sorts the population by penalization, less penalizations solutions go first
    sort_population(population);

    // when we notice that the solution is not improving, it may be the case we have reachead a local minimum
    // so we add some new individuals and expell some of the current population
    if (tc < generations_without_improvements * refreshing_rate)
      refresh_population(population, psum, data);

    // updates the fitness of every individual in the population
    update_fitness(population, psum);

    // if the best individual in our current population is better than the best individual ever found,
    // we update this last and restart the termination conditions. On the other hand, if there hasn't been
    // an improvement in this generation we decrease the termination conditions
    if (improved(best_individual, population[individual_with_less_penalization]))
    {
      write_solution(best_individual.solution, best_individual.penalization, now() - start, out);
      tc = generations_without_improvements; // restart the termination conditions
      mutation_probability = 0;              // restart the mutation probability
    }
    else
      --tc;
  }
  // After finishing the genetic algorithm we apply a simulated annealing heuristics to the best individual
  // ever found to try to improve the solution slighlty
  simulated_annealing(best_individual.solution, best_individual.penalization, data, start, out);

  /*
  for (auto &m : population)
  {
    for (auto &d : m.solution)
      cout << d << " ";
    cout << "pen: " << m.penalization << endl;
    cout << endl;
  }
  cout << endl
       << population.size() << endl;

  */
}

void read_input(is &in, Data &data)
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

  ifstream in(argv[1]);

  Data data;
  read_input(in, data);

  genetic(argv[2], data);
}