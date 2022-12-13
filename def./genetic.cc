#include <iostream>
#include <vector>
#include <climits>
#include <chrono>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <map>
#include <random>
using namespace std;

struct Class;
struct Station;
struct Individual;

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

int penalization(const VI &P, const Data &data)
{

  MVI stations_matrix = createSM(P, data.C, data.M, data.improvements);
  int total_pen = 0;

  for (int i = 0; i < data.M; ++i)
    total_pen += count_penalization(stations_matrix[i], data.w.ce[i], data.w.ne[i]);
  return total_pen;
}

void generate_permutation(VI &parent)
{
  // shuffle(parent.begin(), parent.end(), default_random_engine(parent.size()));
  random_shuffle(parent.begin(), parent.end());
}

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

bool improved(Individual &best_individual, const Individual &P)
{
  if (P.penalization < best_individual.penalization)
  {
    best_individual = P;
    return true;
  }
  return false;
}

bool best_fit(const Individual &P1, const Individual &P2)
{
  return P1.penalization < P2.penalization;
}

void sort_population(MI &population)
{
  sort(population.begin(), population.end(), best_fit);
}

VI generate_random_cut_points(int C)
{
  random_device rd;                            // obtain a random number from hardware
  mt19937 gen(rd());                           // seed the generator
  uniform_int_distribution<> distr1(0, C - 2); // define the range
  int cut1 = distr1(gen);
  uniform_int_distribution<> distr2(cut1 + 1, C - 1);
  int cut2 = distr2(gen);
  return {cut1, cut2};
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

MVI mid_permutation_recombination(const VI &P1, const VI &P2, const VI &cuts, const VI &cleft)
{
  int left = cuts[0];
  int right = cuts[1];
  VI children1 = P1, children2 = P2;

  for (int i = left; i <= right; ++i)
  {
    children1[i] = P2[i];
    children2[i] = P1[i];
  }

  adapt_solution(children1, cleft);
  adapt_solution(children2, cleft);

  return {children1, children2};
}

MVI edge_permutation_recombination(const VI &P1, const VI &P2, const VI &cuts, const VI &cleft)
{
  int C = P1.size();
  int left = cuts[0];
  int right = cuts[1];
  VI children1 = P1, children2 = P2;

  // second pair of children
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

  return {children1, children2};
}

MVI uniform_recombination(const VI &P1, const VI &P2, const VI &cleft)
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

  return {children1, children2};
}

MVI partially_mapped_recombination(const VI &P1, const VI &P2, const VI &cuts, const VI &cleft)
{
  int C = P1.size();
  int left = cuts[0];
  int right = cuts[1];
  map<int, int> P1map, P2map, permutations;
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

  for (int i = left; i <= right; ++i)
  {
    int k = c1[i];
    c1[i] = c2[i];
    c2[i] = k;
    permutations[c1[i]] = c2[i];
    permutations[c2[i]] = c1[i];
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

  return {children1, children2};
}

MVI recombinate(const VI &P1, const VI &P2, const VI &cleft)
{
  int C = P1.size();

  MVI c1 = mid_permutation_recombination(P1, P2, generate_random_cut_points(C), cleft);
  MVI c2 = edge_permutation_recombination(P2, P2, generate_random_cut_points(C), cleft);
  MVI c3 = partially_mapped_recombination(P1, P2, generate_random_cut_points(C), cleft);
  MVI c4 = uniform_recombination(P1, P2, cleft);

  return {c1[0], c1[1], c2[0], c2[1], c3[0], c3[1], c4[0], c4[1]};
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

void update_population(MI &parents, const MI &children, int numChildren, int &psum, const bool &add)
{

  if (add)
  {
    for (int i = 0; i < numChildren; ++i)
    {
      parents.push_back(children[i]);
      psum += children[i].penalization;
    }
  }
  else
  {
    for (int i = 0; i < numChildren; ++i)
    {
      psum -= parents.back().penalization;
      parents.pop_back();
    }
  }
}

MI generate_inicial_population(const Data &data, int &psum)
{
  VI Generator = create_generator_individual(data.C, data.cleft);

  // select the population size (it will remain constant during the iterations)
  int psize = 1000;

  MI inicial_population(psize);

  for (int i = 0; i < psize; ++i)
  {
    generate_permutation(Generator);
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

void refresh_population(MI &population, int &psum, const Data &data)
{
  int psize = population.size();

  for (int i = 0; i < psize - 1; ++i)
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

// probability of accepting a worsening move
double probability(double T, int diff)
{
  return 1 / exp(diff / T);
}

void simulated_annealing(VI solution, int pen, const Data &data, const double &start, const string &out)
{
  int C = data.C;
  double T = 2;
  VI bs = solution;
  int bp = pen;
  while (now() - start < 59)
  {
    VI sp = solution;
    swap(sp[rand() % C], sp[rand() % C]);
    int cp = penalization(sp, data);

    if (cp < pen)
    {
      solution = sp;
      pen = cp;
      if (pen < bp)
      {
        bs = solution;
        bp = pen;
        write_solution(bs, bp, now() - start, out);
      }
    }
    else
    {
      double rndNumber = rand() / (double)RAND_MAX;
      if (rndNumber < probability(T, cp - pen))
      {
        solution = sp;
        pen = cp;
      }
    }
    T = T * 0.9;
  }
}

void genetic(const string &out, const Data &data)
{
  double start = now(); // inicialize the counter

  int individual_with_less_penalization = 0;
  int psum = UNDEF;                                        // total sum of penalization of the enitre population
  MI population = generate_inicial_population(data, psum); // random generated vector of parents (inicial population)
  int numChildren = 8;
  double mutation_probability = 0;

  update_fitness(population, psum);
  sort_population(population);

  Individual best_individual = population[individual_with_less_penalization];

  // termination conditions: we stop if in the last -termination_condition- generations
  // it hasn't been an iMIrovement in the solution
  int termination_conditions = 100000;
  int tc = termination_conditions;

  // do genetic algorithm while termination conditions aren't accoMIlished or until the
  // time has expired
  while (tc > 0 and (now() - start) < 50)
  {

    int pick1 = roulette_wheel_selection(population);
    int pick2 = roulette_wheel_selection(population);

    // selection of the two parents to execute the recombination
    Individual Parent1 = population[pick1];
    Individual Parent2 = population[check(pick1, pick2)]; // we want to make sure to select differt parents

    // cout << pick1 << " " << pick2 << endl;
    //  -numChildren- children are generated by recombination from the previous parents
    MI children = recombination(Parent1.solution, Parent2.solution, data);

    /*
    for (auto &d : Parent1.solution)
      cout << d << " ";
    cout << "pen: " << Parent1.penalization;
    cout << endl;
    for (auto &d : Parent2.solution)
      cout << d << " ";
    cout << "pen: " << Parent2.penalization;
    cout << endl;
    cout << "-------------------" << endl;
    for (auto &m : children)
    {
      for (auto &d : m.solution)
        cout << d << " ";
      cout << "pen: " << m.penalization << endl;
      cout << endl;
    }
    cout << endl;
    */

    // mutation probability increases during the generations in order to avoid reaching a local minimum
    mutation_probability += 0.001;

    // in some cases, a mutation appears in the children in order to allow the appereance of new traits
    mutation(children, mutation_probability, data);

    // the children are introduced into the population
    update_population(population, children, numChildren, psum, true);

    // the full population is sorted by fitness, best individuals appear fisrt
    sort_population(population);

    // the worst -NumChildren- individuals are expelled of the population
    update_population(population, {}, numChildren, psum, false);

    // when we notice that the solution is not improving, it may be the case we have reachead a local minimum
    // so we add some new individuals and expell some of the current population
    if (tc < termination_conditions / 2)
      refresh_population(population, psum, data);

    // the fitness of each individual is updated
    update_fitness(population, psum);

    // if the best individual in our current population is better than the best individual ever found,
    // we update this last and restart the termination conditions. On the other hand, if there hasn't been
    // an improvement in this generation we decrease the termination conditions
    if (improved(best_individual, population[individual_with_less_penalization]))
    {
      write_solution(best_individual.solution, best_individual.penalization, now() - start, out);
      tc = termination_conditions; // restart the termination conditions
      mutation_probability = 0;
    }
    else
      --tc;
  }

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