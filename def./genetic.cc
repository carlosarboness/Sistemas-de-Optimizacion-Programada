#include <iostream>
#include <vector>
#include <climits>
#include <chrono>
#include <cassert>
#include <fstream>
#include <algorithm>
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
    shuffle(parent.begin(), parent.end(), default_random_engine(parent.size()));
    // random_shuffle(parent.begin(), parent.end());
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

void adapt_solution(VI &child, VI &cleft, int left, int right)
{
    int C = child.size();
    int j = 0;

    while (j < C)
    {
        if (j == left)
        {
            j = right + 1;
            if (j == C)
                break;
        }
        if (int a = child[j]; cleft[a] > 0)
            --cleft[a];
        else
        {
            child[j] = most_cleft_index(cleft);
            --cleft[child[j]];
        }
        ++j;
    }
}

void adapt_solution2(VI &child, VI &cleft, int left, int right)
{
    for (int j = left; j <= right; ++j)
    {
        if (int a = child[j]; cleft[a] > 0)
            --cleft[a];
        else
        {
            child[j] = most_cleft_index(cleft);
            --cleft[child[j]];
        }
    }
}

VI uniform(const VI &P1, const VI &P2, VI cleft)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distr(1, 2);

    VI child(P1.size());

    for (int i = 0; i < P1.size(); ++i)
    {
        if (distr(gen) == 0)
        {
            int a = P1[i];
            if (cleft[a] > 0)
            {
                child[i] = a;
                --cleft[a];
            }
            else if (cleft[P2[i]] > 0)
            {
                child[i] = P2[i];
                --cleft[P2[i]];
            }
            else
            {
                int b = most_cleft_index(cleft);
                child[i] = b;
                --cleft[b];
            }
        }
        else
        {
            int a = P2[i];
            if (cleft[a] > 0)
            {
                child[i] = a;
                --cleft[a];
            }
            else if (cleft[P1[i]] > 0)
            {
                child[i] = P1[i];
                --cleft[P1[i]];
            }
            else
            {
                int b = most_cleft_index(cleft);
                child[i] = b;
                --cleft[b];
            }
        }
    }
    return child;
}

MVI recombinate(const VI &P1, const VI &P2, int left, int right, const VI &cleft)
{
    int C = P1.size();

    VI children1, children2, children3, children4;
    children1 = children3 = P1;
    children2 = children4 = P2;
    VI cleft1, cleft2, cleft3, cleft4;
    cleft1 = cleft2 = cleft3 = cleft4 = cleft;

    // fisrt pair of children
    for (int i = left; i <= right; ++i)
    {
        children1[i] = P2[i];
        --cleft1[P2[i]];
        children2[i] = P1[i];
        --cleft2[P1[i]];
    }

    // second pair of children
    for (int i = 0; i < left; ++i)
    {
        children3[i] = P2[i];
        --cleft3[P2[i]];
        children4[i] = P1[i];
        --cleft4[P1[i]];
    }
    for (int i = right + 1; i < C; ++i)
    {
        children3[i] = P2[i];
        --cleft3[P2[i]];
        children4[i] = P1[i];
        --cleft4[P1[i]];
    }

    adapt_solution(children1, cleft1, left, right);
    adapt_solution(children2, cleft2, left, right);
    adapt_solution2(children3, cleft3, left, right);
    adapt_solution2(children4, cleft4, left, right);

    VI children5 = uniform(P1, P2, cleft);

    return {children1, children2, children3, children4, children5};
}

MI recombination(const VI &P1, const VI &P2, const Data &data)
{
    int numChildren = 5;

    VI cuts = generate_random_cut_points(data.C);
    MVI recombination = recombinate(P1, P2, cuts[0], cuts[1], data.cleft);
    MI children;
    children.reserve(numChildren);

    for (int i = 0; i < numChildren; ++i)
        children.push_back(Individual{recombination[i], penalization(recombination[i], data), UNDEFINED});

    return children;
}

/* Adds a mutation to the children. A mutation is a permutation of two elements of each children */
void mutation(MI &children, const Data &data)
{
    double rndNumber = rand() / (double)RAND_MAX; // random number in range [0, 1)

    // add a mutation to the children with propability 0.01
    if (rndNumber < 0.01)
    {

        int C = children[0].solution.size();
        int j = rand() % C;
        int k = rand() % C;

        for (int i = 0; i < children.size(); ++i)
        {
            swap(children[i].solution[j], children[i].solution[k]);
            children[i].penalization = penalization(children[i].solution, data);
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

void update_population(MI &parents, const MI &children, int &psum, const bool &add)
{
    int numChildren = 5;

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

void genetic(const string &out, const Data &data)
{
    double start = now(); // inicialize the counter

    int individual_with_less_penalization = 0;
    int psum = UNDEF;                                        // total sum of penalization of the enitre population
    MI population = generate_inicial_population(data, psum); // random generated vector of parents (inicial population)
    int numChildren = 5;                                     // number of children generated from each 2 parents

    update_fitness(population, psum);
    sort_population(population);

    Individual best_individual = population[individual_with_less_penalization];

    // termination conditions: we stop if in the last -termination_condition- generations
    // it hasn't been an iMIrovement in the solution
    int termination_conditions = 100000;
    int tc = termination_conditions;

    // do genetic algorithm while termination conditions aren't accoMIlished or until the
    // time has expired
    while (tc > 0 and (now() - start) < 60)
    {

        int pick1 = roulette_wheel_selection(population);
        int pick2 = roulette_wheel_selection(population);

        // selection of the two parents to execute the recombination
        Individual Parent1 = population[pick1];
        Individual Parent2 = population[check(pick1, pick2)]; // we want to make sure to select differt parents

        // -numChildren- children are generated by recombination from the previous parents
        MI children = recombination(Parent1.solution, Parent2.solution, data);

        // in some cases, a mutation appears in the children in order to allow the appereance of new traits
        mutation(children, data);

        // the children are introduced into the population
        update_population(population, children, psum, true);

        // the full population is sorted by fitness, best individuals appear fisrt
        sort_population(population);

        // the worst -NumChildren- individuals are expelled of the population
        update_population(population, {}, psum, false);

        // the fitness of each individual is updated
        update_fitness(population, psum);

        // if the best individual in our current population is better than the best individual ever found,
        // we update this last and restart the termination conditions. On the other hand, if there hasn't been
        // an improvement in this generation we decrease the termination conditions
        if (improved(best_individual, population[individual_with_less_penalization]))
        {
            write_solution(best_individual.solution, best_individual.penalization, now() - start, out);
            tc = termination_conditions; // restart the termination conditions
        }
        else
            --tc;
    }
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
