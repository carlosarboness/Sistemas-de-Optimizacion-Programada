#include <iostream>
#include <vector>
#include <climits>
#include <chrono>
#include <cassert>
#include <fstream>
#include <algorithm>
using namespace std;

struct Station;
using VI = vector<int>;
using MI = vector<VI>;
using VS = vector<Station>;
using VD = vector<double>;
using MD = vector<VD>;
using is = ifstream;
using os = ofstream;

int UNDEF = -1;

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
};

// Stores the number of needed requirements in a particular
struct Station
{
    int requirements; // number of requirements in window of size ne in a particular moment of the station
    VI line;          // stores the requirements in order that need to be satisfied by the station in a particular moment
};

struct fit
{
    int id_car = UNDEF, cl = UNDEF;
    double cost = UNDEF;
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

/////// Common funcitons of both greedies ////////
/* Updates each individual line and retuns the penalitzations of it */
int add_car(int bit, Station &st, int ce_i, int ne_i, bool end)
{

    // update the window
    int uw = st.line.size() - ne_i;
    if (uw >= 0)
        st.requirements -= st.line[uw];

    st.line.push_back(bit);
    st.requirements += bit;

    int penalitzation = max(st.requirements - ce_i, 0);
    // we add to the penalitzarions the final windows
    if (end)
    {
        ++uw;
        for (; st.requirements > 0; ++uw)
        {
            st.requirements -= st.line[uw];
            penalitzation += max(st.requirements - ce_i, 0);
        }
    }
    return penalitzation;
}

int update_station_cs(const VI &imp, VS &stations, const Window &w, bool end)
{
    const auto &[ce, ne] = w;

    int M = ce.size();
    int total_penalization = 0;
    for (int i = 0; i < M; ++i)
        total_penalization += add_car(imp[i], stations[i], ce[i], ne[i], end);
    return total_penalization;
}

/* Returns a vector with the inicialized stations of the algorithm.
Firstly the requiremnts are equal to zero and the line is empty */
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

/* Returns the approximate cost (penalization) of improving a car of a class
with improvements imp_i if before or after has improvements imp_j.
The approximate cost is the sum the ne / ce of the stations that both classes
of cars need the improvements*/
double calculate_cost(const Window &w, const VI &imp_i, const VI &imp_j)
{
    const auto &[ce, ne] = w;

    double cost = 0;
    for (int k = 0; k < ce.size(); ++k)
        if (imp_i[k] == 1 and imp_j[k] == 1)
            cost += ne[k] / ce[k];
    return cost;
}

/* Returns the argmax of the VI cleft */
int most_cleft_index(const VI &cleft)
{
    auto it = max_element(cleft.begin(), cleft.end());
    return distance(cleft.begin(), it);
}
//////////////////////////////////////////////////

//////////////////// Greedy 1 ////////////////////
/*
Idea: improve the cars of the class with more cars left as soon as possible
Algorithm: while there are cars left of a particular class, select the available car
        that has more cars left. If there is a tie, choose the one with less aproximate
        cost (penalization). If there is a tie, pick the one with the lowest index.
*/
/* Returns a matrix with the approximate cost (penalization) of improving a
car of the j-th class if the last car improved was from the i-th class. */
MD generate_costs_matrix(const Data &data)
{
    const auto &[C, M, K, w, improvements] = data;

    MD costs(K, VD(K));
    for (int i = 0; i < K; ++i)
        for (int j = i; j < K; ++j)
        {
            double c = calculate_cost(w, improvements[i], improvements[j]);
            costs[i][j] = costs[j][i] = c;
        }
    return costs;
}

/* Returns the best class of cars to be improved next following the greedy algorithm given
the vector of costsof the last car of the solution and the cleft of all the classes of cars */
int best_car(VI &cleft, const VD &costs)
{
    fit bf;

    for (int i = 0; i < cleft.size(); ++i)
        if (cleft[i] > 0 and (cleft[i] > bf.cl or (cleft[i] == bf.cl and costs[i] < bf.cost)))
            bf = fit{i, cleft[i], costs[i]};

    --cleft[bf.id_car];
    return bf.id_car;
}

/* Given the data of the problem, the solution using the greedy algorithm 1 is returned
and the penalization of the solution is stored in pen */
VI gen_sol_greedy1(const Data &data, VI &cleft, int &pen)
{
    const auto &[C, M, K, w, improvements] = data;

    VS stations = inicialize_stations(C, M);
    MD costs = generate_costs_matrix(data);

    VI solution(C);
    // The first car to be improved is the one of the class with mores cars left
    solution[0] = most_cleft_index(cleft);
    --cleft[solution[0]];
    pen += update_station_cs(improvements[solution[0]], stations, w, false);

    for (int k = 1; k < C; ++k)
    {
        solution[k] = best_car(cleft, costs[solution[k - 1]]);
        pen += update_station_cs(improvements[solution[k]], stations, w, k == C - 1);
    }

    return solution;
}
//////////////////////////////////////////////////

//////////////////// Greedy 2 ////////////////////
/*
Idea: improve the cars of the class better fitness as soon as possible
Algorithm: while there are cars left of a particular class, select the available car
        that has better fitness. If there is a tie, choose the one with mone cars
        left of its class. If there is a tie, pick the one with the lowest index.
*/
// returns the argmax of V
int argmax(const VI &vi, const VD &vd)
{
    if (vi.size() != 0)
    {
        auto it = max_element(vi.begin(), vi.end());
        return distance(vi.begin(), it);
    }
    else
    {
        auto it = max_element(vd.begin(), vd.end());
        return distance(vd.begin(), it);
    }
}

/* Given the data of the problem and the cars left, returns the fitness matrix.
The fitness matrix has in the i-th row, the fitness from improving a car of class j
given that the last car was form the class i.
The fitness is a ratio dividing the cars left of the class j and the approximate cost
(penalization) of improving a car of the class j if the last car improved was from the
class i. */
MD generate_fitness_matrix(const Data &data, const VI &cleft)
{
    const auto &[C, M, K, w, improvements] = data;

    MD fitness(K, VD(K));
    for (int i = 0; i < K; ++i)
        for (int j = 0; j < K; ++j)
        {
            double c = calculate_cost(w, improvements[i], improvements[j]);
            fitness[i][j] = (double)cleft[j] / (c == 0 ? 1 : c); // els que queden dividit pel "cost"
        }
    return fitness;
}

/* Given the id_car of the last class of cars to be improved, the fitness matix and
cleft is updated. The ratio of the column id_car is updated with the new cleft.
Prec: cleft[id_car] >= 1 */
void update_fitness(int id_car, MD &fitness, VI &cleft)
{
    int K = cleft.size();
    for (int i = 0; i < K; ++i)
    {
        double c = cleft[id_car] / fitness[i][id_car];
        fitness[i][id_car] = (cleft[id_car] - 1) / c;
    }
    --cleft[id_car];
}

/* Creates the vector of the solution with the first class of cars to be improved and
updates the fitness matrix and cleft. This class is the one with more cars left. */
VI base_case(int C, VI &cleft, MD &fitness)
{
    VI solution(C);
    int base = argmax(cleft, {});
    solution[0] = base;
    update_fitness(base, fitness, cleft);
    return solution;
}

/* Given a vector of fitness and cars left, returs the id of the class of cars
following the greedy algorithm 2 */
int best(const VD &f, VI &cleft)
{
    int id = UNDEF;
    double fit = UNDEF;
    for (int i = 0; i < cleft.size(); ++i)
        if (f[i] > fit or (f[i] == fit and cleft[i] > cleft[id]))
        {
            id = i;
            fit = f[i];
        }

    return id;
}

/* Returns the best fit following the greedy algorithm 2 given the last class of
cars improved is id_last_car and updates the fitness matrix and cleft. */
int best_fit(VI &cleft, MD &fitness, int id_last_car)
{

    int id_car = best(fitness[id_last_car], cleft);
    update_fitness(id_car, fitness, cleft);
    return id_car;
}

/* Given the data of the problem, the solution using the greedy algorithm 2 is returned
and the penalization of the solution is stored in pen */
VI gen_sol_greedy2(const Data &data, VI &cleft, int &pen)
{
    const auto &[C, M, K, w, improvements] = data;

    VS stations = inicialize_stations(C, M);
    MD fitness = generate_fitness_matrix(data, cleft);

    // Calculate solution
    VI solution = base_case(data.C, cleft, fitness);
    pen += update_station_cs(improvements[solution[0]], stations, w, false);

    for (int k = 1; k < C; ++k)
    {
        solution[k] = best_fit(cleft, fitness, solution[k - 1]);
        pen += update_station_cs(improvements[solution[k]], stations, w, k == C - 1);
    }
    return solution;
}
//////////////////////////////////////////////////

//////////////////// Greedy 3 ////////////////////
/*
Idea: improve the cars of the class that supposes less penalization in that moment
Algorithm: while there are cars left of a particular class, select the available car
        that add less penalization. If there is a tie, choose the one with more
        cars left. If there is a tie, pick the one with the lowest index.
*/

/* Returns the penalization of improving next the car with improvements imp
given the current state of the stations */
int try_station(const VI &imp, VS stations, const Window &w, bool end)
{
    const auto &[ce, ne] = w;

    int M = ce.size();
    int total_penalization = 0;
    for (int i = 0; i < M; ++i)
        total_penalization += add_car(imp[i], stations[i], ce[i], ne[i], end);
    return total_penalization;
}

/* Returns the best class of cars to be improved next following the greedy algoritm 3.
It also updates the stations and cleft */
int less_pen(VI &cleft, const Data &data, const VS &stations, bool end)
{
    const auto &[C, M, K, w, improvements] = data;

    int id = 0;
    int pen = INT_MAX;
    for (int i = 0; i < K; ++i)
    {
        if (cleft[i] > 0)
        {
            int c = try_station(improvements[i], stations, w, end);
            if (c < pen or (c == pen and cleft[i] > cleft[id]))
            {
                pen = c;
                id = i;
            }
        }
    }
    --cleft[id];
    return id;
}

/* Given the data of the problem, the solution using the greedy algorithm 1 is returned
and the penalization of the solution is stored in pen */
VI gen_sol_greedy3(const Data &data, VI &cleft, int &pen)
{
    const auto &[C, M, K, w, improvements] = data;

    VI solution(C);
    VS stations = inicialize_stations(C, M);

    // The first car to be improved is the one of the class with mores cars left
    solution[0] = most_cleft_index(cleft);
    --cleft[solution[0]];
    pen += update_station_cs(improvements[solution[0]], stations, w, false);

    for (int k = 1; k < C; ++k)
    {
        solution[k] = less_pen(cleft, data, stations, k == C - 1);
        pen += update_station_cs(improvements[solution[k]], stations, w, k == C - 1);
    }

    return solution;
}
//////////////////////////////////////////////////

/* Given a data and the output file name writes in the output file the best
solution of running two greedy algorithms */
void greedy(const Data &data, VI &cleft1, const string &out)
{
    double start = now();
    VI cleft2 = cleft1; // we make a copy of cleft to use it in the second greedy
    VI cleft3 = cleft1; // we make a copy of cleft to use it in the second greedy

    int pen1 = 0;
    VI sol1 = gen_sol_greedy1(data, cleft1, pen1);
    write_solution(sol1, pen1, now() - start, out);

    int pen2 = 0;
    VI sol2 = gen_sol_greedy2(data, cleft2, pen2);
    if (pen2 < pen1)
        write_solution(sol2, pen2, now() - start, out);

    int pen3 = 0;
    VI sol3 = gen_sol_greedy3(data, cleft3, pen3);
    if (pen3 < min(pen1, pen2))
        write_solution(sol3, pen3, now() - start, out);

    // cout << pen1 << ' ' << pen2 << ' ' << pen3 << endl;
    // cout << min(min(pen3, pen2), pen1) << ' ' << (min(min(pen3, pen2), pen1) == pen1) << ' ' << (min(min(pen3, pen2), pen1) == pen2) << ' ' << (min(min(pen3, pen2), pen1) == pen3) << endl;
}

void read_input(is &in, Data &data, VI &cleft)
{
    auto &[C, M, K, w, improvements] = data;
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

    for (int i = 0; i < K; ++i)
    {
        int aux;
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
    VI cleft; // Stores the number of cars left of the i-th class that still need to be produced
    read_input(in, data, cleft);

    greedy(data, cleft, argv[2]);
}
