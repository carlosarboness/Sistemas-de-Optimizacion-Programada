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

////////////// Functions for both greedies
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

int count_penalization_cs(const VI &sol, const Data &data)
{
    const auto &[C, M, K, w, improvements] = data;

    VS stations = inicialize_stations(C, M);

    int total_penalization = 0;
    for (int i = 0; i < C; ++i)
        total_penalization += update_station_cs(improvements[sol[i]], stations, w, i == C - 1);
    return total_penalization;
}

double calculate_cost(const Window &w, const VI &imp_i, const VI &imp_j)
{
    const auto &[ce, ne] = w;

    double cost = 0;
    for (int k = 0; k < ce.size(); ++k)
        if (imp_i[k] == 1 and imp_j[k] == 1)
            cost += ne[k] / ce[k];
    return cost;
}
////////////////////

////////////////// Greedy1
// returns the argmax of VI
int most_cleft_index(const VI &cleft)
{
    auto it = max_element(cleft.begin(), cleft.end());
    return distance(cleft.begin(), it);
}

int best_car(VI &cleft, const VD &costs)
{
    fit bf;

    for (int i = 0; i < cleft.size(); ++i)
        if (cleft[i] > 0 and (cleft[i] > bf.cl or (cleft[i] == bf.cl and costs[i] < bf.cost)))
            bf = fit{i, cleft[i], costs[i]};

    --cleft[bf.id_car];
    return bf.id_car;
}

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

VI gen_sol_greedy1(const Data &data, VI &cleft)
{

    MD costs = generate_costs_matrix(data);

    // Calculate solution
    VI solution(data.C);

    solution[0] = most_cleft_index(cleft);
    --cleft[solution[0]];
    for (int k = 1; k < data.C; ++k)
        solution[k] = best_car(cleft, costs[solution[k - 1]]);

    return solution;
}
//////////////////

/////////////// Greedy2
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

int best_fit(VI &cleft, MD &fitness, int id_last_car)
{

    int id_car = best(fitness[id_last_car], cleft);
    update_fitness(id_car, fitness, cleft);
    return id_car;
}

double calculate_cost2(const Window &w, const VI &imp_i, const VI &imp_j)
{
    const auto &[ce, ne] = w;
    double cost = 0;

    for (int k = 0; k < ce.size(); ++k)
        if (imp_i[k] == 1 and imp_j[k] == 1)
            cost += ne[k] / ce[k];

    return cost;
}

MD generate_fitness_matrix(const Data &data, const VI &cleft)
{
    const auto &[C, M, K, w, improvements] = data;

    MD fitness(K, VD(K));
    for (int i = 0; i < K; ++i)
        for (int j = 0; j < K; ++j)
        {
            double c = calculate_cost2(w, improvements[i], improvements[j]);
            fitness[i][j] = (double)cleft[j] / (c == 0 ? 1 : c); // els que queden dividit pel "cost"
        }
    return fitness;
}

VI base_case(int C, VI &cleft, MD &fitness)
{
    VI solution(C);
    int base = argmax(cleft, {});
    solution[0] = base;
    update_fitness(base, fitness, cleft);
    return solution;
}

VI gen_sol_greedy2(const Data &data, VI &cleft)
{
    MD fitness = generate_fitness_matrix(data, cleft);

    // Calculate solution
    VI solution = base_case(data.C, cleft, fitness);

    for (int k = 1; k < data.C; ++k)
        solution[k] = best_fit(cleft, fitness, solution[k - 1]);

    return solution;
}
////////////////

void greedy(const Data &data, VI &cleft, const string &out)
{
    double start = now();
    VI cleft2 = cleft;
    VI sol1 = gen_sol_greedy1(data, cleft);
    int pen1 = count_penalization_cs(sol1, data);
    write_solution(sol1, pen1, now() - start, out);
    VI sol2 = gen_sol_greedy2(data, cleft2);
    int pen2 = count_penalization_cs(sol2, data);
    if (pen2 < pen1)
        write_solution(sol2, pen2, now() - start, out);
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
