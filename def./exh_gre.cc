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
using VD = vector<double>;
using MD = vector<VD>;
using is = ifstream;
using os = ofstream;

int UNDEF = -1;

// Stores all the information given by the input
struct Data
{
    int C, M, K; // Number of cars, improvments and classes respectively
    VI ce, ne;   // Stores the ce (ne) of a particular station in the position of the vector
    VC classes;  // Stores each class in the position of the vector
};

struct Class
{
    int ncars;       // Numbers of cars of a particular class that needs to be produces
    VI improvements; // Vector of integers (1/0) that indicate if the class requires i-th improvment
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

/* Restores the stations after a recursive call */
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
addition of a new car in the production line.
The function returns the total sum of penalizations of each line, updated with the new car */
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

int _cs(const VI &line, int ce, int ne)
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

bool add_one(const VI &line, int ones, int zeros, int &req, int ce)
{
    if (ones <= 0)
        return false;

    if (zeros <= 0)
        return true;

    if (req + 1 <= ce)
    {
        req += 1;
        return true;
    }
    return false;
}

void add_bit(VI &line, int bit, int &zeros, int &ones)
{
    line.push_back(bit);
    (bit ? --ones : --zeros);
}

int lb_station(int j, int i, const VI &cs, const VI &cleft, VI line, int ce, int ne, const VC &classes)
{
    int C = cs.size();
    int ones = calculate_ones(j, cleft, classes);
    int zeros = (C - i - 1) - ones;

    int req = 0;

    for (int k = i; k >= 0 and k > i - ne; --k)
        if (line[k])
            ++req;

    for (; zeros > 0 or ones > 0; ++i)
    {
        if (int lw = i - ne + 1; lw >= 0)
            req -= line[lw];

        (add_one(line, ones, zeros, req, ce) ? add_bit(line, 1, zeros, ones) : add_bit(line, 0, zeros, ones));
    }

    return _cs(line, ce, ne);
}

/* Returns a lower bound of the current solution cs.
This lower bound is the sum of the minimum penalization of each station if the improvments were in the best possible way */
int lower_bound(int i, int LowerBound, const VI &cs, const VI &cleft, const VS &st, const VI &ce, const VI &ne, const VC &cl)
{
    if (i == -1)
        return 0;
    for (int j = 0; j < ne.size(); ++j)
        LowerBound += lb_station(j, i, cs, cleft, st[j].line, ce[j], ne[j], cl);

    return LowerBound;
}

VI update(VI generator, int K)
{
    // random_shuffle(generator.begin(), generator.end());

    for (int i = 0; i < K; ++i)
        generator[i] = (generator[i] + 1) % K;

    return generator;
}

/*  */
void exhaustive_search_rec(int i, int cp, int &mp, VI &cs, VI &cleft, VS &stations, const Data &data,
                           const double &start, const string &out, const VI &generator)
{
    const auto &[C, M, K, ce, ne, classes] = data;

    if (now() - start < 60)
    {
        if (i == C)
        {
            mp = cp;
            write_solution(cs, cp, now() - start, out);
        }
        if (lower_bound(i - 1, cp, cs, cleft, stations, ce, ne, classes) < mp)
        {
            for (int cl : generator)
            {
                if (cleft[cl] > 0)
                {
                    --cleft[cl];
                    cs[i] = cl;

                    if (int up = UPL(classes[cl].improvements, stations, ce, ne, i + 1 == C); up + cp < mp)
                        exhaustive_search_rec(i + 1, cp + up, mp, cs, cleft, stations, data, start, out, update(generator, K));

                    restore(stations, ne);
                    ++cleft[cl];
                }
            }
        }
    }
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

/* Returns a vector containing for each index the number of cars
that are left to be produced for that class */
VI count_cleft(const VC &classes)
{
    int K = classes.size();
    VI cleft(K);
    for (int i = 0; i < K; ++i)
        cleft[i] = classes[i].ncars;
    return cleft;
}

//////// Required functions of the greedy algorithm
// returns the argmax of VI
int most_cleft_index(const VI &cleft)
{
    auto it = max_element(cleft.begin(), cleft.end());
    return distance(cleft.begin(), it);
}

int best_fit(VI &cleft, const VD &costs)
{
    fit bf;

    for (int i = 0; i < cleft.size(); ++i)
        if (cleft[i] > 0 and (cleft[i] > bf.cl or (cleft[i] == bf.cl and costs[i] < bf.cost)))
            bf = fit{i, cleft[i], costs[i]};

    --cleft[bf.id_car];
    return bf.id_car;
}

double calculate_cost(const VI &ce, const VI &ne, const VI &imp_i, const VI &imp_j)
{
    double cost = 0;

    for (int k = 0; k < ce.size(); ++k)
        if (imp_i[k] == 1 and imp_j[k] == 1)
            cost += ne[k] / ce[k];

    return cost;
}

MD generate_costs_matrix(const Data &data)
{
    const auto &[C, M, K, ce, ne, classes] = data;

    MD costs(K, VD(K));
    for (int i = 0; i < K; ++i)
        for (int j = i; j < K; ++j)
        {
            double c = calculate_cost(ce, ne, classes[i].improvements, classes[j].improvements);
            costs[i][j] = costs[j][i] = c;
        }
    return costs;
}

VI gen_sol(const Data &data)
{
    const auto &[C, M, K, ce, ne, classes] = data;
    VI cleft = count_cleft(classes);
    MD costs = generate_costs_matrix(data);
    // Calculate solution
    VI sol(C);

    sol[0] = most_cleft_index(cleft);
    --cleft[sol[0]];
    for (int k = 1; k < C; ++k)
        sol[k] = best_fit(cleft, costs[sol[k - 1]]);

    return sol;
}

/* Updates each individual line and retuns the penalitzations of it */
int add_car(int bit, Station &st, int ce, int ne, bool end)
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
        for (; st.requirements > 0; ++uw)
        {
            st.requirements -= st.line[uw];
            penalitzation += max(st.requirements - ce, 0);
        }
    }
    return penalitzation;
}

int update_station_cs(const VI &imp, VS &stations, const VI &ce, const VI &ne, bool end)
{
    int M = ce.size();
    int total_penalization = 0;
    for (int i = 0; i < M; ++i)
        total_penalization += add_car(imp[i], stations[i], ce[i], ne[i], end);
    return total_penalization;
}

int count_penalization_cs(const VI &sol, const Data &data)
{
    const auto &[C, M, K, ce, ne, classes] = data;

    VS stations = inicialize_stations(C, M);

    int total_penalization = 0;
    for (int i = 0; i < C; ++i)
        total_penalization += update_station_cs(classes[sol[i]].improvements, stations, ce, ne, i == C - 1);
    return total_penalization;
}
////////

/* Returns the initial generator */
VI create_generator(int K)
{
    VI generator(K);
    for (int i = 0; i < K; ++i)
        generator[i] = i;
    return generator;
}

/* Given a data and the output file name writes in the output file the solution of the exhaustive search */
void exhaustive_search(const Data &data, const string &out)
{
    const auto &[C, M, K, ce, ne, classes] = data;

    VI cleft = count_cleft(classes); // Stores the number of cars left of the i-th class that still need to be produced
    VS stations = inicialize_stations(C, M);
    VI generator = create_generator(K); // Vector that inticates the order to do the exhausitve search

    double start = now(); // inicialize the counter

    // The first solution is the one using the greedy algorithm
    VI cs = gen_sol(data);                    // current solution
    int mp = count_penalization_cs(cs, data); // minimum penalization
    write_solution(cs, mp, now() - start, out);

    exhaustive_search_rec(0, 0, mp, cs, cleft, stations, data, start, out, generator);
}

void read_input(is &in, Data &data)
{
    auto &[C, M, K, ce, ne, classes] = data;
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
        int aux; // ns com es pot fer per no llegir algo
        in >> aux >> classes[i].ncars;
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
    Data data;

    read_input(in, data);

    exhaustive_search(data, argv[2]);
}
