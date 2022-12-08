#include <iostream>
#include <vector>
#include <climits>
#include <chrono>
#include <cassert>
#include <fstream>
#include <algorithm>
using namespace std;

struct Class;
struct Station;

using VI = vector<int>;
using VD = vector<double>;
using VC = vector<Class>;
using VS = vector<Station>;
using MD = vector<VD>;
using is = ifstream;
using os = ofstream;

int UNDEF = -1;

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

int count_penalization_cs(const VI &sol, int C, const VI &ce, const VI &ne, const VC &classes)
{
    int M = ne.size();
    VS stations = inicialize_stations(C, M);

    int total_penalization = 0;
    for (int i = 0; i < C; ++i)
        total_penalization += update_station_cs(classes[sol[i]].improvements, stations, ce, ne, i == C - 1);
    return total_penalization;
}

int calculate_ones(int j, const VI &cleft, const VC &classes)
{
    int ones = 0;
    for (int i = 0; i < cleft.size(); ++i)
        ones += cleft[i] * classes[i].improvements[j];
    return ones;
}

int optimal(int i, int j, int C, const VI &cleft, const VI &line, int ce, int ne, const VC &classes)
{
    int ones = calculate_ones(j, cleft, classes);

    int req = 0;

    for (int k = i - 1; k >= 0 and k > i - 1 - ne; --k)
        if (line[k])
            ++req;

    return ones > 0 and req + 1 <= ce;
}

// Prec: vectors v and u must be the same size
int alike(const VI &v, const VI &u)
{
    int n = v.size();
    int total = 0;
    for (int i = 0; i < n; ++i)
        total += (v[i] == u[i]);
    return total;
}

// We also update cleft
int most_suitable(int last, const VI &opt, VI &cleft, const VC &classes)
{
    int M = classes.size();
    int class_id, suit = UNDEF;
    for (int i = 0; i < M; ++i)
    {
        int a = alike(opt, classes[i].improvements);
        if (cleft[i] > 0 and class_id != last and (a > suit or (a == suit and cleft[i] > cleft[class_id]))) // si alike és igual, agafem el que li quedin més cotxes
        {
            suit = a;
            class_id = i;
        }
    }
    --cleft[class_id];
    return class_id;
}

void update_stations(const VI &imp, VS &st)
{
    int M = st.size();
    for (int i = 0; i < M; ++i)
        st[i].line.push_back(imp[i]);
}

int best_fit(int i, int C, int last, VI &cleft, VS &st, const VI &ce, const VI &ne, const VC &classes)
{
    int M = ce.size();
    VI opt(M);
    for (int j = 0; j < M; ++j)
        opt[j] = optimal(i, j, C, cleft, st[j].line, ce[j], ne[j], classes);
    int suit = most_suitable(last, opt, cleft, classes);
    update_stations(classes[suit].improvements, st);
    return suit;
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

// returns the argmax of VI
int most_cleft_index(const VI &cleft)
{
    auto it = max_element(cleft.begin(), cleft.end());
    cout << "anem be?" << endl;
    return distance(cleft.begin(), it);
}

VI gen_sol(int C, const VI &ce, const VI &ne, const VC &classes)
{
    VI cleft = count_cleft(classes);
    VI sol(C);
    VS st = inicialize_stations(C, ce.size());
    cout << "ei!" << endl;
    sol[0] = most_cleft_index(cleft);
    --cleft[sol[0]];
    cout << "mmmm" << endl;
    cout << C << endl;
    for (int i = 0; i < C; ++i)
    {
        cout << "hola";
        if (i != 0)
            sol[i] = best_fit(i, C, sol[i - 1], cleft, st, ce, ne, classes);
    }
    return sol;
}

void greedy(const string &s, int C, const VI &ce, const VI &ne, const VC &classes)
{
    double start = now();
    VI solution = gen_sol(C, ce, ne, classes);
    double end = now();
    int penalization = count_penalization_cs(solution, C, ce, ne, classes);
    write_solution(solution, penalization, end - start, s);
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

    greedy(argv[2], C, ce, ne, classes);
}