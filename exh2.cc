#include <iostream>
#include <vector>
#include <climits>
#include <chrono>
#include <queue>
#include <cassert>
#include <fstream>
using namespace std::chrono;
using namespace std;

using VI = vector<int>;
using VD = vector<double>;
using MD = vector<VD>;

double now()
{
    return clock() / double(CLOCKS_PER_SEC);
}

ofstream open(const string &s)
{
    ofstream f;
    f.open(s, ofstream::out | ofstream::trunc);
    f.setf(ios::fixed);
    f.precision(1);
    return f;
}

struct Class
{
    int id, n; // id and # of cars of the Class
    VI imp;    // improvements
};

struct Pen
{
    int sum;      // suma total de 1's que hi ha dins la cua
    queue<int> q; // cua que conté els últims ne elements
};

void write_solution(const VI &current_sol, int pen, const double &time, const string &s)
{
    ofstream f = open(s);
    f << pen << " " << time << endl;
    bool primer = true;
    for (auto &b : current_sol)
    {
        if (primer)
            primer = false;
        else
            f << " ";
        f << b;
    }
    f << endl;
    f.close();
}

int count_pen_millora(int improv, Pen &pena, int ce_i, bool final)
{
    pena.sum -= pena.q.front();
    pena.q.pop();
    pena.q.push(improv);
    pena.sum += improv;
    int pen = max(pena.sum - ce_i, 0);
    if (final and pen > 0)
    {
        while (not pena.q.empty())
        {
            pena.sum -= pena.q.front();
            pena.q.pop();
            pen += max(pena.sum - ce_i, 0);
        }
    }
    return pen;
}

int count_pen_tot(const VI &imp, vector<Pen> &pens, const VI ce, bool final)
{
    int M = pens.size();
    int pen_tot = 0;
    for (int j = 0; j < M; ++j)
        pen_tot += count_pen_millora(imp[j], pens[j], ce[j], final);
    return pen_tot;
}

void exh_rec(int i, int current_pen, int &min_pen, VI &current_sol, vector<Pen> &pens, const VI ce, const VI ne,
             const vector<Class> &classes, VI &mkd, const VI &vec, const double &start, const string &s)
{
    int C = current_sol.size();
    if (current_pen >= min_pen)
    {
        return;
    }
    if (i == C)
    {
        cout << "good" << endl;
        min_pen = current_pen;

        double end = now();
        double elapsed_time = end - start;
        write_solution(current_sol, current_pen, elapsed_time, s);
    }
    else
    {
        for (int j = 0; j < C; ++j)
        {
            if (not mkd[j])
            {
                cout << j << ' ' << vec[j] << endl;
                mkd[j] = true;
                current_sol[i] = vec[j];
                vector<Pen> pens_rec = pens;
                exh_rec(i + 1, current_pen + count_pen_tot(classes[vec[j]].imp, pens, ce, i + 1 == C),
                        min_pen, current_sol, pens, ce, ne, classes, mkd, vec, start, s);
                mkd[j] = false;
                pens = pens_rec;
            }
        }
    }
}

//
double calculate_cost(const VI &ce, const VI &ne, const VI &imp_i, const VI &imp_j)
{
    double cost = 0;
    int M = imp_i.size();
    for (int k = 0; k < M; ++k)
        if (imp_i[k] == imp_j[k])
            cost += ne[k] / ce[k];
    return cost;
}

// returns the argmax of VI
int argmax_VI(const VI &v)
{
    int n = v.size();
    int k = -1;
    int max_elem = 0;
    for (int i = 0; i < n; ++i)
        if (v[i] > max_elem)
        {
            max_elem = v[i];
            k = i;
        }
    return k;
}

int conditions(int last_car, VI &cars_left, const VD &costs)
{
    int K = cars_left.size();
    int cl = 0;
    double cost = -1;
    int id_car = -1;
    for (int i = 0; i < K; ++i)
    {
        if (cars_left[i] > cl or (cars_left[i] == cl and costs[i] < cost))
        {
            id_car = i;
            cl = cars_left[i];
            cost = costs[i];
        }
    }
    --cars_left[id_car];
    return id_car;
}

VI gen_sol(int C, const VI &ce, const VI &ne, const vector<Class> &classes)
{
    // Generate data structues
    int K = classes.size();
    VI solution(C);
    VI cars_left(K);
    for (int i = 0; i < K; ++i)
        cars_left[i] = classes[i].n;
    MD costs(K, VD(K));
    for (int i = 0; i < K; ++i)
        for (int j = i; j < K; ++j)
            costs[i][j] = costs[j][i] = calculate_cost(ce, ne, classes[i].imp, classes[j].imp);

    // Calculate solution
    solution[0] = argmax_VI(cars_left);
    --cars_left[solution[0]];
    for (int k = 1; k < C; ++k)
        solution[k] = conditions(solution[k - 1], cars_left, costs[solution[k - 1]]);
    return solution;
}
//

void exh(const string &s, int C, const VI &ce, const VI &ne, const vector<Class> &classes)
{
    int M = ne.size();
    VI current_sol(C);
    vector<Pen> pens(M);
    for (int i = 0; i < M; ++i)
    {
        pens[i].sum = 0;
        for (int j = 0; j < ne[i]; ++j)
            pens[i].q.push(0);
    }
    VI vec = gen_sol(C, ce, ne, classes);
    VI mkd(C, false);
    double start = now();
    int min_pen = INT_MAX;
    exh_rec(0, 0, min_pen, current_sol, pens, ce, ne, classes, mkd, vec, start, s);
}

void read_input(ifstream &f, int &C, int &M, int &K, VI &ce, VI &ne, vector<Class> &classes)
{
    f >> C >> M >> K;
    ce.resize(M);
    ne.resize(M);
    classes.resize(K);
    for (auto &capacity : ce)
        f >> capacity;
    for (auto &num_cars : ne)
        f >> num_cars;
    for (int i = 0; i < K; ++i)
    {
        f >> classes[i].id >> classes[i].n;
        VI imp(M);
        for (auto &imprv : imp)
            f >> imprv;
        classes[i].imp = imp;
    }
}

int main(int argc, char *argv[])
{

    assert(argc == 3);

    ifstream f(argv[1]);

    int C, M, K;
    VI ce, ne;
    vector<Class> classes;

    read_input(f, C, M, K, ce, ne, classes);

    exh(argv[2], C, ce, ne, classes);
}