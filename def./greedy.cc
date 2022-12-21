#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
using namespace std;

/* GREEDY ALGORITHM
Idea: improve the cars of the class that supposes less penalization in that moment
Algorithm: while there are cars left of a particular class, select the available car
        that addS less penalization. If there is a tie, choose the one with more
        cars left. If there is a tie, pick the one with the lowest index.
*/

struct Station;
using VI = vector<int>;
using MI = vector<VI>;
using VS = vector<Station>;

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

// C
double now()
{
    return clock() / double(CLOCKS_PER_SEC);
}
// C
ofstream open(const string &out)
{
    ofstream f;
    f.open(out, ofstream::out | ofstream::trunc);
    f.setf(ios::fixed);
    f.precision(1);
    return f;
}

// C /*  */
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

/* Updates each individual line with the bit indicating if the improvement
of the station is required and retuns the penalitzations of it */
int update_station(int bit, Station &st, int ce_i, int ne_i, bool end)
{
    // updates the window
    int uw = st.line.size() - ne_i; // uw is the position in the st.line of improvement that is no longer in the station
    if (uw >= 0)
        st.requirements -= st.line[uw];

    st.line.push_back(bit);
    st.requirements += bit;

    int penalitzation = max(st.requirements - ce_i, 0);

    // when is the last car to be added in the production line,
    // we add to the penalization of incomplete final windows
    if (end)
    {
        ++uw;
        while (st.requirements > 0)
        {
            st.requirements -= st.line[uw];
            penalitzation += max(st.requirements - ce_i, 0);
            ++uw;
        }
    }

    return penalitzation;
}

/* Update Production Lines: updates the production lines of all the stations because of the
addition of a new car in the production line.
The function returns the total sum of penalizations of each line, updated with the new car */
int UPL(const VI &imp, VS &stations, const Window &w, bool end)
{
    const auto &[ce, ne] = w; // Get the data

    int M = ce.size();
    int total_penalization = 0;
    for (int i = 0; i < M; ++i)
        total_penalization += update_station(imp[i], stations[i], ce[i], ne[i], end);
    return total_penalization;
}

/* Returns a vector with the inicialized stations of the algorithm.
Firstly the requiremnts are equal to zero and the line is empty */
VS inicialize_stations(int C, int M)
{
    VS stations(M); // empty vector of station wich will be filled during the search

    int requirements = 0;
    VI inicial_line;
    inicial_line.reserve(C); // we reserve size C for the vector (maximum at the end of the search)

    for (Station &st : stations)
        st = Station{requirements, inicial_line};

    return stations;
}

/* Returns the argmax of the VI cleft */
int most_cleft_index(const VI &cleft)
{
    auto it = max_element(cleft.begin(), cleft.end());
    return distance(cleft.begin(), it);
}

/* Returns the penalization of improving next the car with improvements imp
given the current state of the stations */
int try_station(const VI &imp, VS stations, const Window &w, bool end)
{
    const auto &[ce, ne] = w;

    int M = ce.size();
    int total_penalization = 0;
    for (int i = 0; i < M; ++i)
        total_penalization += update_station(imp[i], stations[i], ce[i], ne[i], end);
    return total_penalization;
}

/* Returns the best class of cars to be improved next following the greedy algoritm.
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
VI gen_sol_greedy(const Data &data, VI &cleft, int &pen)
{
    const auto &[C, M, K, w, improvements] = data; //

    VI solution(C);
    VS stations = inicialize_stations(C, M);

    // The first car to be improved is the one of the class with mores cars left
    solution[0] = most_cleft_index(cleft);
    --cleft[solution[0]];
    pen += UPL(improvements[solution[0]], stations, w, false);

    for (int k = 1; k < C; ++k)
    {
        solution[k] = less_pen(cleft, data, stations, k == C - 1);
        pen += UPL(improvements[solution[k]], stations, w, k == C - 1);
    }

    return solution;
}

/* Given a data and the output file name writes in the output file the best
solution of running two greedy algorithms */
void greedy(const Data &data, VI &cleft, const string &out)
{
    double start = now();
    int pen = 0;
    VI sol = gen_sol_greedy(data, cleft, pen);
    write_solution(sol, pen, now() - start, out);
}

// C
void read_input(ifstream &in, Data &data, VI &cleft)
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
    // An error sentence will be printed if there is not the required information in the command line
    if (argc != 3)
    {
        cout << "Syntax: " << argv[0] << " input_file output_file" << endl;
        exit(1);
    }

    ifstream in(argv[1]);
    Data data;
    VI cleft; // Stores the number of cars left of the i-th class that still need to be produced

    read_input(in, data, cleft); // The imput infromation of the imput file in will be readed and stored int data and cleft

    greedy(data, cleft, argv[2]); // Solution given by the greedy algorithm is written in the output file argv[2]
}
