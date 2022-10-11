#include <iostream>
#include <map>
#include <vector>
#include <functional>
#include <cmath>

using func = std::function<double(std::vector<std::string>)>;

int main(int argc, char **argv) {
    using namespace std;
    map<string, func> mapa;
    mapa["sin"] = [](vector<string>x){return sin(stod(x[2]));};
    mapa["add"] = [](vector<string>x){return stod(x[2])+stod(x[3]);};
    mapa["mod"] = [](vector<string>x){return stoi(x[2])%stoi(x[3]);};

    try {
        vector<string> arguments(argv, argv + argc);
        auto SelectedFunc = mapa.at(arguments.at(1));
        cout << SelectedFunc(arguments) << endl;
    } catch (std::out_of_range aor) {
        cout << "Błąd podczas podawania argumentów, poprawna składnia to: 'sin/add/mod [x] [y]'";

    }

    return 0;
}
