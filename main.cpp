#include <iostream>
#include <utility>
#include <vector>
#include <random>

using namespace std;

random_device rd;
mt19937 mt_generator(rd());

vector<vector<int>> genotypeGenerator(int size) {
    vector<vector<int>> genotype;
    vector<int> chromosomes;

    for (int i = 0; i < size; ++i) {
        chromosomes.clear();

        for (int j = 0; j < size ; ++j) {
            uniform_int_distribution<int> uni(0, 1);
            chromosomes.push_back(uni(mt_generator));
        }

        genotype.push_back(chromosomes);
    }

    return genotype;
}

pair<double, double> decodePhenotype(vector<int> chromosomes) {
    double x, y = 0.0;

    for (int i = 0; i < chromosomes.size() / 2; i++) {
        x = x * 2 + chromosomes[i];
    }

    for (int i = (int) chromosomes.size() / 2; i < chromosomes.size(); i++) {
        y = y * 2 + chromosomes[i];
    }

    x = x / pow(2.0, (chromosomes.size() / 2 - 4)) - 8;
    y = y / pow(2.0, (chromosomes.size() / 2 - 4)) - 8;

    return {x, y};
}

auto beale = [](pair<double, double> pair) {
    double x = pair.first;
    double y = pair.second;

    return (pow((1.5 - x + (x * y)),2)) + (pow(2.25 - x + (x * pow(y,2)),2))+(pow(2.625 - x + x * pow(y,3),2));
};


double fitness(vector<int> chromosomes) {
    return 1.0 / (1.0 + abs(beale(decodePhenotype(move(chromosomes)))));
}

int main() {
    for (auto &element: genotypeGenerator(100 + (23132 % 10) * 2)) {
        auto decoded = decodePhenotype(element);
        cout << "X: " << decoded.first << "\nY: " << decoded.second << "\nMaksymalizacja: " << fitness(element) <<"\n"<< endl;
    }

    return 0;
}