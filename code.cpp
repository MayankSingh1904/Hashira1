#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <nlohmann/json.hpp>
#include <gmpxx.h>
using json = nlohmann::json;
using namespace std;
struct Point {
    mpz_class x;
    mpz_class y;
};
mpz_class convertToDecimal(const string& value, int base) {
    mpz_class result;
    mpz_set_str(result.get_mpz_t(), value.c_str(), base);
    return result;
}
mpz_class lagrangeInterpolationAtZero(const vector<Point>& points) {
    mpz_class result = 0;
    size_t k = points.size();
    for (size_t i = 0; i < k; ++i) {
        mpz_class term = points[i].y;
        mpz_class numerator = 1;
        mpz_class denominator = 1;
        for (size_t j = 0; j < k; ++j) {
            if (i == j) continue;
            numerator *= -points[j].x;
            denominator *= (points[i].x - points[j].x);
        }
        mpz_class inv;
        if (mpz_invert(inv.get_mpz_t(), denominator.get_mpz_t(), 0) == 0) {
            term *= numerator;
            term /= denominator;
        } else {
            term *= numerator;
            term *= inv;
        }
        result += term;
    }
    return result;
}
mpz_class solveFromJsonFile(const string& filename) {
    ifstream inFile(filename);
    if (!inFile) {
        cerr << "Unable to open file: " << filename << endl;
        exit(1);
    }
    json j;
    inFile >> j;
    int n = j["keys"]["n"];
    int k = j["keys"]["k"];
    vector<Point> allPoints;
    for (auto& [key, value] : j.items()) {
        if (key == "keys") continue;
        int x = stoi(key);
        int base = stoi(value["base"].get<string>());
        string encoded = value["value"];
        mpz_class y = convertToDecimal(encoded, base);
        allPoints.push_back({x, y});
    }
    vector<Point> subset(allPoints.begin(), allPoints.begin() + k);
    return lagrangeInterpolationAtZero(subset);
}
int main() {
    mpz_class result1 = solveFromJsonFile("testcase1.json");
    mpz_class result2 = solveFromJsonFile("testcase2.json");
    cout << "Secret from Testcase 1: " << result1.get_str() << endl;
    cout << "Secret from Testcase 2: " << result2.get_str() << endl;
    return 0;
}
