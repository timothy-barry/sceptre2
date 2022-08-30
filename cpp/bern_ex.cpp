#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
using namespace std;

int main()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution d(0.01);
    int n = 100000000;

    for(int i=0; i < n; i++) {
        d(gen);
    }
}
