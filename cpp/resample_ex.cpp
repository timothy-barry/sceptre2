#include <iostream>
#include <random>
#include <vector>
using namespace std;

int main() {
  // number of samples to draw (WOR)
  int k = 250;
  // total number of cells
  int n_cells_total = 10000;
  // allocate the binary map
  vector<int> binary_map(n_cells_total, 0);
  // allocate the vector of random indexes
  vector<int> random_idxs(k, 0);
  // current number of successfully sampled random indexes
  int curr_count = 0;
  int random_idx = 0;

  // set up random number generation
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distr(0, n_cells_total - 1);

  // perform sampling WOR
  while(curr_count < k) {
    random_idx = distr(gen);
    if (!binary_map[random_idx]) {
      binary_map[random_idx] = 1;
      random_idxs[curr_count] = random_idx;
      curr_count ++;
    }
  }
  // finally, reset the binary map
  for (int i = 0; i < k; i ++) {
    binary_map[random_idxs[i]] = 0;
  }
  
}
