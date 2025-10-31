#include "seed_annotated_tree.hpp"
#include <fstream>
#include <stdexcept>
#include <vector>
#include <string>
#include <sstream>

// Helper function from panmanUtils
namespace {
void stringSplit(const std::string& s, char delimiter, std::vector<std::string>& result) {
    result.clear();
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        result.push_back(token);
    }
}
}

namespace seed_annotated_tree {

void fillMutationMatricesFromFile(mutationMatrices &mutMat, std::ifstream &inf) {
  std::string line;
  int idx = 0;
  while (getline(inf, line)) {
    std::vector<double> probs;
    std::vector<std::string> fields;
    stringSplit(line, ' ', fields);
    for (const auto &f : fields) {
      probs.push_back(std::stod(f));
    }
    if (probs.size() == 0) {
      break;
    }

    if (idx < 4) {
      if (probs.size() != 4) {
        throw std::invalid_argument(
            "Received invalid mutation matrix (.mm) file");
      }

      std::vector<double> probs;
      for (const auto &f : fields) {
        probs.push_back(std::stod(f));
      }
      mutMat.submat[idx] = std::move(probs);
    } else if (idx == 4) {
      if (probs.size() < 1) {
        throw std::invalid_argument(
            "Received invalid mutation matrix (.mm) file");
      }

      for (const auto &f : fields) {
        std::vector<std::string> subFields;
        stringSplit(f, ':', subFields);
        int64_t size = std::stoll(subFields[0]);
        double prob = std::stod(subFields[1]);
        mutMat.insmat[size] = prob;
      }
    } else if (idx == 5) {
      if (probs.size() < 1) {
        throw std::invalid_argument(
            "Received invalid mutation matrix (.mm) file");
      }

      for (const auto &f : fields) {
        std::vector<std::string> subFields;
        stringSplit(f, ':', subFields);
        int64_t size = std::stoll(subFields[0]);
        double prob = std::stod(subFields[1]);
        mutMat.delmat[size] = prob;
      }
    }
    idx++;
  }

  if (idx != 6) {
    throw std::invalid_argument("Received invalid mutation matrix (.mm) file");
  }
  mutMat.filled = true;
}

void scaleMutationMatrices(mutationMatrices& mutMat, double mutation_rate) {
  // Scale the mutation matrices by the given mutation rate
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      if (i != j) {
        mutMat.submat[i][j] *= mutation_rate;
      }
    }
    // Normalize row to sum to 1
    double rowSum = 0.0;
    for (int j = 0; j < 4; j++) {
      rowSum += mutMat.submat[i][j];
    }
    if (rowSum > 0) {
      for (int j = 0; j < 4; j++) {
        mutMat.submat[i][j] /= rowSum;
      }
    }
  }
}

} // namespace seed_annotated_tree
