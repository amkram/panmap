#include <iostream>
#include <boost/test/unit_test.hpp>
#include "../pmi.hpp"
#include "../seeding.hpp"
#include "../tree.hpp"
#include "../genotype.hpp"
#include "PangenomeMAT.hpp"

static int getIndexFromNucleotide(char nuc) {

  switch (nuc) {
  case 'A':
  case 'a':
    return 0;
  case 'C':
  case 'c':
    return 1;
  case 'G':
  case 'g':
    return 2;
  case 'T':
  case 't':
    return 3;
  case '*':
    return 4;
  default:
    return 5;
  }
  return 5;
}

BOOST_AUTO_TEST_CASE(getNucSeq)
{
  std::cout << "Hello, World!" << std::endl;
  std::string pmat = "/home/alan/panmap/dev/examples/pmats/sars2k.pmat";

  std::cout << "Starting tests with " << pmat << std::endl;

  std::ifstream ifs(pmat);
  boost::iostreams::filtering_streambuf<boost::iostreams::input> b;
  b.push(boost::iostreams::gzip_decompressor());
  b.push(ifs);
  std::istream is(&b);

  std::cout << "Parsing Tree" << std::endl;

  PangenomeMAT::Tree *T = new PangenomeMAT::Tree(is);

  // std::cout << "brute force mutation matrix:" << std::endl;
  // seed_annotated_tree::mutationMatrices mutMat;
  // seed_annotated_tree::fillMutationMatricesFromTree(mutMat, T, 0, 0);

  std::cout << "faster mutation matrix:" << std::endl;
  seed_annotated_tree::mutationMatrices mutMat2;
  seed_annotated_tree::fillMutationMatricesFromTree_test(mutMat2, T, 100, 0.80);


}