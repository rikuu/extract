// Copyright 2016 Riku Walve

#include <extract.hpp>

int main(int argc, char* argv[]) {
  try {
    Extract().run(argc, argv);
  } catch (Exception& e) {
    std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
