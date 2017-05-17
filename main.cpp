#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <iomanip>

#include "helper.h"

double read_determinant(const std::string& filename)
{
  const auto text = file_to_vector(filename);
  const auto line = text[ text.size() - 2];
  const std::string to_remove{"Determinant of expectation matrix = "};
  const auto number = line.substr(to_remove.size(), line.size() - to_remove.size());
  return std::stod(number);
}

int main()
{
  if (!is_regular_file("edible"))
  {
    std::cout << "Executable 'edible' not found\n";
    return 1;
  }

  //const std::string tree_filename{"~/GitHubs/edible/test/20170517_alignment_information_3.tree"};
  const std::string tree_filename{"~/GitHubs/edible/test/primate.tree"};
  const std::string out_filename{"out.txt"};


  {
    std::string cmd = "./edible -o " + tree_filename + " " + out_filename + " >/dev/null";
    const int error = std::system(cmd.c_str());
    if (error)
    {
      std::cout << "Command '" << cmd << "' failed\n";
      return error;
    }
  }
  std::cout << std::setprecision(99) << read_determinant(out_filename) << '\n';
}
