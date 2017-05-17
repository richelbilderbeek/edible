#ifndef HELPER_H
#define HELPER_H

#include <string>
#include <vector>

/// FileToVector reads a file and converts it to a std::vector<std::string>
/// From http://www.richelbilderbeek.nl/CppFileToVector.htm
std::vector<std::string> file_to_vector(const std::string& filename);

////Determines if a filename is a regular file
/// From http://www.richelbilderbeek.nl/CppIsRegularFile.htm
bool is_regular_file(const std::string& filename) noexcept;


#endif // HELPER_H
