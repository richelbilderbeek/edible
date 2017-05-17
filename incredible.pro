SOURCES += $$PWD/main.cpp
include(incredible.pri)

# C++14
CONFIG += c++14
QMAKE_CXXFLAGS += -std=c++14

# High warnings levels
QMAKE_CXXFLAGS += -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -Weffc++ -Werror

# Allow debug and release mode
CONFIG += debug_and_release

# In debug mode, turn on gcov, memcache and UBSAN
CONFIG(debug, debug|release) {

  # gcov
  QMAKE_CXXFLAGS += -fprofile-arcs -ftest-coverage
  LIBS += -lgcov

  # helgrind, for helgrind and memcheck
  QMAKE_LFLAGS += -pthread -Wl,--no-as-needed

  # UBSAN
  QMAKE_CXXFLAGS += -fsanitize=undefined
  QMAKE_LFLAGS += -fsanitize=undefined
  LIBS += -lubsan
}

message(Host name: $$QMAKE_HOST.name)
contains(QMAKE_HOST.name,fwn-biol-132-102) {
  message("Host is university computer in my office")
  QMAKE_CXX = g++-5
  QMAKE_LINK = g++-5
  QMAKE_CC = gcc-5
}
