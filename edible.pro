# Original code uses void main
QMAKE_CFLAGS += -Wno-main

# Original code uses implicit return values
QMAKE_CFLAGS += -Wno-implicit-int

# Original code uses unused variables
QMAKE_CFLAGS += -Wno-unused-variable
QMAKE_CFLAGS += -Wno-unused-but-set-variable


# Original code uses void main
QMAKE_CXXFLAGS += -Wno-main

# Original code uses implicit return values
QMAKE_CXXFLAGS += -Wno-implicit-int

# Original code uses unused variables
QMAKE_CXXFLAGS += -Wno-unused-variable
QMAKE_CXXFLAGS += -Wno-unused-but-set-variable


INCLUDEPATH += ..
INCLUDEPATH += . $$PWD/src

include(../Meschach/Meschach.pri)

# Input
HEADERS += \
    $$PWD/src/edible.h \
    $$PWD/src/gtr.h \
    $$PWD/src/new_models.h \
    $$PWD/src/variables.h

SOURCES += \
    $$PWD/src/edible.c \
    $$PWD/src/gtr.c \
    $$PWD/src/llh.c \
    $$PWD/src/matrix.c \
    $$PWD/src/new_models.c \
    $$PWD/src/options.c \
    $$PWD/src/partial.c \
    $$PWD/src/random.c \
    $$PWD/src/read.c \
    $$PWD/src/tree.c \
    $$PWD/src/utility.c
