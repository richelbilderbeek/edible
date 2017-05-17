# Original code uses implicit return values
#QMAKE_CFLAGS += -Wno-implicit-int

# Original code uses implicit return values
QMAKE_CFLAGS += -Wno-sign-compare

# Original code uses unused variables
QMAKE_CFLAGS += -Wno-unused-parameter
#QMAKE_CFLAGS += -Wno-unused-variable
#QMAKE_CFLAGS += -Wno-unused-but-set-variable

# Original code uses implicit return values
#QMAKE_CXXFLAGS += -Wno-implicit-int

# Original code uses unused variables
#QMAKE_CXXFLAGS += -Wno-unused-variable
#QMAKE_CXXFLAGS += -Wno-unused-but-set-variable

#QMAKE_CFLAGS += -std=c11
#QMAKE_CFLAGS += -Werror

INCLUDEPATH += ..
INCLUDEPATH += . $$PWD/src

include(../Meschach/Meschach.pri)
include(edible.pri)

SOURCES += $$PWD/src/edible.c
