# Thanks to Meschach
QMAKE_CFLAGS += -Wno-sign-compare
QMAKE_CFLAGS += -Wno-unused-parameter
QMAKE_CFLAGS += -Wno-unused-variable
QMAKE_CFLAGS += -Wno-unused-function
QMAKE_CFLAGS += -Wno-parentheses
QMAKE_CFLAGS += -Wno-unused-but-set-variable
QMAKE_CFLAGS += -Wno-sequence-point

QMAKE_CFLAGS += -Werror

INCLUDEPATH += ..
INCLUDEPATH += . $$PWD/src

include(../Meschach/Meschach.pri)
include(edible.pri)

SOURCES += $$PWD/main.c
