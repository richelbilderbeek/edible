SOURCES += \
    qtmaindialog.cpp \
    qtmain.cpp

HEADERS  += qtmaindialog.h

FORMS    += qtmaindialog.ui

# Meschach
QMAKE_CFLAGS += -Wno-sign-compare
QMAKE_CFLAGS += -Wno-unused-parameter
QMAKE_CFLAGS += -Wno-unused-variable
QMAKE_CFLAGS += -Wno-unused-function
QMAKE_CFLAGS += -Wno-parentheses
QMAKE_CFLAGS += -Wno-unused-but-set-variable
QMAKE_CFLAGS += -Wno-sequence-point
include(../Meschach/Meschach.pri)

QMAKE_CFLAGS += -Werror

INCLUDEPATH += .. $$PWD $$PWD/src

include(edible.pri)

# Qt
QT += core gui widgets
