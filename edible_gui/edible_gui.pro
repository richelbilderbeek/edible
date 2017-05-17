QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TEMPLATE = app

SOURCES += main.cpp\
        qtmaindialog.cpp

HEADERS  += qtmaindialog.h

FORMS    += qtmaindialog.ui

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

#QMAKE_CFLAGS += -std=c11


INCLUDEPATH += ..
INCLUDEPATH += . $$PWD/src

include(../Meschach/Meschach.pri)
include(edible.pri)

SOURCES += $$PWD/src/edible.c
