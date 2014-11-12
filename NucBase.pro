#-------------------------------------------------
#
# Project created by QtCreator 2011-05-16T11:16:58
#
#-------------------------------------------------

QT       += core gui widgets

QMAKE_CXXFLAGS +=  -s -Wall -ansi -pedantic -std=c++0x -Werror -fopenmp

LIBS += -ldivsufsort -fopenmp

TARGET = NucBase
TEMPLATE = app
RC_FILE = dna.rc

SOURCES += main.cpp\
        mainwindow.cpp \
    nucbase.cpp \
    computethread.cpp \
    nucsequences.cpp \
    convertdialog.cpp

HEADERS  += \
    nucbase.hxx \
    nucbase.hpp \
    computethread.hpp \
    mainwindow.hpp \
    nucsequences.hpp \
    nucsequences.hxx \
    convertdialog.hpp

FORMS    += mainwindow.ui \
    convertdialog.ui

ICON     = dna.icns






