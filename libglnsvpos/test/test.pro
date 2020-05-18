include(gtest_dependency.pri)

INCLUDEPATH += "googletest/include"
INCLUDEPATH += "googletest/src"
INCLUDEPATH += "../include/libglnsvpos"

TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG += thread
CONFIG -= qt

HEADERS += \
        check_suites.h \
        tst_test2.h

SOURCES += \
        check_main.cpp \
        check_position.cpp \
        main.cpp
