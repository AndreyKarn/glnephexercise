TEMPLATE = app
CONFIG -= qt
CONFIG -= app_bundle
CONFIG += console

INCLUDEPATH += "../include/libglnsvpos"
INCLUDEPATH += "../src"

isEmpty(BOOST_INCLUDE_DIR): BOOST_INCLUDE_DIR=$$(BOOST_INCLUDE_DIR)
# set by Qt Creator wizard
isEmpty(BOOST_INCLUDE_DIR): BOOST_INCLUDE_DIR="D:/Repository/glnephexercise/libglnsvpos/Boost_test/boost"
!isEmpty(BOOST_INCLUDE_DIR): INCLUDEPATH *= $${BOOST_INCLUDE_DIR}

isEmpty(BOOST_INCLUDE_DIR): {
    message("BOOST_INCLUDE_DIR is not set, assuming Boost can be found automatically in your system")
}

SOURCES += \
    ../src/func.cpp \
    ../src/glnsvpos.cpp \
    ../src/rungekutta.cpp \
    main.cpp

HEADERS += \
    ../include/libglnsvpos/func.h \
    ../include/libglnsvpos/glnsvpos.h \
    ../include/libglnsvpos/rungekutta.h \
    ../include/libglnsvpos/structures.h
