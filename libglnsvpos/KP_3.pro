TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        src/diffs.cpp \
        src/glnsvpos.cpp \
        src/main.cpp \
        src/rungekutta.cpp

HEADERS += \
    include/libglnsvpos/diffs.h \
    include/libglnsvpos/glnsvpos.h \
    include/libglnsvpos/rungekutta.h

#win32:CONFIG(release, debug|release): LIBS += -LD:/Programs/QT/Tools/dwarfstack/lib/ -ldwarfstack
#else:win32:CONFIG(debug, debug|release): LIBS += -LD:/Programs/QT/Tools/dwarfstack/lib/ -ldwarfstackd
#else:unix: LIBS += -LD:/Programs/QT/Tools/dwarfstack/lib/ -ldwarfstack

#INCLUDEPATH += D:/Programs/QT/Tools/dwarfstack/include
#DEPENDPATH += D:/Programs/QT/Tools/dwarfstack/include
