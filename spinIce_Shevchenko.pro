#mpi
LIBS += "C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x86\msmpi.lib"
LIBS += "C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x86\msmpifec.lib"
LIBS += "C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x86\msmpifmc.lib"
INCLUDEPATH += "C:\Program Files (x86)\Microsoft SDKs\MPI\Include"
DEPENDPATH += "C:\Program Files (x86)\Microsoft SDKs\MPI\Include"
#

QT += core
QT -= gui

TARGET = spinIce_Shevchenko
CONFIG += console c++11
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp


LIBS += -L$$PWD/../partsEngine -lPartsEngine
INCLUDEPATH += ../partsEngine/
