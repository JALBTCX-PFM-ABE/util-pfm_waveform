INCLUDEPATH += /c/PFM_ABEv7.0.0_Win64/include
LIBS += -L /c/PFM_ABEv7.0.0_Win64/lib -lCHARTS -lnvutility -lpfm -lgdal -lxml2 -lpoppler -lz -liconv -lstdc++ -lm
DEFINES += NVWIN3X
CONFIG += console
CONFIG -= qt
QMAKE_LFLAGS += 
######################################################################
# Automatically generated by qmake (2.01a) Wed Jan 22 13:47:14 2020
######################################################################

TEMPLATE = app
TARGET = pfm_waveform
DEPENDPATH += .
INCLUDEPATH += .

# Input
HEADERS += pfm_waveform.h version.h
SOURCES += main.c process_waveforms.c
