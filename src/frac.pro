TEMPLATE	      = app
CONFIG		     += qt warn_on force_debug_info thread c++17
FORMS		      = mainwindow.ui gradeditor.ui
HEADERS		      = include/genkernel.h \
			include/bitarray.h \
			include/fpvec.h \
			include/fractal.h \
			include/gpuhandler.h \
                        include/gradeditor.h \
                        include/mainwindow.h \
                        include/util-widgets.h

SOURCES		      = bitarray.cc \
			colors.cc \
                        genkernel.cc \
                        gradeditor.cc \
			mainfrac.cc \
                        util-widgets.cc \
                        vparith.cc

isEmpty(PREFIX) {
PREFIX = /usr/local
}

TARGET                = frac
DATADIR               = $$PREFIX/share/frac
DOCDIR                = $$PREFIX/share/doc/frac

*-g++ {
QMAKE_CXXFLAGS += -fno-diagnostics-show-caret
}

unix:INCLUDEPATH      += include
win32:INCLUDEPATH     += include

!win32:DEFINES       += "DATADIR=\\\"$$DATADIR\\\""
!win32:DEFINES       += "DOCDIR=\\\"$$DOCDIR\\\""
release:DEFINES      += NO_CHECK
win32:DEFINES        += QT_DLL QT_THREAD_SUPPORT HAVE_CONFIG_H _USE_MATH_DEFINES

target.path           = $$PREFIX/bin
INSTALLS += target

unix {
  CUDA_DIR = $$system(which nvcc | sed 's,/bin/nvcc,,')

  INCLUDEPATH += $$CUDA_DIR/include
  QMAKE_LIBDIR += $$CUDA_DIR/lib

  LIBS += -lcuda
}

QT += widgets gui

RESOURCES += \
    frac.qrc
