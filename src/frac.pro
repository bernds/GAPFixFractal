TEMPLATE	      = app
CONFIG		     += qt warn_on force_debug_info thread c++17
FORMS		      = mainwindow.ui \
			location.ui \
			rotation.ui \
			gradeditor.ui \
			batchrender.ui \
			hybriddialog.ui \
                        prefs.ui

HEADERS		      = include/genkernel.h \
			include/batchrender.h \
			include/bitarray.h \
			include/fpvec.h \
			include/fractal.h \
			include/gpuhandler.h \
                        include/gradeditor.h \
                        include/hybriddialog.h \
                        include/locationdialog.h \
                        include/mainwindow.h \
                        include/renderer.h \
                        include/rotationdialog.h \
                        include/settings.h \
                        include/util-widgets.h

SOURCES		      = batchrender.cc \
                        bitarray.cc \
			colors.cc \
                        genkernel.cc \
                        gradeditor.cc \
                        gpuhandler.cc \
                        hybriddialog.cc \
                        locationdialog.cc \
                        mainfrac.cc \
                        renderer.cc \
                        rotationdialog.cc \
                        settings.cc \
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
}
win32 {
  CUDA_DIR=$$(CUDA_PATH)
}
INCLUDEPATH += $$CUDA_DIR/include
QMAKE_LIBDIR += $$CUDA_DIR/lib
QMAKE_LIBDIR += $$CUDA_DIR/lib64
QMAKE_LIBDIR += $$CUDA_DIR/lib/x64
LIBS += -lcuda

QT += widgets gui

RESOURCES += \
    frac.qrc
