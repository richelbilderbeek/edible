# Original code uses void main
QMAKE_CFLAGS += -Wno-main

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

INCLUDEPATH += . src

INCLUDEPATH += ../Meschach
HEADERS += ../Meschach/*.h

SOURCES += \
    ../Meschach/arnoldi.c \
    ../Meschach/bdfactor.c \
    ../Meschach/bkpfacto.c \
    ../Meschach/chfactor.c \
    ../Meschach/conjgrad.c \
    ../Meschach/copy.c \
    #../Meschach/dmacheps.c \
    ../Meschach/err.c \
    ../Meschach/extras.c \
    ../Meschach/fft.c \
    #../Meschach/fmacheps.c \
    ../Meschach/givens.c \
    ../Meschach/hessen.c \
    ../Meschach/hsehldr.c \
    ../Meschach/init.c \
    #../Meschach/iotort.c \
    ../Meschach/iter0.c \
    ../Meschach/iternsym.c \
    ../Meschach/itersym.c \
    #../Meschach/itertort.c \
    ../Meschach/ivecop.c \
    ../Meschach/lanczos.c \
    ../Meschach/lufactor.c \
    ../Meschach/machine.c \
    ../Meschach/matlab.c \
    ../Meschach/matop.c \
    ../Meschach/matrixio.c \
    #../Meschach/maxint.c \
    ../Meschach/meminfo.c \
    ../Meschach/memory.c \
    ../Meschach/memstat.c \
    #../Meschach/memtort.c \
    ../Meschach/mfunc.c \
    #../Meschach/mfuntort.c \
    ../Meschach/norm.c \
    ../Meschach/otherio.c \
    ../Meschach/pxop.c \
    ../Meschach/qrfactor.c \
    ../Meschach/schur.c \
    ../Meschach/solve.c \
    ../Meschach/sparse.c \
    ../Meschach/sparseio.c \
    ../Meschach/spbkp.c \
    ../Meschach/spchfctr.c \
    ../Meschach/splufctr.c \
    ../Meschach/sprow.c \
    ../Meschach/spswap.c \
    #../Meschach/sptort.c \
    ../Meschach/submat.c \
    ../Meschach/svd.c \
    ../Meschach/symmeig.c \
    #../Meschach/torture.c \
    #../Meschach/tutadv.c \
    #../Meschach/tutorial.c \
    ../Meschach/update.c \
    ../Meschach/vecop.c \
    ../Meschach/version.c \
    ../Meschach/zcopy.c \
    ../Meschach/zfunc.c \
    ../Meschach/zgivens.c \
    ../Meschach/zhessen.c \
    ../Meschach/zhsehldr.c \
    ../Meschach/zlufctr.c \
    ../Meschach/zmachine.c \
    ../Meschach/zmatio.c \
    ../Meschach/zmatlab.c \
    ../Meschach/zmatop.c \
    ../Meschach/zmemory.c \
    ../Meschach/znorm.c \
    ../Meschach/zqrfctr.c \
    ../Meschach/zschur.c \
    ../Meschach/zsolve.c \
    #../Meschach/ztorture.c \
    ../Meschach/zvecop.c

# Input
HEADERS += src/edible.h \
           src/gtr.h \
           src/new_models.h \
           src/variables.h

SOURCES += src/edible.c \
           src/gtr.c \
           src/llh.c \
           src/matrix.c \
           src/new_models.c \
           src/options.c \
           src/partial.c \
           src/random.c \
           src/read.c \
           src/tree.c \
           src/utility.c \
           src/pcalc/pcalc.c
