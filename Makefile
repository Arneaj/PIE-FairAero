file := path_to_petsc.config
export PETSC_DIR := $(shell cat ${file})

export PETSC_ARCH=arch-linux-c-debug

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

basic: src/basic.cpp
	$(info    PETSC_DIR is $(PETSC_DIR))
	${CXX} ${PETSC_CC_INCLUDES} src/basic.cpp -o basic.out ${PETSC_LIB}
	$(info    Compilation done)

grid: src/grid.cpp
	$(info    PETSC_DIR is $(PETSC_DIR))
	${CXX} ${PETSC_CC_INCLUDES} src/grid.cpp -o grid.out ${PETSC_LIB}
	$(info    Compilation done)

clean::
	rm *.out

