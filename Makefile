CC= g++ 

EXE1=gstoreu
EXE2=gstored

#CCFLAGS+= -g -Wall -march=native -fopenmp  -DGENEARATOR_USE_PACKED_TYPE 
#CCFLAGS+= -g -Wall -march=native -fopenmp  
#CCFLAGS = -O3 -Wall -march=native -fgcse-sm -fgcse-las -fgcse-after-reload -floop-strip-mine -ftree-loop-im -fivopts -funswitch-loops -fopenmp -DGENEARATOR_USE_PACKED_TYPE 
CCFLAGS = -O3 -Wall -march=native -fgcse-sm -fgcse-las -fgcse-after-reload -floop-strip-mine -ftree-loop-im -fivopts -funswitch-loops -fopenmp 

#OBJS=  	grid.o util.o pr.o bfs.o  wcc.o kcore.o wcc2.o bfs2.o
DEPS=gstore.cpp\
		gstore.h\
		pr.cpp\
		bfs.cpp\
		bfs2.cpp\
		kcore.cpp\
		util.cpp\
		wcc.cpp\
		wcc.h\
		algo.h\
		kcore.h\
		pr.h\
		bfs.h\
		bfs2.h\
		Makefile

SRC=gstore.cpp\
		pr.cpp\
		bfs.cpp\
		bfs2.cpp\
		kcore.cpp\
		util.cpp\
		wcc.cpp\
	

${EXE1}: $(DEPS)
	${CC} $(CCFLAGS) -DHALF_GRID ${SRC} -o ${EXE1} -laio

${EXE2}: $(DEPS)
	${CC} $(CCFLAGS) ${SRC} -o ${EXE2} -laio
	
clean:
	rm -rf *.o ${EXE1} ${EXE2}

