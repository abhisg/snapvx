all: libSolver.a

CFLAGS += -I../CVXCanon/include/Eigen/ -std=c++0x -O2 -DCTRLC=1 -fPIC -x c++

SOLVER_OBJS = Edge.o ModSquare.o NetLasso.o Node.o ProximalMap.o Square.o
libSolver.a: $(SOLVER_OBJS)
	ar rcs $@ $^
	ranlib $@


%.o: %.hpp
	g++ $(CFLAGS) -c $< -o $@

.PHONY:clean
clean: 
	rm -rf *.o *.a
