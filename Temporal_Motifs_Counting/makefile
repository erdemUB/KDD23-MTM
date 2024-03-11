tmc: main.o tmc.o
	g++ -std=c++11 main.o tmc.o -o TMC

main.o: main.cpp
	g++ -c -std=c++11 main.cpp

tmc.o: tmc.cpp
	g++ -c -std=c++11 tmc.cpp

clean:
	rm *.o TMC