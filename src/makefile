mtm: main.o mtm.o
	g++ -std=c++11 main.o mtm.o -o MTM

main.o: main.cpp
	g++ -c -std=c++11 main.cpp

mtm.o: mtm.cpp
	g++ -c -std=c++11 mtm.cpp

clean:
	rm *.o MTM
