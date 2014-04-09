simulation.run: LibSimulation.so main.C
	`root-config --cxx --cflags --libs` -O2 -o $@ $^

LibSimulation.so: code/*.C code/*.h
	`root-config --cxx --cflags --libs` -O2 --shared -o $@ code/*.C

plots: simulation.run scan.root plot.C
	./simulation.run
	root -l scan.root plot.C

plot.C: