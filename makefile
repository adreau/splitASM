molsplit: molsplit.o mol_stat.o stat_per_interval.o outliers_det.o create_bed_file.o
	g++ -std=c++17 -o splitASM molsplit.o mol_stat.o stat_per_interval.o outliers_det.o create_bed_file.o -lstdc++fs -pthread

molsplit.o: molsplit.cpp mol_stat.h stat_per_interval.h outliers_det.h
	g++ -std=c++17 -c molsplit.cpp -lstdc++fs -pthread

mol_stat.o: mol_stat.cpp mol_stat.h
	g++ -std=c++17 -c mol_stat.cpp

stat_per_interval.o: stat_per_interval.cpp stat_per_interval.h
	g++ -std=c++17 -c stat_per_interval.cpp

outliers_det.o: outliers_det.cpp outliers_det.h
		g++ -std=c++17 -c outliers_det.cpp

create_bed_file.o: create_bed_file.cpp create_bed_file.h
	  g++ -std=c++17 -c create_bed_file.cpp

clean:
	rm -f *~ *.o
