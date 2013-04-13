CPPUNIT_INCLUDE = 
CPPUNIT_LIB = 


all_tests: meshtest
	@echo "Start all tests!"

meshtest: all_tests.cpp ../src/mesh.cpp 
	g++ -L$(CPPUNIT_LIB) -I./include/ -I$(CPPUNIT_INCLUDE)\
	 	all_tests.cpp ../src/mesh.cpp -o mesh_test.exe && ./mesh_test.exe


