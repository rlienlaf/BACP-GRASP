#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>     	
#include <ctime>
#include <cmath>
#include <algorithm>
int main(int argc, char* argv[]){
	clock_t tStart = clock();
	int random_seed = time(NULL);
	if(argc == 7){
		random_seed = atoi(argv[6]);
	}
	int sum = 0;
	srand (random_seed);
	int restart = atoi(argv[3]);
	int max_iter = atoi(argv[4]);
	p = atoi(argv[2]);

	Course* courses;
	Precedence* precedences;
	
	//Leer los datos
	string aux_a, aux_b;
  string line, filename = argv[1];
  ifstream myfile (filename);
  if (myfile.is_open())
  {
  	myfile >> years;
  	myfile >> num_periods;
  	myfile >> min_load >> max_load;
  	myfile >> min_courses >> max_courses;
  	myfile >> num_courses;
  	myfile >> num_precedences;
  	courses = new Course[num_courses];
		precedences = new Precedence[num_precedences];
  	for (int i = 0; i < num_courses; ++i)
  	{
  		myfile >> courses[i].name >> courses[i].credits;
  		sum += courses[i].credits;
  	}
  	for (int i = 0; i < num_precedences; ++i)
  	{
  		myfile >> aux_b >> aux_a;
  		for (int j = 0; j < num_courses; ++j)
  		{
  			if (courses[j].name == aux_a)
  			{
  				precedences[i].first = courses[j];
  			}
  			if (courses[j].name == aux_b)
  			{
  				precedences[i].second = courses[j];
  			}
  		}
  	}
    myfile.close();
  }
  else cout << "Unable to open file"; 