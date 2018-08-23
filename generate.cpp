#include <iostream>
#include <fstream>
#include <cstdlib>     	
#include <ctime>
using namespace std;

int main(){
	ofstream myfile ("seeds.txt");
  if (myfile.is_open())
  {
  	for (int i = 0; i < 100 ; ++i)
		{
			myfile << time(NULL)+i << "\n";
		}	
    myfile.close();
  }
  else cout << "Unable to open file";

}