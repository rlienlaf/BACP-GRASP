#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>     	
#include <ctime>
#include <cmath>
#include <algorithm>
using namespace std;

int main(int argc, char* argv[]){

	float tiempo, brokens, eval, p, restart, iter, len;
	float eval_min, eval_max, eval_av, tiempo_min, tiempo_max, tiempo_av, brokens_min;
	string aux_a, aux_b;
	string files[4] = {"toy.txt","bacp8-2.txt","bacp10-2.txt","bacp12-2.txt"};
  string line;
  string file;

  ofstream writefile ("tables.txt");
  if (writefile.is_open())
  {
  	for (int m = 0; m < 4 ; ++m)
		{

  			for (int i = 0; i < 2+1; ++i)
  			{
  				for (int j = 10; j < 51; j+=40)
  				{
  					for (int k = 200; k < 1001; k+=800)
  					{
  						for (int l = 1; l < 3+1 ; ++l)
  						{
  							file = "results/results_"+files[m]+"_"+to_string(i)+"_"+to_string(j)+"_"+to_string(k)+"_"+to_string(l)+".txt";
  							ifstream myfile (file);
  							if (myfile.is_open())
  							{
  								eval_min=100000;
  								eval_max=0;
  								eval_av=0; 
  								tiempo_min=100000;
  								tiempo_max=0;
  								tiempo_av = 0;
  								brokens_min=100000;
  								for (int n = 0; n < 20; ++n)
  								{
  									myfile >> tiempo >> brokens >> eval >> p >> restart >> iter >> len;
  									eval_min=min(eval_min,eval);
  									eval_max=max(eval_max,eval);
  									eval_av+=eval;
  									tiempo_min=min(tiempo_min,tiempo);
  									tiempo_max=max(tiempo_max,tiempo);
  									tiempo_av+=tiempo;
  									brokens_min=min(brokens_min,brokens);
  								}
  								writefile << files[m] << " & " << p << " & " << restart << " & " << iter << " & " << len << " & " << tiempo_min << " & " << tiempo_av/20 << " & " << tiempo_max << " & " << brokens_min << " & " << eval_min << " & " << eval_av/20 << " & " << eval_max << "\\\\ \n"; 
  								myfile.close();
  							}
  							else cout << "Unable to open file"; 
  						}
  					}
  				}
  			}
  			
		}	
    writefile.close();
  }
  else cout << "Unable to open file";
}
  