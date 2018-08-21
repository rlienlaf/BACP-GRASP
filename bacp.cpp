#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>     /* srand, rand */	
#include <ctime>
#include <algorithm>
#include <cmath>
using namespace std;

int length;
int years, periods, num_periods, min_courses, max_courses, min_load, max_load, num_courses, num_precedences, p;
float average = 0;

struct Course{
	string name;
	int credits;
};

struct Precedence{
	Course first;
	Course second;
};

struct Solution{
	int period;
	Course course;
};

struct Load{
	int	period;
	int total;
};

struct Candidate{
	int period;
	int brokens[3] = {0};
	int broken_constraints;
	float eval;
};


//Solution* buildSolution(Course* cou, Precedence* pre);
bool sortLoadByTotal(const Load &left, const Load &right) { return left.total < right.total;}
bool sortLoadByPeriod(const Load &left, const Load &right) { return left.period < right.period;}
bool sortCandidateByBrokenConstraints(const Candidate &left, const Candidate &right) { return left.broken_constraints < right.broken_constraints;}
bool sortCandidateByEvaluation(const Candidate &left, const Candidate &right) { return left.eval < right.eval;}
bool sortSolutionByPeriod(const Solution &left, const Solution &right) { return left.period < right.period;}


void setToZero(int* array, int len){
	for (int i = 0; i < len; ++i)
	{
		array[i] = 0;
	}
}

void countLoadPerPeriod(Solution* solution_to_count, int* load_per_period){
	//cout << "begin countLoadPerPeriod \n";
	//cout << "for \n";
	for (int i = 0; i < num_courses; ++i)
	{
		load_per_period[solution_to_count[i].period]+=solution_to_count[i].course.credits;
	}
	//cout << "return countLoadPerPeriod \n";
}

void countCoursesPerPeriod(Solution* solution_to_count, int* courses_per_period){
	//cout << "begin countCoursesPerPeriod \n";
	for (int i = 0; i < num_courses; ++i)
	{
		courses_per_period[solution_to_count[i].period]+=1;
	}
	//cout << "return countCoursesPerPeriod \n";
}

void showSolution(Solution* solution_to_print){
	int total_load[periods]={0};
	int total_courses[periods]={0};
	int count = 0;
	Solution* aux = new Solution[num_courses];
	copy(solution_to_print, solution_to_print + num_courses, aux);
	
	//cout << "call countLoadPerPeriod \n";
	countLoadPerPeriod(aux, total_load);
	//cout << "call countCoursesPerPeriod \n";
	countCoursesPerPeriod(aux, total_courses);

	int* max = max_element(total_courses, total_courses + periods);

	sort(aux, aux+num_courses, sortSolutionByPeriod);

	cout << "Solution: \n";
	for (int i = 0; i < periods; ++i)
	{
		cout << i << "\t|";
		for (int j = 0; j < *max ; ++j)
		{
			if(aux[count].period == i){
				cout << aux[count].course.name << "\t";	
				count++;
			}
			else{
				cout << "\t";
			}
		}
		cout << "| " << total_load[i] << "\n"; 
	}

}

int countBrokenPrecedences(Solution* sol, Precedence* pre){
	int broken = 0, first_period, second_period;
	for (int i = 0; i < num_precedences; ++i)
	{
		for (int j = 0; j < num_courses; ++j)
		{
			if (sol[j].course.name == pre[i].first.name)
			{
				first_period = sol[j].period;
			}
		}
		for (int j = 0; j < num_courses; ++j)
		{
			if (sol[j].course.name == pre[i].second.name)
			{
				second_period = sol[j].period;
			}
		}
		if (second_period <= first_period)
		{
			broken++;
		}
	}
	return broken;
}

int countBrokenLoad(int* load_per_period){
	int broken = 0;
	for (int i = 0; i < periods; ++i)
	{
		if (load_per_period[i] < min_load || load_per_period[i] > max_load)
		{
			broken++;
		}
	}
	return broken;
}


int countBrokenCourses(int* courses_per_period){
	int broken = 0;
	for (int i = 0; i < periods; ++i)
	{
		if (courses_per_period[i] < min_courses || courses_per_period[i] > max_courses)
		{
			broken++;
		}
	}
	return broken;
}


float calculateAverage (Course* c){
	float count = 0; 
	for (int i = 0; i < num_courses; ++i)
	{
		count += c[i].credits;
	}
	return (count/periods);
}


void showBrokens(Solution* sol, Precedence* pre){
	int total_load_aux[periods]={0};
	int total_courses_aux[periods]={0};
	countLoadPerPeriod(sol, total_load_aux);
	countCoursesPerPeriod(sol, total_courses_aux);
	int broken_load = countBrokenLoad(total_load_aux);
	int broken_courses = countBrokenCourses(total_courses_aux);
	int broken_precedences = countBrokenPrecedences(sol, pre);
	cout << "broken loads: " << broken_load << "\n";
	cout << "broken courses: " << broken_courses << "\n";
	cout << "broken precedences: " << broken_precedences << "\n";
}


float evaluateSolution(Solution* solution_to_evaluate){
	float evaluation = 0;
	int total_load_aux[periods]={0};
	//int total_courses_aux[periods]={0};
	countLoadPerPeriod(solution_to_evaluate, total_load_aux);
	//countCoursesPerPeriod(solution_to_evaluate, total_courses_aux);
	if(p == 0)
	{
		int* max = max_element(total_load_aux, total_load_aux + periods);
		evaluation = *max;
	}
	else
	{
		for (int i = 0; i < periods; ++i)
		{
			evaluation += pow(abs((total_load_aux[i] - average)), p);
		}
	}
	return evaluation;
}

void hillClimbing(Solution* best, Solution* initial, Precedence* pre, int iter){
	int hc_total_load[periods]={0};
	int hc_total_courses[periods]={0};
	int elem;
	int best_candidate;

	//crear alternative_solution = best_solution
	Solution* alternative_solution = new Solution[num_courses];
	cout << "HC-initial: \n";
	showSolution(initial);
	showBrokens(initial, pre);
	cout << "evaluation: " << evaluateSolution(initial) << "\n";
	copy(initial, initial+num_courses,alternative_solution);

	//*alternative_solution = *initial;
	/*cout << "HC-alternative: \n";
	showSolution(alternative_solution);
	showBrokens(alternative_solution, pre);
	*/

	Candidate* candidate_periods = new Candidate[periods];
	
	for (int j = 0; j < iter; ++j)
	{
		//elegir elemento de alternative_solution
		elem = rand()%num_courses;
		cout << "elem: " << elem << "\n";
		best_candidate = initial[elem].period;
		cout << "initial - name: " << initial[elem].course.name << "\n";
		cout << "initial - credits: " << initial[elem].course.credits << "\n";
		cout << "initial - period: " << initial[elem].period << "\n";
		cout << "best_candidate: " << best_candidate << "\n";
		//cout << "best candidate1 " << best_candidate;
		Candidate initial_candidate;
		//modificar alternative_solution
		for (int i = 0; i < periods; ++i)
		{
			
			alternative_solution[elem].period = i;
			cout << "movement: " << alternative_solution[elem].course.name << " " << alternative_solution[elem].period << "\n";

			setToZero(hc_total_load, periods);
			setToZero(hc_total_courses, periods);
			showSolution(alternative_solution);
			showBrokens(alternative_solution, pre);
			cout << "evaluation: " << evaluateSolution(alternative_solution) << "\n";
			candidate_periods[i].period = i;

			countLoadPerPeriod(alternative_solution, hc_total_load);
			countCoursesPerPeriod(alternative_solution, hc_total_courses);
			for (int k = 0; k < periods; ++k)
			{
				cout << hc_total_load[k] <<" ";
			}
			//cout << "\n";
			candidate_periods[i].brokens[0] = countBrokenLoad(hc_total_load);
			//cout << "b0: " <<candidate_periods[i].brokens[0]<< "\n";
			candidate_periods[i].brokens[1] = countBrokenCourses(hc_total_courses);
			//cout << "b1: " <<candidate_periods[i].brokens[1]<< "\n";
			candidate_periods[i].brokens[2] = countBrokenPrecedences(alternative_solution, pre);
			//cout << "b2: " <<candidate_periods[i].brokens[2]<< "\n";
			candidate_periods[i].broken_constraints = candidate_periods[i].brokens[0] + candidate_periods[i].brokens[1] + candidate_periods[i].brokens[2];
			//cout << "bt: " <<candidate_periods[i].broken_constraints<< "\n";
			candidate_periods[i].eval = evaluateSolution(alternative_solution);
			cout << "\ne: " <<candidate_periods[i].eval<< "\n";
			if (i==best_candidate)
			{
				//copy(candidate_periods[i], candidate_periods[i] + sizeof(struct Candidate),initial_candidate);
				initial_candidate = candidate_periods[i];
			}
		}

		sort(candidate_periods, candidate_periods+periods, sortCandidateByBrokenConstraints);

		cout << "candidates: \n";
		for (int i = 0; i < periods; ++i)
		{
			cout << candidate_periods[i].period << " " << candidate_periods[i].broken_constraints << "\n";
		}

		cout << "periodos: " << initial_candidate.period << " -> " << candidate_periods[0].period << "\n";
		cout << "broken: " << initial_candidate.broken_constraints << " -> " << candidate_periods[0].broken_constraints << "\n";
		


		// DE AQUI PA ABAJO ESTA TODO MALO
		if (candidate_periods[0].broken_constraints <= initial_candidate.broken_constraints)
		{
			if (candidate_periods[0].broken_constraints == 0)
			{
				if(evaluateSolution(alternative_solution) <= evaluateSolution(initial))
				{
					alternative_solution[elem].period = candidate_periods[0].period;
					copy(alternative_solution, alternative_solution+num_courses,initial);
				}
			}
			else 
			{
				alternative_solution[elem].period = candidate_periods[0].period;
				copy(alternative_solution, alternative_solution+num_courses,initial);
			}
		}
		//cout << " " << best_candidate << "\n";
	}



  setToZero(hc_total_load, periods);
  setToZero(hc_total_courses, periods);
	countLoadPerPeriod(best, hc_total_load);
	countCoursesPerPeriod(best, hc_total_courses);
	int best_solution_broken_constraints = countBrokenLoad(hc_total_load) + countBrokenCourses(hc_total_courses) + countBrokenPrecedences(best, pre);
	cout << "best brokens:" << best_solution_broken_constraints << "\n";

  setToZero(hc_total_load, periods);
  setToZero(hc_total_courses, periods);
	countLoadPerPeriod(initial, hc_total_load);
	countCoursesPerPeriod(initial, hc_total_courses);
	int initial_solution_broken_constraints = countBrokenLoad(hc_total_load) + countBrokenCourses(hc_total_courses) + countBrokenPrecedences(initial, pre);
	cout << "initial brokens: " << initial_solution_broken_constraints << "\n";


	cout << "\nThis round best: \n";
	showSolution(initial);
	showBrokens(initial, pre);
	cout << "evaluation: " << evaluateSolution(initial) << "\n";


	cout << "\nGlobal best: \n";
	showSolution(best);
	showBrokens(best, pre);
	cout << "evaluation: " << evaluateSolution(best) << "\n";
	if (initial_solution_broken_constraints <= best_solution_broken_constraints)
	{
		if (initial_solution_broken_constraints == 0)
		{
			if (evaluateSolution(initial) <= evaluateSolution(best))
			{
				copy(initial, initial + num_courses, best);
			}
		}
		else 
		{
			copy(initial, initial + num_courses, best);
		}
	}
	cout << "\nNew global best: \n";
	showSolution(best);
	showBrokens(best, pre);
	cout << "evaluation: " << evaluateSolution(best) << "\n";
}



Solution* buildSolution(Course* courses, Precedence* precedences){
	int rcl[length] = {0};
	int chosen;
	Load load[periods];
	Solution* posible_solution = new Solution[num_courses];
	//fill load
	for (int i = 0; i < periods; ++i)
	{
		load[i].period = i;
		load[i].total = 0;
	}

	for (int i = 0; i < num_courses; ++i)
	{
		posible_solution[i].course = courses[i];
		//update load
		for (int j = 0; j < periods; ++j)
		{
			load[j].total += courses[i].credits;
		}
		sort(load, load+periods, sortLoadByTotal);
		for (int j = 0; j < length; ++j)
		{
			rcl[j] = load[j].period;
		}
		chosen = rand()%length;
		posible_solution[i].period = rcl[chosen];
		for (int j = 0; j < periods; ++j)
		{
			if(load[j].period!=rcl[chosen]){
				load[j].total -= courses[i].credits;
			}
		}
	}
	/*
	cout << "solution: \n";
	for (int i = 0; i < num_courses; ++i)
  {
  	cout << posible_solution[i].period << ' ' << posible_solution[i].course.name << ' ' << posible_solution[i].course.credits <<"\n";
  }

  sort(load, load+periods, sortLoadByPeriod);
  for (int i = 0; i < periods; ++i)
  {
  	cout << "periodo: " << load[i].period << "\ttotal: " << load[i].total << '\n';
  }*/

	//cout << "return solution \n";
	return posible_solution;
}




int main(int argc, char* argv[]){

	int random_seed = time(NULL);
	int sum = 0;
	cout << random_seed << '\n';
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
  cout << "===================================START===================================\n";
  if (myfile.is_open())
  {
  	myfile >> years;
  	myfile >> num_periods;
  	myfile >> min_load >> max_load;
  	myfile >> min_courses >> max_courses;
  	myfile >> num_courses;
  	myfile >> num_precedences;
  	
  	cout << years << '\n';
  	cout << num_periods << '\n';
  	cout << min_load << " " << max_load << '\n';
  	cout << min_courses << " " << max_courses << '\n';
  	cout << num_courses << '\n';
  	cout << num_precedences << '\n';

  	courses = new Course[num_courses];
		precedences = new Precedence[num_precedences];

  	for (int i = 0; i < num_courses; ++i)
  	{
  		myfile >> courses[i].name >> courses[i].credits;
  		sum += courses[i].credits;
  		cout << courses[i].name << ' ' << courses[i].credits << '\n';
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
  		cout << precedences[i].first.name << ' ' << precedences[i].second.name << '\n' ;
  	}
    myfile.close();
  }
  else cout << "Unable to open file"; 

  periods = years*num_periods;
	length = (atoi(argv[5])<periods)?atoi(argv[5]):periods;
  average = calculateAverage(courses);

  cout << "total credits : " << sum <<'\n';
  cout << "periods : " << periods << "\n";
	cout << "average : " << average << "\n";

  /*Solution* initial_solution = buildSolution(courses, precedences);
  cout << "\nInitial solution: \n";
  showSolution(initial_solution);
  showBrokens(initial_solution, precedences);
*/
  Solution* best_solution = buildSolution(courses, precedences);
  Solution* initial_solution;
  cout << "=================================\n";
  cout << "\nBest solution: \n";
  showSolution(best_solution);
  showBrokens(best_solution, precedences);
  cout << "evaluation: " << evaluateSolution(best_solution) << "\n";
  cout << "=================================\n";


  for (int i = 0; i < restart; ++i)
  {
  	initial_solution = buildSolution(courses, precedences);
  	cout << "-------------------------------------------------------------------------\n"; 
  	cout << "Initial solution iteration " << i << "\n";
  	showSolution(initial_solution);
  	showBrokens(initial_solution, precedences);
  	cout << "evaluation: " << evaluateSolution(initial_solution) << "\n";
  	hillClimbing(best_solution, initial_solution, precedences, max_iter);
  }
	cout << "=================================\n";
  cout << "\nBest solution: \n";
  showSolution(best_solution);
  showBrokens(best_solution, precedences);
  cout << "evaluation: " << evaluateSolution(best_solution) << "\n";
  cout << "=================================\n";
  return 0;
}


