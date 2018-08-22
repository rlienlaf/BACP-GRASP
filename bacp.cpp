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
	float total;
	int flag = 0;
};

struct Candidate{
	int period;
	float brokens[3] = {0};
	float broken_constraints;
	float eval;
};

struct PreEvaluation{
	Course course;
	float value = 0;
};


//Solution* buildSolution(Course* cou, Precedence* pre);
bool sortLoadByTotal(const Load &left, const Load &right) { return left.total < right.total;}
bool sortLoadByPeriod(const Load &left, const Load &right) { return left.period < right.period;}
bool sortLoadByFlag(const Load &left, const Load &right){ return left.flag < right.flag;}

bool sortCandidateByBrokenConstraints(const Candidate &left, const Candidate &right) { return left.broken_constraints < right.broken_constraints;}
bool sortCandidateByEvaluation(const Candidate &left, const Candidate &right) { return left.eval < right.eval;}
bool sortSolutionByPeriod(const Solution &left, const Solution &right) { return left.period < right.period;}
bool sortPreProcessingByEval(const PreEvaluation &left, const PreEvaluation &right) { return left.value < right.value;}


void setToZero(int* array, int len){
	for (int i = 0; i < len; ++i)
	{
		array[i] = 0;
	}
}

void preProcessing(Course* co, Precedence* pre){
	PreEvaluation* eval = new PreEvaluation[num_courses];
	for (int i = 0; i < num_courses; ++i)
	{
		eval[i].course = co[i];
		//cout << "course: " << co[i].name << "\n";
		for (int j = 0; j < num_precedences; ++j)
		{
			if (pre[j].first.name == eval[i].course.name)
			{
				eval[i].value -=1;
			}
			if (pre[j].second.name == eval[i].course.name)
			{
				eval[i].value +=1;
			}
		}
	}
	sort(eval, eval + num_courses, sortPreProcessingByEval);
	for (int i = 0; i < num_courses; ++i)
	{
		co[i] = eval[i].course;
		//cout << "new course: " << co[i].name << "\n";
	}
}

void countLoadPerPeriod(Solution* solution_to_count, int* load_per_period, int stop = num_courses){
	//cout << "begin countLoadPerPeriod \n";
	//cout << "for \n";
	for (int i = 0; i < stop; ++i)
	{
		//cout <<" " << solution_to_count[i].period << "->" << load_per_period[solution_to_count[i].period] <<" ";
		load_per_period[solution_to_count[i].period]+=solution_to_count[i].course.credits;
		//cout <<" " << solution_to_count[i].period << "->" << load_per_period[solution_to_count[i].period] <<" ";
	}
	//cout <<" countedload\n";
	//cout << "return countLoadPerPeriod \n";
}

void countCoursesPerPeriod(Solution* solution_to_count, int* courses_per_period, int stop = num_courses){
	//cout << "begin countCoursesPerPeriod \n";
	for (int i = 0; i < stop; ++i)
	{
		courses_per_period[solution_to_count[i].period]+=1;
	}
	//cout << "return countCoursesPerPeriod \n";
}

void showSolution(Solution* solution_to_print, int stop = num_courses){
	int total_load[periods]={0};
	int total_courses[periods]={0};
	int count = 0;
	Solution* aux = new Solution[num_courses];
	copy(solution_to_print, solution_to_print + stop, aux);
	
	//cout << "call countLoadPerPeriod \n";
	countLoadPerPeriod(aux, total_load, stop);
	//cout << "call countCoursesPerPeriod \n";
	countCoursesPerPeriod(aux, total_courses, stop);

	int* max = max_element(total_courses, total_courses + periods);

	sort(aux, aux+stop, sortSolutionByPeriod);

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
//delete(aux);
}


int countBrokenPrecedences(Solution* sol, Precedence* pre, int stop = num_courses){
	int broken = 0, first_period=-1, second_period=periods+1;
	for (int i = 0; i < num_precedences; ++i)
	{
		first_period=-1, second_period=periods+1;
		for (int j = 0; j < stop; ++j)
		{
			if (sol[j].course.name == pre[i].first.name)
			{
				first_period = sol[j].period;
			}
		}
		for (int j = 0; j < stop; ++j)
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



float countBrokenLoad(int* load_per_period, int stop = periods){
	float broken = 0;
	for (int i = 0; i < stop; ++i)
	{
		if (load_per_period[i] < min_load || load_per_period[i] > max_load)
		{
			broken += pow(abs((load_per_period[i] - average)), 1);
			//broken++;
		}
	}
	return broken;
}


float countBrokenCourses(int* courses_per_period, int stop = periods){
	float broken = 0;
	float average_courses = num_courses/periods;
	for (int i = 0; i < stop; ++i)
	{
		if (courses_per_period[i] < min_courses || courses_per_period[i] > max_courses)
		{
			broken += pow(abs((courses_per_period[i] - average_courses)), 1);
			//broken++;
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
	float broken_load = countBrokenLoad(total_load_aux);
	float broken_courses = countBrokenCourses(total_courses_aux);
	float broken_precedences = countBrokenPrecedences(sol, pre);
	cout << "broken loads: " << broken_load << "\n";
	cout << "broken courses: " << broken_courses << "\n";
	cout << "broken precedences: " << broken_precedences << "\n";
}


float evaluateSolution(Solution* solution_to_evaluate, int stop = num_courses){
	float evaluation = 0;
	int total_load_aux[periods]={0};
	//int total_courses_aux[periods]={0};
	countLoadPerPeriod(solution_to_evaluate, total_load_aux, stop);
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
	int elem, elem2;
	int best_candidate;
	int aux;
	int counter;

	//crear alternative_solution = best_solution
	Solution* alternative_solution = new Solution[num_courses];
	//cout << "HC-initial: \n";
	//showSolution(initial);
	//showBrokens(initial, pre);
	//cout << "evaluation: " << evaluateSolution(initial) << "\n";
	copy(initial, initial+num_courses,alternative_solution);

	//*alternative_solution = *initial;
	/*cout << "HC-alternative: \n";
	showSolution(alternative_solution);
	showBrokens(alternative_solution, pre);
	*/

	Candidate* candidate_periods = new Candidate[periods];
	Candidate initial_candidate;
	
	for (int j = 0; j < iter; ++j)
	{
		//elegir elemento de alternative_solution
		elem = rand()%num_courses;
		//elem2 = rand()%num_courses;
		
		//cout << "elem: " << elem << "\n";
		best_candidate = initial[elem].period;
		//cout << "initial - name: " << initial[elem].course.name << "\n";
		//cout << "initial - credits: " << initial[elem].course.credits << "\n";
		//cout << "initial - period: " << initial[elem].period << "\n";
		//cout << "best_candidate: " << best_candidate << "\n";
		



		//modificar alternative_solution
		for (int i = 0; i < periods; ++i)
		{
			
			alternative_solution[elem].period = i;
			//cout << "movement: " << alternative_solution[elem].course.name << " " << alternative_solution[elem].period << "\n";

			setToZero(hc_total_load, periods);
			setToZero(hc_total_courses, periods);
			//showSolution(alternative_solution);
			//showBrokens(alternative_solution, pre);
			//cout << "evaluation: " << evaluateSolution(alternative_solution) << "\n";
			candidate_periods[i].period = i;

			countLoadPerPeriod(alternative_solution, hc_total_load);
			countCoursesPerPeriod(alternative_solution, hc_total_courses);
			/*for (int k = 0; k < periods; ++k)
			{
				cout << hc_total_load[k] <<" ";
			}*/
			candidate_periods[i].brokens[0] = countBrokenLoad(hc_total_load);
			candidate_periods[i].brokens[1] = countBrokenCourses(hc_total_courses);
			candidate_periods[i].brokens[2] = countBrokenPrecedences(alternative_solution, pre);
			candidate_periods[i].broken_constraints = candidate_periods[i].brokens[0] + candidate_periods[i].brokens[1] + candidate_periods[i].brokens[2];
			candidate_periods[i].eval = evaluateSolution(alternative_solution);
			//cout << "\ne: " <<candidate_periods[i].eval<< "\n";
			if (i==best_candidate)
			{
				//copy(candidate_periods[i], candidate_periods[i] + sizeof(struct Candidate),initial_candidate);
				initial_candidate = candidate_periods[i];
			}
		}

		/*		for (int i = 0; i < periods; ++i)
		{
			cout << candidate_periods[i].period << " " << candidate_periods[i].broken_constraints << " " << candidate_periods[i].eval << "\n";
		}*/

		sort(candidate_periods, candidate_periods+periods, sortCandidateByBrokenConstraints);
		counter = 0;
 /*
		for (int i = 0; i < periods; ++i)
		{
			cout << candidate_periods[i].period << " " << candidate_periods[i].broken_constraints << " " << candidate_periods[i].eval << "\n";
		}
		*/
		for (int i = 0; i < periods; ++i)
		{
			if(candidate_periods[i].broken_constraints == candidate_periods[0].broken_constraints)
			{
				counter++;
			}
		}
		//cout << "-> " << counter << "\n";

		sort(candidate_periods, candidate_periods + counter, sortCandidateByEvaluation);
		
/*
		for (int i = 0; i < periods; ++i)
		{
			cout << candidate_periods[i].period << " " << candidate_periods[i].broken_constraints << " " << candidate_periods[i].eval << "\n";
		}*/


		/*cout << "candidates: \n";
		for (int i = 0; i < periods; ++i)
		{
			cout << candidate_periods[i].period << " " << candidate_periods[i].broken_constraints << "\n";
		}*/

		//cout << "periodos: " << initial_candidate.period << " -> " << candidate_periods[0].period << "\n";
		//cout << "broken: " << initial_candidate.broken_constraints << " -> " << candidate_periods[0].broken_constraints << "\n";
		
		//Metodo 3
/*
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
		*/

		// Metodo 4
		//cout << candidate_periods[0].broken_constraints << " " << initial_candidate.broken_constraints << "\n";
		
		if (candidate_periods[0].broken_constraints < initial_candidate.broken_constraints)
		{
			//cout << "if\n";
			alternative_solution[elem].period = candidate_periods[0].period;
			copy(alternative_solution, alternative_solution + num_courses, initial);
		}
		else if (candidate_periods[0].broken_constraints == initial_candidate.broken_constraints)
		{
			//cout << "else\n";
			aux = alternative_solution[elem].period;
			alternative_solution[elem].period = candidate_periods[0].period;
			//cout << evaluateSolution(alternative_solution) << " " << evaluateSolution(initial) << "\n";
			if (evaluateSolution(alternative_solution) <= evaluateSolution(initial))
			{
				//cout << "else if\n";
				//alternative_solution[elem].period = candidate_periods[0].period;
				copy(alternative_solution, alternative_solution + num_courses, initial);
			}
			else
			{
				//cout << "else else\n";
				alternative_solution[elem].period = aux;
			}
		}

		//cout << " " << best_candidate << "\n";
	}



  setToZero(hc_total_load, periods);
  setToZero(hc_total_courses, periods);
	countLoadPerPeriod(best, hc_total_load);
	countCoursesPerPeriod(best, hc_total_courses);
	float best_solution_broken_constraints = countBrokenLoad(hc_total_load) + countBrokenCourses(hc_total_courses) + countBrokenPrecedences(best, pre);
	cout << "best brokens:" << best_solution_broken_constraints << "\n";

  setToZero(hc_total_load, periods);
  setToZero(hc_total_courses, periods);
	countLoadPerPeriod(initial, hc_total_load);
	countCoursesPerPeriod(initial, hc_total_courses);
	float initial_solution_broken_constraints = countBrokenLoad(hc_total_load) + countBrokenCourses(hc_total_courses) + countBrokenPrecedences(initial, pre);
	cout << "initial brokens: " << initial_solution_broken_constraints << "\n";

/*
	cout << "\nThis round best: \n";
	showSolution(initial);
	showBrokens(initial, pre);
	cout << "evaluation: " << evaluateSolution(initial) << "\n";


	cout << "\nGlobal best: \n";
	showSolution(best);
	showBrokens(best, pre);
	cout << "evaluation: " << evaluateSolution(best) << "\n";
*/

//metodo 1
	/*
	if (evaluateSolution(initial) <= evaluateSolution(best))
	{
		copy(initial, initial + num_courses, best);
	}*/
//metodo 2
	/*
	if (initial_solution_broken_constraints <= best_solution_broken_constraints)
	{

			if (evaluateSolution(initial) <= evaluateSolution(best))
			{
				copy(initial, initial + num_courses, best);
			}
	}*/
//metodo 3
/*
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
	}*/

	//metodo 4
	if (initial_solution_broken_constraints < best_solution_broken_constraints)
	{
		copy(initial, initial + num_courses, best);
	}
	else if (initial_solution_broken_constraints == best_solution_broken_constraints)
	{
			if (evaluateSolution(initial) <= evaluateSolution(best))
			{
				copy(initial, initial + num_courses, best);
			}
	}

	cout << "\nNew global best: \n";
	showSolution(best);
	showBrokens(best, pre);
	cout << "evaluation: " << evaluateSolution(best) << "\n";

	//delete(candidate_periods);
	//delete(alternative_solution);
}



Solution* buildSolution(Course* courses, Precedence* pre){
	int bs_total_load[periods]={0};
	int bs_total_courses[periods]={0};
	int rcl[length] = {0};
	int chosen;
	int min_period = -1;
	Load load[periods] = {0};
	Solution* posible_solution = new Solution[num_courses];
	//fill load
	
	for (int j = 0; j < periods; ++j)
		{
			load[j].period = j;
			load[j].total = 0;
			load[j].flag = 0;
		}

	for (int i = 0; i < num_courses; ++i)
	{
		for (int j = 0; j < periods; ++j)
		{
			load[j].period = j;
			load[j].total = 0;
			load[j].flag = 0;
		}
		posible_solution[i].course = courses[i];
		//update load

		min_period = -1;
		for (int j = 0; j < num_precedences; ++j)
		{
			if (posible_solution[i].course.name == pre[j].second.name)
			{
				for (int k = 0; k < i+1; ++k)
				{	
					if (posible_solution[k].course.name == pre[j].first.name)
					{
						min_period = max(min_period, posible_solution[k].period);
					}
				}	
			}
		}


		for (int j = 0; j < periods; ++j)
		{
					posible_solution[i].period = load[j].period;
					setToZero(bs_total_load, periods);
					setToZero(bs_total_courses, periods);
					/*for (int l = 0; l < periods; ++l)
					{
						cout << bs_total_load[l];
					}*/
					countLoadPerPeriod(posible_solution, bs_total_load, i+1);
					countCoursesPerPeriod(posible_solution, bs_total_courses, i+1);
					/*for (int l = 0; l < periods; ++l)
					{
						cout << bs_total_load[l];
					}*/
					
					/*
					cout << "\titeration " << load[j].period << "\n\tvalors: \n";
					cout << "\t" << countBrokenLoad(bs_total_load) << "\n";
					cout << "\t" << countBrokenCourses(bs_total_courses) << "\n";
					cout << "\t" << countBrokenPrecedences(posible_solution, pre, i+1) << "\n";
					*/

					load[j].total = countBrokenLoad(bs_total_load) + countBrokenCourses(bs_total_courses) + countBrokenPrecedences(posible_solution, pre, i+1);
					//load[j].total = evaluateSolution(posible_solution, i+1);
		}
		

		//cout << "load :\n";
		/*for (int i = 0; i < periods; ++i)
		{
			/* code 
		} */
		
		//cout << "ITERATION " << i << "\n";
		
		if (min_period >=0)
		{
			for (int j = 0; j <= min_period; ++j)
			{
				load[j].flag = 1;
			}
			sort(load, load + periods , sortLoadByFlag);
			
			sort(load, load + (periods - (min_period + 1)), sortLoadByTotal);
			
			sort(load + (periods - (min_period + 1)) , load + periods, sortLoadByTotal);
			
		}
		else
		{
			sort(load, load + periods, sortLoadByTotal);
		}
		
		for (int j = 0; j < length; ++j)
		{
			//cout << "iteration " << j << "\n";
			rcl[j] = load[j].period;
		}
		/*for (int j = 0; j < periods; ++j)
		{
			cout << "load: " << load[j].period << " " << load[j].total << " " << load[j].flag << "\n" ;	
		}
		for (int j = 0; j < length; ++j)
		{
			cout << "rcl: " << rcl[j] << "\n" ;
		}*/
		chosen = rand()%length;
		//cout << "rcl->chosen: "<< rcl[chosen] << "\n";
		
		/*for (int j = 0; j < periods; ++j)
		{
			if(load[j].period!=rcl[chosen]){
				posible_solution[i].period = load[j].period;
				setToZero(bs_total_load, periods);
				setToZero(bs_total_courses, periods);
				countLoadPerPeriod(posible_solution, bs_total_load, i+1);
				countCoursesPerPeriod(posible_solution, bs_total_courses, i+1);
				load[j].total -= countBrokenLoad(bs_total_load, i+1) + countBrokenCourses(bs_total_courses, i+1) + countBrokenPrecedences(posible_solution, pre, i+1);
				//load[j].total -= courses[i].credits;
			}
		}*/
		posible_solution[i].period = rcl[chosen];
		//showSolution(posible_solution, i+1);
	}
	
	/*cout << "solution: \n";
	for (int i = 0; i < num_courses; ++i)
  {
  	cout << posible_solution[i].period << ' ' << posible_solution[i].course.name << ' ' << posible_solution[i].course.credits <<"\n";
  }*/

  /*sort(load, load+periods, sortLoadByPeriod);
  for (int i = 0; i < periods; ++i)
  {
  	cout << "periodo: " << load[i].period << "\ttotal: " << load[i].total << '\n';
  }*/

	//cout << "return solution \n";
	//cout << "Built solution: \n";
	//showSolution(posible_solution);
	//showBrokens(posible_solution, pre);
	return posible_solution;
}




int main(int argc, char* argv[]){
	int random_seed = time(NULL);
	//int random_seed = time(NULL);
	cout << argc << '\n';
	if(argc == 7){
		random_seed = atoi(argv[6]);
	}
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

  preProcessing(courses, precedences);
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
  	cout << "-------------------------------------------------------------------------\n"; 
  	hillClimbing(best_solution, initial_solution, precedences, max_iter);
  	}
  
	cout << "=================================\n";
  cout << "\nBest solution: \n";
  showSolution(best_solution);
  showBrokens(best_solution, precedences);
  cout << "evaluation: " << evaluateSolution(best_solution) << "\n";
  cout << "=================================\n";
  cout << random_seed << '\n';

  return 0;
}


