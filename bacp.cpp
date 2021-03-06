#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>     	
#include <ctime>
#include <cmath>
#include <algorithm>
using namespace std;

int length;
int years, periods, num_periods, min_courses, max_courses, min_load, max_load, num_courses, num_precedences, p;
float average = 0;

//Estructura para guardar un curso
struct Course{
	string name;
	int credits;
};

//Estructuras para guardar los prerrequisitos
struct Precedence{
	Course first;
	Course second;
};

//Estructura que representa un nodo de la solución
struct Solution{
	int period;
	Course course;
};

//Estructuras auxiliares que ayudan a los cálculos y al acceso de datos
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

//Funciones de ordenamiento
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

//Función de preprocesamiento que ordena los cursos dependiendo un ratio de "ser prerrequsito" 
//y "tener prerrequisitos"
void preProcessing(Course* co, Precedence* pre){
	PreEvaluation* eval = new PreEvaluation[num_courses];
	for (int i = 0; i < num_courses; ++i)
	{
		eval[i].course = co[i];
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
	}
}

//Funciones auxiliares que ayudan a contar la carga y la cantidad de cursos por periodo
void countLoadPerPeriod(Solution* solution_to_count, int* load_per_period, int stop = num_courses){
	for (int i = 0; i < stop; ++i)
	{
		load_per_period[solution_to_count[i].period]+=solution_to_count[i].course.credits;
	}
}
void countCoursesPerPeriod(Solution* solution_to_count, int* courses_per_period, int stop = num_courses){
	for (int i = 0; i < stop; ++i)
	{
		courses_per_period[solution_to_count[i].period]+=1;
	}
}

//Función para contar cuántos prerrequisitos rotos tiene la solución sol
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

//Función para contar cuántas restricciones de carga por periodo se rompen
float countBrokenLoad(int* load_per_period, int stop = periods){
	float broken = 0;
	for (int i = 0; i < stop; ++i)
	{
		if (load_per_period[i] < min_load || load_per_period[i] > max_load)
		{
			broken += pow(abs((load_per_period[i] - average)), 1);
		}
	}
	return broken;
}

//Función para contar cuántas restricciones de cantidad de cursos por periodo se rompen
float countBrokenCourses(int* courses_per_period, int stop = periods){
	float broken = 0;
	float average_courses = num_courses/periods;
	for (int i = 0; i < stop; ++i)
	{
		if (courses_per_period[i] < min_courses || courses_per_period[i] > max_courses)
		{
			broken += pow(abs((courses_per_period[i] - average_courses)), 1);
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

//Función de Evaluación, parámetro p se ingresa al ejecutar el programa
float evaluateSolution(Solution* solution_to_evaluate, int stop = num_courses){
	float evaluation = 0;
	int total_load_aux[periods]={0};
	countLoadPerPeriod(solution_to_evaluate, total_load_aux, stop);
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

//Post Procesamiento: HillClimbing + mejor mejora
void hillClimbing(Solution* best, Solution* initial, Precedence* pre, int iter){
	int hc_total_load[periods]={0};
	int hc_total_courses[periods]={0};
	int elem;
	int best_candidate;
	int aux;
	int counter;
	//crear alternative_solution
	Solution* alternative_solution = new Solution[num_courses];
	copy(initial, initial+num_courses,alternative_solution);
	Candidate* candidate_periods = new Candidate[periods];
	Candidate initial_candidate;
	for (int j = 0; j < iter; ++j)
	{
		//elegir elemento de alternative_solution
		elem = rand()%num_courses;
		best_candidate = initial[elem].period;
		//modificar alternative_solution
		for (int i = 0; i < periods; ++i)
		{
			alternative_solution[elem].period = i;
			setToZero(hc_total_load, periods);
			setToZero(hc_total_courses, periods);
			candidate_periods[i].period = i;
			countLoadPerPeriod(alternative_solution, hc_total_load);
			countCoursesPerPeriod(alternative_solution, hc_total_courses);
			candidate_periods[i].brokens[0] = countBrokenLoad(hc_total_load);
			candidate_periods[i].brokens[1] = countBrokenCourses(hc_total_courses);
			candidate_periods[i].brokens[2] = countBrokenPrecedences(alternative_solution, pre);
			candidate_periods[i].broken_constraints = candidate_periods[i].brokens[0] + candidate_periods[i].brokens[1] + candidate_periods[i].brokens[2];
			candidate_periods[i].eval = evaluateSolution(alternative_solution);
			if (i==best_candidate)
			{
				initial_candidate = candidate_periods[i];
			}
		}
		sort(candidate_periods, candidate_periods+periods, sortCandidateByBrokenConstraints);
		counter = 0;
		for (int i = 0; i < periods; ++i)
		{
			if(candidate_periods[i].broken_constraints == candidate_periods[0].broken_constraints)
			{
				counter++;
			}
		}
		sort(candidate_periods, candidate_periods + counter, sortCandidateByEvaluation);
		if (candidate_periods[0].broken_constraints < initial_candidate.broken_constraints)
		{
			alternative_solution[elem].period = candidate_periods[0].period;
			copy(alternative_solution, alternative_solution + num_courses, initial);
		}
		else if (candidate_periods[0].broken_constraints == initial_candidate.broken_constraints)
		{
			aux = alternative_solution[elem].period;
			alternative_solution[elem].period = candidate_periods[0].period;
			if (evaluateSolution(alternative_solution) <= evaluateSolution(initial))
			{
				copy(alternative_solution, alternative_solution + num_courses, initial);
			}
			else
			{
				alternative_solution[elem].period = aux;
			}
		}
	}
  setToZero(hc_total_load, periods);
  setToZero(hc_total_courses, periods);
	countLoadPerPeriod(best, hc_total_load);
	countCoursesPerPeriod(best, hc_total_courses);
	float best_solution_broken_constraints = countBrokenLoad(hc_total_load) + countBrokenCourses(hc_total_courses) + countBrokenPrecedences(best, pre);

  setToZero(hc_total_load, periods);
  setToZero(hc_total_courses, periods);
	countLoadPerPeriod(initial, hc_total_load);
	countCoursesPerPeriod(initial, hc_total_courses);
	float initial_solution_broken_constraints = countBrokenLoad(hc_total_load) + countBrokenCourses(hc_total_courses) + countBrokenPrecedences(initial, pre);

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
}

//Función para construir una solución. Greedy + RCL
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
					countLoadPerPeriod(posible_solution, bs_total_load, i+1);
					countCoursesPerPeriod(posible_solution, bs_total_courses, i+1);
					load[j].total = countBrokenLoad(bs_total_load) + countBrokenCourses(bs_total_courses) + countBrokenPrecedences(posible_solution, pre, i+1);
		}
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
			rcl[j] = load[j].period;
		}
		chosen = rand()%length;
		posible_solution[i].period = rcl[chosen];
	}
	return posible_solution;
}

int main(int argc, char* argv[]){
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
  preProcessing(courses, precedences);
  periods = years*num_periods;
	length = (atoi(argv[5])<periods)?atoi(argv[5]):periods;
  average = calculateAverage(courses);
  Solution* best_solution = buildSolution(courses, precedences);
  Solution* initial_solution;

 //Se ejecuta <restart> veces un HC, incluyendo crear su solución inicial
  for (int i = 0; i < restart; ++i)
  {
  	initial_solution = buildSolution(courses, precedences);
  	hillClimbing(best_solution, initial_solution, precedences, max_iter);
  	}
  return 0;
}