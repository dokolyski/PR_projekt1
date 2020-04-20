#include <vector>
#include <omp.h>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <utility>
#include <chrono>

using Time = double;
template<typename Func, typename ... Vargs>
auto measureTime(Func f, Vargs ... args) -> std::pair<Time, std::vector<int>> {
	using namespace std;

	auto startTimePoint = chrono::high_resolution_clock::now();
	auto results = f(args...);
	auto endTimePoint = chrono::high_resolution_clock::now();

	auto time = chrono::duration_cast<chrono::duration<double>>(endTimePoint - startTimePoint);
	return std::pair< Time, decltype(results)>(time.count(), results);
}

void printResult(std::vector<int> primes) {
	std::cout << "znaleziono " << primes.size() << " liczb pierwszych w zadanym przedziale" << std::endl;
	
	for (int i = 0; i < primes.size(); i++) {
		std::cout << primes[i] << " ";
		if (i > 0 && i % 10 == 0) {
			std::cout << std::endl;
		}
	}
}

//TODO test efektywnoœci czy wyznaczaæ liczby pierwsze jako dzielnika - kiedy op³acalne

std::vector<int> getPrimeFactorsFromEratostenes(int maxFactor) {
	std::vector<int> primeFactors;
	std::vector<bool> isComposite(maxFactor + 1, false);

	int sqrt = floor(std::sqrt(maxFactor));
	
	for (int i = 2; i <= sqrt; i++) { // bez wyznaczenia liczb pierw. wœród dzielników
		if (isComposite[i]) {
			continue;
		} else {
			primeFactors.push_back(i);
			for (int j = i; j <= maxFactor; j += i) {
				isComposite[j] = true;
			}
		}
	}

	for (int i = sqrt + 1; i <= maxFactor; i++) {
		if (!isComposite[i]) {
			primeFactors.push_back(i);
		}
	}
	
	return primeFactors;
}

std::vector<int> seqByDivide(int min, int max) {
	std::vector<int> primes;
	int maxFactor = floor(sqrt(max));

	for (int i = min; i <= max; i++) {
		bool isPrime = true;

		for (int j = 2; j <= maxFactor; j++) { // bez wyznaczenia liczb pierw. wœród dzielników
			if (i % j == 0) {
				isPrime = false;
				break;
			}
		}

		if (isPrime) {
			primes.push_back(i);
		}
	}
	return primes;
}

//sito w przedziale, dzielniki pierwsze w przedziale <2; sqrt(max)>
std::vector<int> seqByEratostenes(int min, int max) {
	std::vector<int> primes;
	std::vector<bool> isComposite(max - min + 1, false);

	int maxFactor = floor(std::sqrt(max));
	std::vector<int> primeFactors = getPrimeFactorsFromEratostenes(maxFactor); //wyznaczenie pierwiastków

	int pFactorMultiple; //wielokrotnoœæ pFactor >= min, otrzymujemy liczbe zlozona, czyli wartoœæ startow¹ sita

	for (int pFactor: primeFactors) {
		// wyznaczenie
		if (pFactor >= min) {
			primes.push_back(pFactor);
		}

		if (min % pFactor) {
			pFactorMultiple = min - min % pFactor + pFactor; //dodajemy pFactor aby by³a liczba wiêksza od min oraz wielokrotnoœci¹ pFactor
		} else {
			pFactorMultiple = min;
		}

		for (int j = pFactorMultiple; j <= max; j += pFactor) { //usuwamy wielokrotnoœci liczb pierwszych w zakresie <min; max>
			isComposite[j - min] = true;
		}
	}

	for (int i = min; i <= max; i++) {
		if (!isComposite[i - min]) {
			primes.push_back(i);
		}
	}
	
	return primes;
}

//sito w przedziale, wszystkie dzielniki w przedziale <2; sqrt(max)>
std::vector<int> seqByEratostenesWithoutFactors(int min, int max) {
	std::vector<int> primes;
	std::vector<bool> isComposite(max + 1, false);

	int maxFactor = floor(std::sqrt(max));
	int primeMultiple;

	for (int i = 2; i <= maxFactor; i++) { // bez wyznaczenia liczb pierw. wœród dzielników
		if (i >= min && !isComposite[i]) {
			primes.push_back(i);
		}

		if (min % i) {
			primeMultiple = min - min % i + i;
		}
		else {
			primeMultiple = min;
		}

		for (int j = primeMultiple; j <= max; j += i) {
			isComposite[j] = true;
		}
	}

	for (int i = min; i <= max; i++) {
		if (!isComposite[i]) {
			primes.push_back(i);
		}
	}
	
	return primes;
}

//procesy dziel¹ siê potencjalnymi dzielnikami z podzbioru liczb pierwszych
std::vector<int> parByEratostenes(int min, int max) {
	std::vector<int> primes;
	std::vector<std::vector<bool>> isComposite(omp_get_num_threads(), std::vector<bool>(max - min + 1, false)); //tablice do "wykreœlania"
	
	int maxFactor = floor(std::sqrt(max));
	std::vector<int> primeFactors = getPrimeFactorsFromEratostenes(maxFactor); // wyznaczenie pierwiastków
	
	#pragma omp parallel
	{
		// ka¿dy proces pracuje na w³asnej tablicy
		int thread_num = omp_get_thread_num();
		int primeMultiple;


		#pragma omp for
		for (int i = 0; i < primeFactors.size(); i++) {
			if (primeFactors[i] >= min) {
				primes.push_back(primeFactors[i]);
			}

			if (min % primeFactors[i]) {
				primeMultiple = min - min % primeFactors[i] + primeFactors[i];
			} 
			else {
				primeMultiple = min;
			}

			for (int j = primeMultiple; j <= max; j += primeFactors[i]) {
				isComposite[thread_num][j - min] = true;
			}
		}
	}

	// zebranie wyników do kupy
	for (int i = min; i <= max; i++) {
		
		for (int j = 1; j < omp_get_num_threads(); j++) {
			isComposite[0][i - min] = isComposite[0][i - min] | isComposite[j][i - min];
		}

		if (!isComposite[0][i - min]) {
			primes.push_back(i);
		}
	}
	
	return primes;
}

//procesy szukaj¹ dzielników dla w³asnych zakresów
std::vector<int> parByEratostenesInRanges(int min, int max) {
	int range = (max - min + 1) / omp_get_num_threads();
	
	std::vector<int> primes;
	std::vector<std::vector<int>> threadsPrimes(omp_get_num_threads()); //tablice do "wykreœlania"


	#pragma omp parallel
	{
		// ka¿dy proces poszukuje pierwiastków w wyznaczonym dla niego przedziale
		int thread_num = omp_get_thread_num();

		int localMin = min + range * thread_num;
		int localMax = localMin + range - 1;


		int maxFactor = floor(std::sqrt(localMax));
		std::vector<int> primes = seqByEratostenes(localMin, localMax); // wyznaczenie liczb pierwszych
		threadsPrimes[thread_num] = primes;
	}

	// zebranie wyników do kupy
	for (int thread = 0; thread < omp_get_num_threads(); thread++) {
		for (int prime : threadsPrimes[thread]) {
			primes.push_back(prime);
		}
	}

	return primes;
}


//procesy szukaj¹ dzielników dla w³asnych zakresów
std::vector<int> parByEratostenesInRangesViaFor(int min, int max) {
	std::vector<int> primes;
	std::vector<bool> isComposite(max + 1, false); //tablice do "wykreœlania"
	int maxFactor = floor(std::sqrt(max));
	std::vector<int> primeFactors = getPrimeFactorsFromEratostenes(maxFactor); // wyznaczenie pierwiastków

#pragma omp parallel num_threads(4)
	{
		std::cout << omp_get_thread_num()<<std::endl;
		for (int i = 0; i < primeFactors.size(); i++) {
			if ((primeFactors[i] >= min) && (omp_get_thread_num() == 0)) {
				primes.push_back(primeFactors[i]);
			}

			int primeMultiple;
			bool start = true;
			int fakePrimeFactors = 1;
			#pragma omp for
			for (int j = min; j <= max; j += fakePrimeFactors) {
				if (start) {
					if(j % primeFactors[i] != 0) {
						continue;
					}
					start = false;
					fakePrimeFactors = primeFactors[i];
				}
				isComposite[j] = true;
			}
		}

		for (int i = min; i <= max; i++) {
			if (!isComposite[i]) {
				primes.push_back(i);
			}
		}
	}
	return primes;
}

//procesy szukaj¹ dzielników dla w³asnych zakresów, wspolny zbior pierwiastkow
std::vector<int> parByEratostenesCommonPrimesSet(int min, int max) {
	using Range = std::pair<int, int>;

	std::vector<int> primes;
	std::vector<std::vector<int>> threadsPrimes(omp_get_num_threads()); //tablice do "wykreœlania"

	std::vector<int> commonPrimes = getPrimeFactorsFromEratostenes(max); // wyznaczenie pierwiastków

	std::vector<Range> subranges(pow(omp_get_num_threads(), 2));
	int range = (max - min + 1) / omp_get_num_threads();

	//wyznaczanie podprzedzia³ów
	for (int thread = 0; thread < omp_get_num_threads(); thread++) {
		int localMin = min + range * thread;
		int localMax = localMin + range - 1;

		int subRange = (localMax - localMin + 1) / omp_get_num_threads();

		for (int i = 0; i < omp_get_num_threads(); i++) {
			int subRangeLocalMin = localMin + subRange * i;
			int subRangeLocalMax = localMin + subRange - 1;
			
			subranges[thread * omp_get_num_threads() + i] = Range(subRangeLocalMin, subRangeLocalMax);
		}
	}


	#pragma omp parallel
	{	
		std::vector<int> primes;

		// ka¿dy proces poszukuje pierwiastków w wyznaczonym dla niego przedziale
		#pragma omp for dynamic
		for (int i = 0; i < subranges.size(); i++) {
			auto localPrimes = seqByEratostenes(subranges[i].first, subranges[i].second);
			primes.insert(primes.end(), localPrimes.begin(), localPrimes.end());
		}
		
		threadsPrimes[omp_get_thread_num()] = primes;
	}

	// zebranie wyników do kupy
	for (int thread = 0; thread < omp_get_num_threads(); thread++) {
		for (int prime : threadsPrimes[thread]) {
			primes.push_back(prime);
		}
	}
	return primes;
}

int main() {
	const int MAX = 200;

	auto func = { 
		seqByDivide,
		seqByEratostenes,
		seqByEratostenesWithoutFactors,
		parByEratostenes,
		parByEratostenesInRanges,
		parByEratostenesCommonPrimesSet
	};

	const auto test = [&func](int start, int stop) {
		std::vector<std::vector<int>> results(func.size);
		for (int i = 0; i < func.size(); i++) {
			auto pack = measureTime(*std::next(func.begin(), i), start, stop);
			std::cout 
				<< "<" << start << " , " << stop << "> " 
				<< "Function " << "[" << i << "]" << std::endl 
				<< "Time: " << pack.first << std::endl;
			
			printResult(pack.second);
			std::cout << std::endl << std::endl;
			results[i] = pack.second;
		}

		return results;
	};

	/*
	const int SEQUENTIAL = 1;
	const int PARALLEL_MAX_LOGIC_THREADS = 4;
	const int PARALLEL_MAX_PHYSICAL_THREADS = 4;
	const int PARALLEL_HALF_PHYSICAL_THREADS = 2;

	omp_set_num_threads(SEQUENTIAL);
	if (PARALLEL_MAX_PHYSICAL_THREADS != PARALLEL_MAX_LOGIC_THREADS) {
		omp_set_num_threads(PARALLEL_MAX_PHYSICAL_THREADS);
	}
	
	if (PARALLEL_HALF_PHYSICAL_THREADS != SEQUENTIAL &&
		PARALLEL_HALF_PHYSICAL_THREADS != PARALLEL_MAX_LOGIC_THREADS &&
		PARALLEL_HALF_PHYSICAL_THREADS != PARALLEL_MAX_PHYSICAL_THREADS) {
		omp_set_num_threads(PARALLEL_HALF_PHYSICAL_THREADS);
	}
	*/

	test(2, MAX);
	test(MAX / 2, MAX);
	test(2, MAX / 2);
}