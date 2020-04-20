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

std::vector<int> seqByDivide(int min, int max, int thread_num) {
	std::vector<int> primes;

	for (int i = min; i <= max; i++) {
		bool isPrime = true;

		for (int j = 2; j <= floor(sqrt(i)); j++) { // bez wyznaczenia liczb pierw. wœród dzielników
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
std::vector<int> seqByEratostenes(int min, int max, int thread_num) {
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
std::vector<int> seqByEratostenesWithoutFactors(int min, int max, int thread_num) {
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
std::vector<int> parByEratostenes(int min, int max, int thread_num) {
	std::vector<int> primes;
	std::vector<std::vector<bool>> isComposite(thread_num, std::vector<bool>(max - min + 1, false)); //tablice do "wykreœlania"
	
	int maxFactor = floor(std::sqrt(max));
	std::vector<int> primeFactors = getPrimeFactorsFromEratostenes(maxFactor); // wyznaczenie pierwiastków
	
	#pragma omp parallel num_threads(thread_num)
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
		
		for (int j = 1; j < thread_num; j++) {
			isComposite[0][i - min] = isComposite[0][i - min] | isComposite[j][i - min];
		}

		if (!isComposite[0][i - min]) {
			primes.push_back(i);
		}
	}
	
	return primes;
}

//procesy szukaj¹ dzielników dla w³asnych zakresów
std::vector<int> parByEratostenesInRanges(int min, int max, int thread_num) {
	int range = (max - min + 1) / thread_num;
	
	std::vector<int> primes;
	std::vector<std::vector<int>> threadsPrimes(thread_num); //tablice do "wykreœlania"


	#pragma omp parallel num_threads(thread_num)
	{
		// ka¿dy proces poszukuje pierwiastków w wyznaczonym dla niego przedziale
		int thread_num = omp_get_thread_num();

		int localMin = min + range * thread_num;
		int localMax = localMin + range - 1;
		if (thread_num == omp_get_num_threads() - 1) {
			localMax = max;
		}

		int maxFactor = floor(std::sqrt(localMax));
		std::vector<int> primes = seqByEratostenes(localMin, localMax, NULL); // wyznaczenie liczb pierwszych
		threadsPrimes[thread_num] = primes;
	}

	// zebranie wyników do kupy
	for (int thread = 0; thread < thread_num; thread++) {
		for (int prime : threadsPrimes[thread]) {
			primes.push_back(prime);
		}
	}

	return primes;
}


//procesy szukaj¹ dzielników dla w³asnych zakresów
//std::vector<int> parByEratostenesInRangesViaFor(int min, int max,  int thread_num) {
//	std::vector<int> primes;
//	std::vector<bool> isComposite(max + 1, false); //tablice do "wykreœlania"
//	int maxFactor = floor(std::sqrt(max));
//	std::vector<int> primeFactors = getPrimeFactorsFromEratostenes(maxFactor); // wyznaczenie pierwiastków
//
//#pragma omp parallel num_threads(thread_num)
//	{
//		std::cout << omp_get_thread_num()<<std::endl;
//		for (int i = 0; i < primeFactors.size(); i++) {
//			if ((primeFactors[i] >= min) && (omp_get_thread_num() == 0)) {
//				primes.push_back(primeFactors[i]);
//			}
//
//			int primeMultiple;
//			bool start = true;
//			int fakePrimeFactors = 1;
//			#pragma omp for
//			for (int j = min; j <= max; j += fakePrimeFactors) {
//				if (start) {
//					if(j % primeFactors[i] != 0) {
//						continue;
//					}
//					start = false;
//					fakePrimeFactors = primeFactors[i];
//				}
//				isComposite[j] = true;
//			}
//		}
//
//		for (int i = min; i <= max; i++) {
//			if (!isComposite[i]) {
//				primes.push_back(i);
//			}
//		}
//	}
//	return primes;
//}

//procesy szukaj¹ dzielników dla w³asnych zakresów, wspolny zbior pierwiastkow
std::vector<int> parByEratostenesCommonPrimesSet(int min, int max, int thread_num) {
	using Range = std::pair<int, int>;

	std::vector<int> primes;
	std::vector<std::vector<int>> threadsPrimes(thread_num); //tablice do "wykreœlania"

	std::vector<int> commonPrimes = getPrimeFactorsFromEratostenes(max); // wyznaczenie pierwiastków

	std::vector<Range> subranges(pow(thread_num, 2));
	int range = (max - min + 1) / thread_num;

	//wyznaczanie podprzedzia³ów
	for (int thread = 0; thread < thread_num; thread++) {
		int localMin = min + range * thread;
		int localMax = localMin + range - 1;
		if (thread == thread_num - 1) {
			localMax = max;
		}

		int subRange = (localMax - localMin + 1) / thread_num;

		for (int i = 0; i < thread_num; i++) {
			int subRangeLocalMin = localMin + subRange * i;
			int subRangeLocalMax = subRangeLocalMin + subRange - 1;
			if (thread == thread_num - 1) {
				subRangeLocalMax = localMax;
			}
			std::cout << "min: " << subRangeLocalMin << " max: " << subRangeLocalMax << std::endl;
			subranges[thread * thread_num + i] = Range(subRangeLocalMin, subRangeLocalMax);
		}
	}


	#pragma omp parallel num_threads(thread_num)
	{	
		std::vector<int> primes;
		// ka¿dy proces poszukuje pierwiastków w wyznaczonym dla niego przedziale
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < subranges.size(); i++) {
			auto localPrimes = seqByEratostenes(subranges[i].first, subranges[i].second, NULL);
			primes.insert(primes.end(), localPrimes.begin(), localPrimes.end());
		}
		
		threadsPrimes[omp_get_thread_num()] = primes;
	}

	// zebranie wyników do kupy
	for (int thread = 0; thread < thread_num; thread++) {
		for (int prime : threadsPrimes[thread]) {
			primes.push_back(prime);
		}
	}
	return primes;
}

int main() {
	const int MAX = 200;

	auto func = { 
		//seqByDivide,
		//seqByEratostenes,
		//seqByEratostenesWithoutFactors,
		//parByEratostenes,
		//parByEratostenesInRanges,
		parByEratostenesCommonPrimesSet
	};

	const auto test = [&func](int start, int stop, int num) {
		std::vector<std::vector<int>> results(func.size());
		for (int i = 0; i < func.size(); i++) {
			auto pack = measureTime(*std::next(func.begin(), i), start, stop, num);
			std::cout 
				<< "<" << start << " , " << stop << "> " 
				<< "Function " << "[" << i << "]" << std::endl 
				<< "Threads " << "[" << num << "]" << std::endl
				<< "Time: " << pack.first << std::endl;
			
			printResult(pack.second);
			std::cout << std::endl << std::endl;
			results[i] = pack.second;
		}

		return results;
	};

	const std::vector<int> thread_num = { 4,2,1 };

	for (int num : thread_num) {
		test(2, MAX, num);
		test(MAX / 2, MAX, num);
		test(2, MAX / 2, num);
	}

	return 0;
}