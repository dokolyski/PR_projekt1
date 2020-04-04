#include <vector>
#include <omp.h>
#include <cstdio>
#include <cmath>

void printResult(std::vector<int> result) {
	printf("znaleziono %d liczb pierwszych w zadanym przedziale", result.size());
	for (int i = 0; i < result.size(); i++) {
		if (i % 10 == 0) {
			printf("\n");
		}
		printf("%d ", result[i]);
	}
}

//TODO test efektywnoœci czy wyznaczaæ liczby pierwsze jako dzielnika - kiedy op³acalne

std::vector<int> getFactorsFromEratostenes(int maxFactor) {
	std::vector<int> factors;
	bool* isComposite = new bool[maxFactor + 1]();
	int sqrt = floor(std::sqrt(maxFactor));
	for (int i = 2; i <= sqrt; i++) { // bez wyznaczenia liczb pierw. wœród dzielników
		if (isComposite[i]) continue;
		factors.push_back(i);
		for (int j = i; j <= maxFactor; j += i) {
			isComposite[j] = true;
		}
	}
	for (int i = sqrt + 1; i <= maxFactor; i++) {
		if (!isComposite[i]) {
			factors.push_back(i);
		}
	}
	delete[] isComposite;
	return factors;
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

std::vector<int> seqByEratostenes(int min, int max) {
	std::vector<int> primes;
	bool* isComposite = new bool[max - min + 1]();
	int maxFactor = floor(std::sqrt(max));
	std::vector<int> factors = getFactorsFromEratostenes(maxFactor); //wyznaczenie pierwiastków
	for (int i = 0; i < factors.size(); i++) {
		// wyznaczenie
		int j;
		if (min % factors[i]) {
			j = min + factors[i] - min % factors[i];
		}
		else {
			j = min;
		}

		for (j; j <= max; j += factors[i]) {
			isComposite[j - min] = true;
		}
	}
	for (int i = 0; i < factors.size(); i++) {
		if (factors[i] >= min) {
			primes.push_back(factors[i]);
		}
	}
	for (int i = min; i <= max; i++) {
		if (!isComposite[i - min]) {
			primes.push_back(i);
		}
	}
	delete[] isComposite;
	return primes;
}

std::vector<int> seqByEratostenesWithoutFactors(int min, int max) {
	std::vector<int> primes;
	bool* isComposite = new bool[max+1]();
	int maxFactor = floor(std::sqrt(max));
	for (int i = 2; i <= maxFactor; i++) { // bez wyznaczenia liczb pierw. wœród dzielników
		if ((i >= min)&&(!isComposite[i])) {
			primes.push_back(i);
		}
		int j;
		if (min % i) {
			j = min + i - (min % i);
		}
		else {
			j = min;
		}
		for (j; j <= max; j += i) {
			isComposite[j] = true;
		}
	}
	for (int i = min; i <= max; i++) {
		if (!isComposite[i]) {
			primes.push_back(i);
		}
	}
	delete[] isComposite;
	return primes;
}

std::vector<int> parByEratostenes(int min, int max) { // procesy dziel¹ siê potencjalnymi dzielnikami
	std::vector<int> primes;
	bool** isCompositeArray = new bool* [omp_get_num_threads()]; // tablice do "wykreœlania"
	int maxFactor = floor(std::sqrt(max));
	std::vector<int> factors = getFactorsFromEratostenes(maxFactor); // wyznaczenie pierwiastków
	#pragma omp parallel
	{
		// ka¿dy proces pracuje na w³asnej tablicy
		int thread_num = omp_get_thread_num();
		isCompositeArray[thread_num] = new bool[max - min + 1]();
		#pragma omp for
		for (int i = 0; i < factors.size(); i++) {
			int j;
			if (min % factors[i]) {
				j = min + factors[i] - (min % factors[i]);
			}
			else {
				j = min;
			}
			for (j; j <= max; j += factors[i]) {
				isCompositeArray[thread_num][j - min] = true;
			}
		}
	}
	for (int i = 0; i < factors.size(); i++) {
		if (factors[i] >= min) {
			primes.push_back(factors[i]);
		}
	}
	// zebranie wyników do kupy
	for (int i = min; i <= max; i++) {
		for (int j = 1; j < omp_get_num_threads(); j++) {
			isCompositeArray[0][i - min] |= isCompositeArray[j][i - min];
		}
		if (!isCompositeArray[0][i - min]) {
			primes.push_back(i);
		}
	}
	// zwolnienie pamiêci
	for (int j = 0; j < omp_get_num_threads(); j++) {
		delete[] isCompositeArray[j];
	}
	delete[] isCompositeArray;
	return primes;
}



int main() {
	printResult(seqByEratostenes(2, 200));
	//printResult(parByEratostenes(2, 200));
}