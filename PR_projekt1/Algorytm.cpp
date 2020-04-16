#include <vector>
#include <omp.h>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <algorithm>

void printResult(std::vector<int> primes) {
	std::cout << "znaleziono " << primes.size() << " liczb pierwszych w zadanym przedziale" << std::endl;
	
	for (int i = 0; i < primes.size(); i++) {
		std::cout << primes[i] << " ";
		if (i > 0 && i % 10 == 0) {
			std::cout << std::endl;
		}
	}
}

//TODO test efektywno�ci czy wyznacza� liczby pierwsze jako dzielnika - kiedy op�acalne

std::vector<int> getPrimeFactorsFromEratostenes(int maxFactor) {
	std::vector<int> primeFactors;
	std::vector<bool> isComposite(maxFactor + 1, false);

	int sqrt = floor(std::sqrt(maxFactor));
	
	for (int i = 2; i <= sqrt; i++) { // bez wyznaczenia liczb pierw. w�r�d dzielnik�w
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

		for (int j = 2; j <= maxFactor; j++) { // bez wyznaczenia liczb pierw. w�r�d dzielnik�w
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
	std::vector<int> primeFactors = getPrimeFactorsFromEratostenes(maxFactor); //wyznaczenie pierwiastk�w

	int pFactorMultiple; //wielokrotno�� pFactor >= min, otrzymujemy liczbe zlozona, czyli warto�� startow� sita

	for (int pFactor: primeFactors) {
		// wyznaczenie
		if (pFactor >= min) {
			primes.push_back(pFactor);
		}

		if (min % pFactor) {
			pFactorMultiple = min - min % pFactor + pFactor; //dodajemy pFactor aby by�a liczba wi�ksza od min oraz wielokrotno�ci� pFactor
		} else {
			pFactorMultiple = min;
		}

		for (int j = pFactorMultiple; j <= max; j += pFactor) { //usuwamy wielokrotno�ci liczb pierwszych w zakresie <min; max>
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

	for (int i = 2; i <= maxFactor; i++) { // bez wyznaczenia liczb pierw. w�r�d dzielnik�w
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

//procesy dziel� si� potencjalnymi dzielnikami z podzbioru liczb pierwszych
std::vector<int> parByEratostenes(int min, int max) {
	std::vector<int> primes;
	std::vector<std::vector<bool>> isComposite(omp_get_num_threads(), std::vector<bool>(max - min + 1, false)); //tablice do "wykre�lania"
	
	int maxFactor = floor(std::sqrt(max));
	std::vector<int> primeFactors = getPrimeFactorsFromEratostenes(maxFactor); // wyznaczenie pierwiastk�w
	
	#pragma omp parallel
	{
		// ka�dy proces pracuje na w�asnej tablicy
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

	// zebranie wynik�w do kupy
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


//procesy szukaj� dzielnik�w dla w�asnych zakres�w
std::vector<int> parByEratostenesInRanges(int min, int max) {
	int range = (max - min) / omp_get_num_threads();
	
	std::vector<int> primes;
	std::vector<std::vector<int>> threadsPrimes(omp_get_num_threads()); //tablice do "wykre�lania"


	#pragma omp parallel
	{
		// ka�dy proces poszukuje pierwiastk�w w wyznaczonym dla niego przedziale
		int thread_num = omp_get_thread_num();

		int localMin = min + range * thread_num;
		int localMax = localMin + range;


		int maxFactor = floor(std::sqrt(localMax));
		std::vector<int> primes = seqByEratostenes(localMin, localMax); // wyznaczenie liczb pierwszych
		threadsPrimes[thread_num] = primes;
	}

	// zebranie wynik�w do kupy
	for (int thread = 0; thread < omp_get_num_threads(); thread++) {
		for (int prime : threadsPrimes[thread]) {
			primes.push_back(prime);
		}
	}

	return primes;
}

int main() {
	printResult(parByEratostenesInRanges(2, 200));
}