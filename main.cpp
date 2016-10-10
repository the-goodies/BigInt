#include "lib/utilities/error.h"
#include "lib/containers/Array.h"
#include "lib/utilities/utility.h"
#include "BigInt.h"
#include <ctime>

using std::endl;
using std::cout;
using std::cin;

/*
void karatsuba_mult_correctness_test()
{
	const int LIMIT = 200;
	Array<u64> num_a_storage;
	Array<u64> num_b_storage;
	for (int i = 0; i <= LIMIT; ++i)
	{
		num_a_storage.add(getRandomInt(0, 1000000000 - 1));
		num_b_storage.add(getRandomInt(0, 1000000000 - 1));
	}

	int TESTS = 20;
	for (int test = 0; test < TESTS; ++test)
	{
		shuffle(num_a_storage);
		shuffle(num_b_storage);

		BigInt num_a(num_a_storage);
		BigInt num_b(num_b_storage);

		if (test % 2)
		{
			num_a *= -1;
		}

		if (num_a.kar(num_b) != num_a * num_b)
		{
			cout << test << ") " << "KARATSUBA NUMBER: " << num_a.kar(num_b) << " CORRECT NUMBER: " << num_a * num_b << endl;
		}
		if (num_a.kar("0") != num_a * "0")
		{
			cout << test << ") " << "KARATSUBA NUMBER: " << num_a.kar("0") << " CORRECT NUMBER: " << num_a * "0" << endl;
		}

		// unary - test
		if (num_b * -1 != -num_b)
		{
			cout << "num_b * -1 = " << num_b * -1 << " -num_b = " << -num_b << endl;
		}
	}
}
*/


/* not DONE 
void karatsuba_mult_efficiency_test()
{
	const int LIMIT = 200;
	Array<u64> num_a_storage;
	Array<u64> num_b_storage;
	for (int i = 0; i <= LIMIT; ++i)
	{
		num_a_storage.add(getRandomInt(1, 1000000000 - 1));
		num_b_storage.add(getRandomInt(1, 1000000000 - 1));
	}

	clock_t start_time;
	float time_elapsed;
	float total_time_elapsed = 0.0f;

	int TESTS = 20;
	for (int test = 0; test < TESTS; ++test)
	{
		shuffle(num_a_storage);
		shuffle(num_b_storage);

		BigInt num_a(num_a_storage);
		BigInt num_b(num_b_storage);

		//
		start_time = clock();
		
		time_elapsed = float(clock() - start_time) / CLOCKS_PER_SEC;
		//

		total_time_elapsed += time_elapsed;
	}
}
*/

int main()
{
	/*
	const int T = 1000;
	int count = 0;
	for (int t = 0; t < T; ++t)
	{
		int r1 = getRandomInt(-40000, 41474);
		int r2 = getRandomInt(-40000, 41474);
		int r1_minus_r2_int = r1 - r2;
		int r2_minus_r1_int = r2 - r1;
		int r1_plus_r2_int = r1 + r2;
		int r1_times_r2 = r1 * r2;
		int r1_divide_r2 = r1 / r2;
		int r2_divide_r1 = r2 / r1;
		BigInt r1_bigint(r1);
		BigInt r2_bigint(r2);

		BigInt r1_minus_r2_bigint(r1_bigint - r2_bigint);
		BigInt r2_minus_r1_bigint(r2_bigint - r1_bigint);
		BigInt r1_plus_r2_bigint(r1_bigint + r2_bigint);
		BigInt r2_plus_r1_bigint(r2_bigint + r1_bigint);
		BigInt r1_times_r2_bigint(r1_bigint * r2_bigint);
		BigInt r2_times_r1_bigint(r2_bigint * r1_bigint);
		BigInt r1_divide_r2_bigint(r1_bigint / r2_bigint);
		BigInt r2_divide_r1_bigint(r2_bigint / r1_bigint);

		if (r1_times_r2 != r1_times_r2_bigint.convertToInt()) cout << "r1_times_r2_bigint: " << r1_times_r2_bigint << " r1_times_r2_int: " << r1_times_r2 << endl;
		if (r1_minus_r2_int != r1_minus_r2_bigint.convertToInt()) cout << "r1_minus_r2_bigint: " << r1_minus_r2_bigint << " r1_minus_r2_int: " << r1_minus_r2_int <<  endl;
		if (r2_minus_r1_int != r2_minus_r1_bigint.convertToInt()) cout << "r2_minus_r1_bigint: " << r2_minus_r1_bigint << " r2_minus_r1_int: " << r2_minus_r1_int <<  endl;
		if (r1_plus_r2_int != r1_plus_r2_bigint.convertToInt()) cout << "r1_plus_r2_bigint: " << r1_plus_r2_bigint << endl;
		if (r1_divide_r2 != r1_divide_r2_bigint.convertToInt()) cout << "T: " << t << "   " << r1_divide_r2_bigint.convertToInt() << "\n";
		if (r2_divide_r1 != r2_divide_r1_bigint.convertToInt()) cout << "T: " << t << "   " << r2_divide_r1_bigint.convertToInt() << "\n";
	}

	BigInt c("100010001000100045646548949816516548941000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498498946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548941000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498498946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816100010001000100045646548949816516548941000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498498946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548941000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498498946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849516548948946516848948115611000100010001000456465489498165165489410001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984989465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489410001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984989465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498161000100010001000456465489498165165489410001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984989465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489410001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984910001000100010004564654894981651654894894651684894811561654984989465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498495165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498491000100010001000456465489498165165489489465168489481156165498496549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849");
	BigInt d("-100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849100010001000100045646548949816516548948946516848948115616549849");
	BigInt a("1234500000000000424236363424");
	BigInt b("-12345678901234243634636436365");

	clock_t start_time;
	float time_elapsed;

	start_time = clock();
	c * c;
	time_elapsed = float(clock() - start_time) / CLOCKS_PER_SEC;
	cout << "TIME ELAPSED: " << time_elapsed << endl;

	start_time = clock();
	c.kar(c);
	time_elapsed = float(clock() - start_time) / CLOCKS_PER_SEC;
	cout << "TIME ELAPSED: " << time_elapsed << endl;

	system("pause");
	*/
	return 0;
}
