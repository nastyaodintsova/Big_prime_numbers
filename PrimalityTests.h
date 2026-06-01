#pragma once

#include "LongInt.h"

NTResult Trial_Compositeness_Check(const long_int &f, const long_int &max);

NTResult Strong_Pseudoprime(long_int &&base, const long_int &mod);

NTResult Strong_Pseudoprime(const long_int &base, const long_int &mod);

NTResult Mersenne_Primes(const long_int &n);

NTResult Lucas_Prime_Test(const long_int &n, const std::vector<long_int> &prime_factors,
                          const long_int &b);

NTResult Pepin_Test(long_int &&n);
NTResult Pepin_Test(const long_int &n);
