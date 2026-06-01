#include <boost/multiprecision/gmp.hpp>
#include <iostream>

#pragma once

using long_int = boost::multiprecision::mpz_int;

enum class NTStatus {
    Composite,
    ProbablyPrime,
    Prime,

    FactorFound,
    CompleteFactorization,

    PrimeGenerated,

    Failure,
    InvalidInput
};

struct NTResult {
    NTStatus status;

    long_int divisor = 0;
    std::vector<long_int> primes;
    std::vector<long_int> degs;
};

struct FactorPower {
    long_int cofactor;
    long_int factor;
    long_int multiplicity;
};
