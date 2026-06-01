#include "PrimalityTests.h"
#include "Arithmetic.h"

NTResult Trial_Compositeness_Check(const long_int &f, const long_int &max) {
    if (f <= 1) {
        return {NTStatus::Composite};
    }

    if (f % 2 == 0) {
        if (f == 2) {
            return {NTStatus::Prime};
        }
        return {NTStatus::FactorFound, 2};
    }

    long_int d = 3;

    while (d * d <= f && d <= max) {
        if (f % d == 0) {
            return {NTStatus::FactorFound, d};
        }

        d += 2;
    }

    if (d * d > f) {
        return {NTStatus::ProbablyPrime};
    }

    return {NTStatus::ProbablyPrime};
}

NTResult Strong_Pseudoprime(long_int &&base, const long_int &mod) {
    auto fp = Extract_Factor_Power(mod - 1, 2);

    long_int t = fp.cofactor;
    long_int a = fp.multiplicity;

    long_int test = Fast_Pow(std::move(base), std::move(t), mod);

    if (test == 1 || test == mod - 1) {
        return {NTStatus::ProbablyPrime};
    }

    for (long_int i = 0; i < a; ++i) {
        test *= test;
        test %= mod;

        if (test == mod - 1) {
            return {NTStatus::ProbablyPrime};
        }
    }

    return {NTStatus::Composite};
}

NTResult Strong_Pseudoprime(const long_int &safe_base, const long_int &mod) {
    long_int base = safe_base;
    return Strong_Pseudoprime(std::move(base), mod);
}

NTResult Mersenne_Primes(const long_int &n) {
    long_int m = Fast_Pow(2, n);
    long_int s = 4;

    for (long_int i = 2; i <= n - 1; ++i) {
        s = (s * s - 2) % m;

        if (s == 0) {
            return {NTStatus::Prime};
        }
    }

    return {NTStatus::Composite};
}

NTResult Lucas_Prime_Test(const long_int &n, const std::vector<long_int> &prime_factors,
                          const long_int &b) {

    long_int k = n - 1;
    if (Fast_Pow(b, k, n) != 1) {
        return {NTStatus::Composite};
    }
    for (const auto &p : prime_factors) {
        if (Fast_Pow(b, k / p, n) == 1) {
            return {NTStatus::Composite};
        }
    }

    return {NTStatus::Prime};
}

NTResult Pepin_Test(long_int &&n) {
    long_int num = Fast_Pow(2, std::move(n)) + 1;

    long_int exp = (num - 1) / 2;
    long_int res = Fast_Pow(3, std::move(exp), num);

    if (res == num - 1) {
        return {NTStatus::Prime};
    }

    return {NTStatus::Composite};
}

NTResult Pepin_Test(const long_int &n) {
    long_int safe_n = n;
    return Pepin_Test(std::move(safe_n));
}
