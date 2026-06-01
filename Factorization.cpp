#include "Factorization.h"
#include "Arithmetic.h"
#include "Sieves.h"

NTResult Trial_Division(long_int &&f, const uint64_t &max) {
    NTResult fp;
    long long i = 0;

    for (int k = 2; k < 4; ++k) {
        if (f % k == 0) {
            fp.primes.push_back(k);
            fp.degs.push_back(1);
            f /= k;

            auto r = Extract_Factor_Power(f, k);
            f = r.cofactor;
            fp.degs[i] = r.multiplicity;
            ++i;
        }
    }

    long_int d = 5;
    int add = 2;

    while (d <= max && d * d <= f) {
        if (f % d == 0) {
            fp.primes.push_back(d);
            fp.degs.push_back(1);
            f /= d;

            auto r = Extract_Factor_Power(f, d);
            f = r.cofactor;
            fp.degs[i] = r.multiplicity;
            ++i;
        }

        d += add;
        add = 6 - add;
    }

    if (d * d > f) {
        fp.primes.push_back(f);
        fp.degs.push_back(1);
        ++i;
    }

    if (i == 1) {
        fp.status = NTStatus::Prime;
        return fp;
    } else {
        fp.status = NTStatus::CompleteFactorization;
        return fp;
    }
}

NTResult Trial_Division(const long_int &f, const long_int &max) {
    long_int safe_f = f;
    return Trial_Division(std::move(safe_f), max);
}

NTResult Pollard_Rho(const long_int &n, const long_int &c, const long_int &max) {

    long_int x1 = 2;
    long_int x2 = c + 4;
    long_int range = 1;
    long_int terms = 0;
    long_int product = 1;

    while (terms <= max) {
        for (long_int j = 1; j < range && terms < max; ++j) {
            x2 = (x2 * x2 + c) % n;
            product *= (x1 - x2) % n;
            ++terms;

            if (terms % 10 == 0) {
                long_int g = Gcd(n, product);

                if (g > 1 && g < n) {
                    return {NTStatus::FactorFound, g};
                }

                if (g == n) {
                    return {NTStatus::Failure};
                }

                product = 1;
            }
        }

        x1 = std::move(x2);
        range *= 2;
    }

    return {NTStatus::Failure};
}

NTResult Pollard_P_Optimised(const long_int &n, const long_int &c, int64_t limit) {

    long_int m = c % n;
    auto res = segmented_sieve(limit);
    std::vector<int64_t> primes;
    primes.reserve(res.primes.size());

    for (auto &x : res.primes) {
        primes.push_back(x.convert_to<int64_t>());
    }

    for (int64_t p : primes) {

        if (p < 2) {
            continue;
        }

        int64_t e = 0;
        long long pp = p;

        while (pp <= limit) {
            pp *= p;
            ++e;
        }

        long_int pe = Fast_Pow(p, e);
        m = Fast_Pow(m, pe, n);
    }

    long_int g = Gcd(m - 1, n);

    if (g > 1 && g < n) {
        return {NTStatus::FactorFound, g};
    }

    return {NTStatus::Failure};
}

NTResult CFRAC(long_int &&n, int limit = 10000) {
    long_int sqrt = boost::multiprecision::sqrt(n);
    long_int b = 1;
    long_int prev_a = 1;
    long_int prev_b = 0;

    std::vector<long_int> x_list;
    std::vector<long_int> y_list;

    for (int i = 0; i < limit; ++i) {
        long_int ai = (sqrt + prev_a) / b;
        long_int next_a = ai * b - prev_a;
        long_int next_b = (n - next_a * next_a) / b;

        if (next_b == 0)
            break;

        long_int xi = next_a;
        long_int yi = (xi * xi) % n;

        long_int g = Gcd(std::move(yi), std::move(n));

        if (g > 1 && g < n) {
            return {NTStatus::FactorFound, g};
        }

        prev_a = std::move(next_a);
        b = std::move(next_b);
    }

    return {NTStatus::Failure};
}

NTResult CFRAC(const long_int &n, int limit = 10000) {
    long_int safe_n = n;
    return CFRAC(std::move(safe_n), limit);
}
