#include "Factorization.h"
#include "LongInt.h"
#include "PrimalityTests.h"
#include "Sieves.h"

#include <future>
#include <iostream>
#include <mutex>

namespace {

void Print_Status(const NTResult &res) {
    switch (res.status) {
    case NTStatus::Composite:
        std::cout << "Composite\n";
        break;

    case NTStatus::ProbablyPrime:
        std::cout << "ProbablyPrime\n";
        break;

    case NTStatus::Prime:
        std::cout << "Prime\n";
        break;

    case NTStatus::FactorFound:
        std::cout << "FactorFound\n";
        if (res.divisor != 0) {
            std::cout << res.divisor << '\n';
        }
        break;

    case NTStatus::CompleteFactorization:
        std::cout << "CompleteFactorization\n";
        for (size_t i = 0; i < res.primes.size(); ++i) {
            std::cout << res.primes[i] << '^' << (i < res.degs.size() ? res.degs[i] : 1) << '\n';
        }
        break;

    case NTStatus::PrimeGenerated:
        std::cout << "PrimeGenerated\n";
        std::cout << res.primes.size() << '\n';
        for (auto &p : res.primes) {
            std::cout << p << ' ';
        }
        std::cout << '\n';
        break;

    case NTStatus::Failure:
        std::cout << "Failure\n";
        break;

    case NTStatus::InvalidInput:
        std::cout << "InvalidInput\n";
        break;
    }
}

NTResult ParallelAnalysis(long_int n) {
    std::mutex mtx;
    NTResult best{NTStatus::Failure};

    auto update = [&](NTResult r) {
        std::lock_guard<std::mutex> lock(mtx);

        auto better = [](NTStatus a, NTStatus b) {
            if (a == NTStatus::FactorFound)
                return true;
            if (a == NTStatus::Prime && b != NTStatus::FactorFound)
                return true;
            if (a == NTStatus::ProbablyPrime && b == NTStatus::Failure)
                return true;
            return false;
        };

        if (better(r.status, best.status))
            best = std::move(r);
    };

    auto f1 =
        std::async(std::launch::async, [&]() { update(Trial_Compositeness_Check(n, 1000000)); });

    auto f2 = std::async(std::launch::async, [&]() { update(Strong_Pseudoprime(n, n)); });

    auto f3 = std::async(std::launch::async, [&]() { update(Pollard_Rho(n, 1, 10000)); });

    auto f4 = std::async(std::launch::async, [&]() { update(Pollard_P_Optimised(n, 1, 10000)); });

    auto f5 = std::async(std::launch::async, [&]() { update(CFRAC(n, 10000)); });

    f1.get();
    f2.get();
    f3.get();
    f4.get();
    f5.get();

    return best;
}

NTResult Run() {
    int mode;

    if (!(std::cin >> mode)) {
        return {NTStatus::InvalidInput};
    }

    switch (mode) {

    case 1: {
        long_int n;
        if (!(std::cin >> n)) {
            return {NTStatus::InvalidInput};
        }

        return ParallelAnalysis(n);
    }

    case 2: {
        long_int base, mod;
        if (!(std::cin >> base >> mod)) {
            return {NTStatus::InvalidInput};
        }

        return Strong_Pseudoprime(base, mod);
    }

    case 3: {
        long_int n;
        if (!(std::cin >> n)) {
            return {NTStatus::InvalidInput};
        }
        return Mersenne_Primes(n);
    }

    case 4: {
        long_int f;
        uint64_t max = 1000000;

        if (!(std::cin >> f)) {
            return {NTStatus::InvalidInput};
        }
        if (!(std::cin >> max)) {
            max = 1000000;
        }

        return Trial_Division(f, max);
    }

    case 5: {
        long_int n, c = 1, max = 10000;

        if (!(std::cin >> n)) {
            return {NTStatus::InvalidInput};
        }
        if (!(std::cin >> c)) {
            c = 1;
        }
        if (!(std::cin >> max)) {
            max = 10000;
        }

        return Pollard_Rho(n, c, max);
    }

    case 6: {
        uint64_t limit;

        if (!(std::cin >> limit)) {
            return {NTStatus::InvalidInput};
        }

        return segmented_sieve(limit);
    }

    case 7: {
        long_int n, c = 1;
        int64_t limit = 10000;

        if (!(std::cin >> n)) {
            return {NTStatus::InvalidInput};
        }
        if (!(std::cin >> c)) {
            c = 1;
        }
        if (!(std::cin >> limit)) {
            limit = 10000;
        }

        return Pollard_P_Optimised(n, c, limit);
    }

    case 8: {
        long_int n, b = 2;
        size_t k;

        if (!(std::cin >> n >> b >> k)) {
            return {NTStatus::InvalidInput};
        }

        std::vector<long_int> primes(k);
        for (auto &p : primes)
            if (!(std::cin >> p)) {
                return {NTStatus::InvalidInput};
            }

        return Lucas_Prime_Test(n, primes, b);
    }

    case 9: {
        long_int n;

        if (!(std::cin >> n)) {
            return {NTStatus::InvalidInput};
        }

        return Pepin_Test(n);
    }

    case 10: {
        long_int n;
        int limit = 10000;

        if (!(std::cin >> n)) {
            return {NTStatus::InvalidInput};
        }

        if (!(std::cin >> limit)) {
            limit = 10000;
        }

        return CFRAC(n, limit);
    }

    default:
        return {NTStatus::InvalidInput};
    }
}

} // namespace

int main() {
    NTResult res = Run();
    Print_Status(res);

    return 0;
}
