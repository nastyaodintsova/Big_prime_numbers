#include "Sieves.h"

#include <algorithm>
#include <cmath>

// оптимизированное под размер l1-кэша решето Эратосфена, проверяющая простоту сегментами.
const constexpr uint64_t l1d_cache_size = 32768;

NTResult segmented_sieve(uint64_t limit) {
    uint64_t sqrt = static_cast<uint64_t>(std::sqrt(limit) + 1);
    uint64_t segment_size = std::max(sqrt, l1d_cache_size);

    std::vector<bool> is_prime(sqrt + 1, true);

    std::vector<uint64_t> base_primes;
    std::vector<uint64_t> multiples;

    std::vector<long_int> result_primes;

    for (uint64_t i = 2; i * i <= sqrt; ++i) {
        if (is_prime[i]) {
            for (uint64_t j = i * i; j <= sqrt; j += i) {
                is_prime[j] = false;
            }
        }
    }

    for (uint64_t i = 2; i <= sqrt; ++i) {
        if (is_prime[i]) {
            base_primes.push_back(i);
            multiples.push_back(i * i);
        }
    }

    std::vector<char> sieve(segment_size);

    for (uint64_t low = 0; low <= limit; low += segment_size) {
        std::fill(sieve.begin(), sieve.end(), true);
        uint64_t high = std::min(low + segment_size - 1, limit);
        for (size_t i = 0; i < base_primes.size(); ++i) {
            uint64_t p = base_primes[i];

            uint64_t start = std::max(p * p, (low + p - 1) / p * p);

            for (uint64_t j = start; j <= high; j += p) {
                sieve[j - low] = false;
            }
        }

        for (uint64_t n = std::max<uint64_t>(low, 2); n <= high; ++n) {
            if (sieve[n - low]) {
                result_primes.push_back(n);
            }
        }
    }

    return {NTStatus::PrimeGenerated, 0, result_primes, {}};
}
