#include <iostream>
#include <boost/multiprecision/gmp.hpp>

using boost::multiprecision::mpz_int;

mpz_int Gcd(mpz_int&& a, mpz_int&& b) {
    while (b != 0) {
        mpz_int temp = b;
        b = a % b;
        a = temp;
    }
    if (a >= 0) {
        return a;
    }
    else {
        return -a;
    }
}

mpz_int Safe_Gcd(const mpz_int& safe_a, const mpz_int& safe_b) {
    mpz_int a = safe_a;
    mpz_int b = safe_b;
    return Gcd(std::move(a), std::move(b));
}
mpz_int Fast_Pow(mpz_int&& num, mpz_int&& pow, const mpz_int& mod) {
    mpz_int res = 1;
    num = num % mod;
    while (pow > 0) {
        if (pow % 2 == 1) {
            res = (res * num) % mod;
        }
        num = num * num % mod;
        pow /= 2;
    }
    return res;
}

mpz_int Safe_Fast_Pow(const mpz_int& safe_num, const mpz_int& safe_pow, const mpz_int& mod) {
    mpz_int num = safe_num;
    mpz_int pow = safe_pow;
    return Fast_Pow(std::move(num), std::move(pow), mod);
}

mpz_int Fast_Pow_No_Mod(mpz_int&& num, mpz_int&& pow) {
    mpz_int res = 1;
    while (pow > 0) {
        if (pow % 2 == 1) {
            res = res * num;
        }
        num = num * num;
        pow /= 2;
    }
    return res;
}

mpz_int Safe_Fast_Pow_No_Mod(const mpz_int& num, const mpz_int& pow) {
    mpz_int safe_num = num;
    mpz_int safe_pow = pow;
    return Fast_Pow_No_Mod(std::move(safe_num), std::move(safe_pow));
}

bool Strong_Pseudoprime(mpz_int&& base, const mpz_int& mod) {
    mpz_int t = mod - 1;
    mpz_int a = 0;
    while (t % 2 == 0) {
        t /= 2;
        ++a;
    }
    mpz_int test = Fast_Pow(std::move(base), std::move(t), mod);
    if (test == 1 || test == mod - 1) {
        return true;
    }
    for (mpz_int i = 0; i < a; ++i) {
        test = test * test % mod;
        if (test == mod - 1) {
            return true;
        }
    }
    return false;
    //вероятно стоит вместе с boolом выводить результат,
    //чтобы пользователь понял, что именно получилось
}

bool Safe_Strong_Pseudoprime(const mpz_int& safe_base, const mpz_int& mod) {
    mpz_int base = safe_base;
    return Strong_Pseudoprime(std::move(base), mod);
}

bool Safe_Mersenne_Primes(const mpz_int& n) {
    mpz_int m = Safe_Fast_Pow_No_Mod(2, n);
    mpz_int s = 4;
    for (mpz_int i = 2; i <= n - 1; ++i) {
        s = (s * s - 2) % m;
        if (s == 0) {
            return true;
        }
    }
    return false;
}

//true если простое, иначе false и чекать разложение в векторах
//опять же, мб стоит вывести это в консоль для понимания
bool Trial_Division(mpz_int&& f, const mpz_int& max, std::vector<mpz_int>& primes, std::vector<mpz_int>& degs) {
    long long i = 0;
    for (int k = 2; k < 4; ++k) {
        if (f % k == 0) {
            primes.push_back(k);
            degs.push_back(1);
            f /= k;
            while (f % k == 0) {
                ++degs[i];
                f /= k;
            }
            ++i;
        }
    }

    mpz_int d = 5;
    int add = 2;

    while (d <= max && d * d <= f) {
        if (f % d == 0) {
            primes.push_back(d);
            degs.push_back(1);
            f /= d;
            while (f % d == 0) {
                ++degs[i];
                f /= d;
            }
            ++i;
        }
        d += add;
        add = 6 - add;
    }

    if (d * d > f) {
        primes.push_back(f);
        degs.push_back(1);
        ++i;
    }

    if (i == 1) {
        return true;
    } else {
        return false;
    }
}

bool Safe_Trial_Division(const mpz_int& f, const mpz_int& max, std::vector<mpz_int>& primes, std::vector<mpz_int>& degs) {
    mpz_int safe_f = f;
    return Trial_Division(std::move(safe_f), max, primes, degs);
}

//возвращает найденный делитель n или 0 если неудача
//сначала strong_pseudoprime!
mpz_int Safe_Pollard_Rho(const mpz_int& n, const mpz_int& c, const mpz_int& max = 10000) {
    mpz_int x1 = 2;
    mpz_int x2 = c + 4;
    mpz_int range = 1;
    mpz_int terms = 0;
    mpz_int product = 1;
    while (terms <= max) {
        for (mpz_int j = 1; j < range && terms < max; ++j) {
            x2 = (x2 * x2 + c) % n;
            product *= (x1 - x2) % n;
            ++terms;
            if (terms % 10 == 0) {
                mpz_int g = Safe_Gcd(n, product);
                if (g > 1 && g < n) {
                    return g;
                }
                if (g == n) {
                    return 0;
                }
                product = 1;
            }
        }
        x1 = x2;
        range *= 2;
    }
    return 0;
}

const int64_t L1D_CACHE_SIZE = 32768;

std::vector<int64_t> segmented_sieve(int64_t limit)
{
    int64_t sqrt = (int64_t) std::sqrt(limit);
    int64_t segment_size = std::max(sqrt, L1D_CACHE_SIZE);
    int64_t count = (limit < 2) ? 0 : 1;

    int64_t i = 3;
    int64_t n = 3;
    int64_t s = 3;

    std::vector<char> sieve(segment_size);
    std::vector<char> is_prime(sqrt + 1, true);
    std::vector<int64_t> primes;
    std::vector<int64_t> multiples;

    for (int64_t low = 0; low <= limit; low += segment_size)
    {
        std::fill(sieve.begin(), sieve.end(), true);
        int64_t high = low + segment_size - 1;
        high = std::min(high, limit);


        for (; i * i <= high; i += 2) {
            if (is_prime[i]) {
                for (int64_t j = i * i; j <= sqrt; j += i) {
                    is_prime[j] = false;
                }
            }
        }

        for (; s * s <= high; s += 2)
        {
            if (is_prime[s])
            {
                primes.push_back(s);
                multiples.push_back(s * s - low);
            }
        }

        for (std::size_t i = 0; i < primes.size(); i++)
        {
            int64_t j = multiples[i];
            for (int64_t k = primes[i] * 2; j < segment_size; j += k) {
                sieve[j] = false;
            }
            multiples[i] = j - segment_size;
        }

        for (; n <= high; n += 2) {
            if (sieve[n - low]) {
                count++;
            }
        }
    }

    return primes;
}

//вероятно не optimised будет быстрее при малом limit
mpz_int Safe_Pollard_P_Optimised(const mpz_int& n, const mpz_int& c, int64_t limit) {
    mpz_int m = c % n;
    std::vector<int64_t> primes = segmented_sieve(limit);

    for (int64_t p : primes) {
        int64_t e = std::floor(std::log(limit) / std::log(p));
        mpz_int pe = Safe_Fast_Pow_No_Mod(p, e);
        m = Safe_Fast_Pow(m, pe, n);
    }

    mpz_int g = gcd(m - 1, n);
    if (g > 1 && g < n) {
        return g;
    }
    return 0;

}

bool Lucas_Prime_Test(const mpz_int& n, const std::vector<mpz_int>& prime_factors, const mpz_int& b) {
    mpz_int k = n - 1;

    if (Safe_Fast_Pow(b, k, n) != 1) {
        return false;
    }

    for (const auto& p : prime_factors) {
        if (Safe_Fast_Pow(b, k / p, n) == 1) {
            return false;
        }
    }

    return true;
}

bool Pepin_Test(mpz_int&& n) {
    mpz_int num = Fast_Pow_No_Mod(2, std::move(n)) - 1;

    mpz_int exp = (num - 1) / 2;
    mpz_int res = Fast_Pow(3, std::move(exp), num);

    return (res == num - 1);
}

bool Safe_Pepin_Test(const mpz_int& n) {
    mpz_int safe_n = n;
    return Pepin_Test(std::move(safe_n));
}

//только 1 делитель, если нашли, иначе 0
mpz_int CFRAC(mpz_int&& n, int limit = 10000) {
    mpz_int sqrt = boost::multiprecision::sqrt(n);
    mpz_int a = sqrt;
    mpz_int b = 1;
    mpz_int prev_a = 1;
    mpz_int prev_b = 0;

    std::vector<mpz_int> x_list;
    std::vector<mpz_int> y_list;

    for (int i = 0; i < limit; ++i) {
        mpz_int ai = (sqrt + prev_a) / b;
        mpz_int next_a = ai * b - prev_a;
        mpz_int next_b = (n - next_a * next_a) / b;

        if (next_b == 0) {
            break;
        }

        mpz_int xi = next_a;
        mpz_int yi = (xi * xi) % n;
        x_list.push_back(xi);
        y_list.push_back(yi);
        mpz_int g = Gcd(std::move(yi), std::move(n));

        if (g > 1 && g < n) {
            return g;
        }

        prev_a = next_a;
        b = next_b;
    }

    return 0;
}

mpz_int Safe_CFRAC(const mpz_int& n, int limit = 10000) {
    mpz_int safe_n = n;
    return CFRAC(std::move(safe_n), limit);
}
int main() {

}