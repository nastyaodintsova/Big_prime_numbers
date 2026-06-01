#include "Arithmetic.h"

FactorPower Extract_Factor_Power(long_int n, const long_int &d) {
    long_int k = 0;

    while (n % d == 0) {
        n /= d;
        ++k;
    }

    return {n, d, k};
}

long_int Gcd(long_int &&a, long_int &&b) {
    if (a < 0) {
        a = -a;
    }
    while (b != 0) {
        a %= b;
        std::swap(a, b);
    }
    return a;
}

long_int Gcd(const long_int &safe_a, const long_int &safe_b) {
    long_int a = safe_a;
    long_int b = safe_b;
    return Gcd(std::move(a), std::move(b));
}

long_int Fast_Pow(long_int &&num, long_int &&pow, const long_int &mod) {
    long_int res = 1;
    num %= mod;
    while (pow > 0) {
        if (pow % 2 == 1) {
            res = (res * num) % mod;
        }
        num *= num;
        num %= mod;
        pow /= 2;
    }
    return res;
}

long_int Fast_Pow(const long_int &safe_num, const long_int &safe_pow, const long_int &mod) {
    long_int num = safe_num;
    long_int pow = safe_pow;
    return Fast_Pow(std::move(num), std::move(pow), mod);
}

long_int Fast_Pow(long_int &&num, long_int &&pow) {
    long_int res = 1;
    while (pow > 0) {
        if (pow % 2 == 1) {
            res *= num;
        }
        num *= num;
        pow /= 2;
    }
    return res;
}

long_int Fast_Pow(const long_int &num, const long_int &pow) {
    long_int safe_num = num;
    long_int safe_pow = pow;
    return Fast_Pow(std::move(safe_num), std::move(safe_pow));
}
