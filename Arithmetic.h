#pragma once

#include "LongInt.h"

FactorPower Extract_Factor_Power(long_int n, const long_int &d);

long_int Gcd(long_int &&a, long_int &&b);
long_int Gcd(const long_int &a, const long_int &b);

long_int Fast_Pow(long_int &&num, long_int &&pow, const long_int &mod);

long_int Fast_Pow(const long_int &num, const long_int &pow, const long_int &mod);

long_int Fast_Pow(long_int &&num, long_int &&pow);

long_int Fast_Pow(const long_int &num, const long_int &pow);
