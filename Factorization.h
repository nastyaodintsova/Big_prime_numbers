#pragma once

#include "LongInt.h"

NTResult Trial_Division(long_int &&f, const uint64_t &max);

NTResult Trial_Division(const long_int &f, const long_int &max);

NTResult Pollard_Rho(const long_int &n, const long_int &c, const long_int &max);

NTResult Pollard_P_Optimised(const long_int &n, const long_int &c, int64_t limit);

NTResult CFRAC(long_int &&n, int limit);

NTResult CFRAC(const long_int &n, int limit);
