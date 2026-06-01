cmake_minimum_required(VERSION 3.28)
project(Big_prime_numbers)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(Big_prime_numbers
        main.cpp

        Arithmetic.cpp
        PrimalityTests.cpp
        Factorization.cpp
        Sieves.cpp
)

target_include_directories(Big_prime_numbers PRIVATE
        include
)

target_link_libraries(Big_prime_numbers
        gmp
        gmpxx
)