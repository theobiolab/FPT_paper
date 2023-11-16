#include <iostream>
#include <gmpxx.h>

int main() {
    // Initialize GMP variables
    mpz_class a, b, result;

    // Set values for a and b
    a = "123456789012345678901234567890";
    b = "987654321012345678901234567890";

    // Perform arithmetic operations
    result = a + b;

    // Print the result
    std::cout << "Result: " << result << std::endl;

    return 0;
}
