#include <iostream>
#include <mpc.h>

int main() {
    // Initialize MPC variables
    mpc_t a, b, result;
    mpc_init2(a, 256);  // Initialize a with precision 256 bits
    mpc_init2(b, 256);
    mpc_init2(result, 256);

    // Set values for a and b
    mpc_set_str(a, "1.23", 10, MPC_RNDNN);  // Set a to 1.23
    mpc_set_d(b, 2.0, MPC_RNDNN);       // Set b to 2.0

    // Perform arithmetic operation (e.g., addition)
    mpc_add(result, a, b, MPC_RNDNN);

    // Print the result
    gmp_printf("Result: %.Ff\n", result);

    // Free allocated memory
    mpc_clear(a);
    mpc_clear(b);
    mpc_clear(result);

    return 0;
}
