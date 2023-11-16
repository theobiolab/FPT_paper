#!/bin/bash

if [ "$(uname)" = "Darwin" ]; then
    export DYLD_LIBRARY_PATH=$(pwd)/lib/gmp-6.3.0/.libs:$DYLD_LIBRARY_PATH 
    export DYLD_LIBRARY_PATH=$(pwd)/lib/mpfr-4.2.1/src/.libs:$DYLD_LIBRARY_PATH 
    export DYLD_LIBRARY_PATH=$(pwd)/lib/mpc-1.3.1/src/.libs:$DYLD_LIBRARY_PATH
elif [ "$(uname)" = "Linux" ]; then
    export LD_LIBRARY_PATH=$(pwd)/lib/gmp-6.3.0/.libs:$LD_LIBRARY_PATH 
    export LD_LIBRARY_PATH=$(pwd)/lib/mpfr-4.2.1/src/.libs:$LD_LIBRARY_PATH 
    export LD_LIBRARY_PATH=$(pwd)/lib/mpc-1.3.1/src/.libs:$LD_LIBRARY_PATH 
else
    echo "Unsupported operating system"
fi