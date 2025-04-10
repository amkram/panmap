#include "coordinates.hpp"
#include "logging.hpp"
#include <iostream>

using namespace logging;

int main() {
    // Print header
    msg("Starting simple integration test");
    
    // Test coordinates basic functionality
    coordinates::tupleCoord_t coord(0, 0, 0, 0);
    
    if (coord.isValidBasic()) {
        msg("Basic coordinate validation passed");
        return 0;
    } else {
        err("Basic coordinate validation failed");
        return 1;
    }
} 