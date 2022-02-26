#include <iostream>

#define CHECK(cond) \
    if (!(cond)) { \
        std::cerr << "CHECK FAILED: " #cond << \
        " \n" << __FUNCTION__ << " (" << __LINE__ << ")" << std::endl; \
        return EXIT_FAILURE; \
    }
