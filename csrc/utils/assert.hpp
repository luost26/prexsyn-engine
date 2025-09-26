#pragma once

#define Ensures(condition)                                                     \
    do {                                                                       \
        if (!(condition)) {                                                    \
            throw std::runtime_error("Precondition failed: " #condition        \
                                     ", at " __FILE__ ":" +                    \
                                     std::to_string(__LINE__));                \
        }                                                                      \
    } while (0)
