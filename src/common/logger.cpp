#include "logger.h"

Delta2::common::Logger::Logger() {
    disable = false;
}

void Delta2::common::Logger::printf(const char* fmt, ...) {
    if (!disable) {
        std::lock_guard guard(print_lock);
        std::printf("<Logger> ");
        va_list args;
        va_start(args, fmt);
        vprintf(fmt, args);
        va_end(args);
    }
}

void Delta2::common::Logger::disablePrinting() {
    disable = true;
}
