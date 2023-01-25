#include <stdarg.h>
#include <stdio.h>
#include <mutex>

namespace Delta2 {
    namespace common {
        class Logger {
            public:
            Logger();
            void printf(const char* fmt, ...);
            void disablePrinting();
            private:
            std::mutex print_lock;
            bool disable;
        };
    }
}
