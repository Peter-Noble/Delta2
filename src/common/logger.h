#include <stdarg.h>
#include <stdio.h>
#include <mutex.h>

namespace Delta2 {
    namespace common {
        class Logger {
            public:
            Logger();
            void printf(char* fmt, ...);
            void disablePrinting();
            private:
            std::mutex print_lock;
            bool disable;
        };
    }
}
