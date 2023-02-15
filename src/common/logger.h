#include <stdarg.h>
#include <stdio.h>
#include <mutex>

namespace Delta2 {
    namespace common {
        class Logger {
            public:
            Logger();
            void printf(int prio, const char* fmt, ...);
            void disablePrinting();
            int priority;
            private:
            std::mutex print_lock;
            bool disable;
        };
    }
}
