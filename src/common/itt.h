#pragma once

#include <ittnotify.h>

namespace Delta2 {
    namespace common {
        struct IttHandles {
            __itt_domain* detailed_domain;
            __itt_domain* phases_domain;
            __itt_domain* step_domain;

            IttHandles() {
                printf("IttHandles init\n");
                detailed_domain = __itt_domain_create("Delta2.detailed");
                phases_domain = __itt_domain_create("Delta2.phases");
                step_domain = __itt_domain_create("Delta2.step");
            }

            void disable_detailed_domain() {
                if (detailed_domain != 0x0) {
                    detailed_domain->flags = 0;
                }
            }

            void disable_phases_domain() {
                if (phases_domain != 0x0) {
                    phases_domain->flags = 0;
                }
            }

            void disable_step_domain() {
                if (step_domain != 0x0) {
                    step_domain->flags = 0;
                }

            }
        };
    }
}
