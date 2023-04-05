#pragma once

#include "contact.h"
#include "cluster.h"
#include "contact_bundle.h"

namespace Delta2 {
    namespace collision {
        std::vector<ContactBundle> bundle_hits(Cluster& cluster, std::vector<Contact<double>>& hits);
        std::vector<std::vector<ContactBundle>> colour_hits(Cluster& cluster, std::vector<Contact<double>>& hits);
    }
}
