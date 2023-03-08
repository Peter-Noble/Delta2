#pragma once

#include "../model/particle.h"
#include <vector>

namespace Delta2 {
    class D2Writer {
        public:
        D2Writer(std::vector<Delta2::Particle>& particles);
        void capture(std::vector<Delta2::Particle>& particles);
        void write(std::vector<Delta2::Particle>& particles, int frame_rate, std::string file_name);
        private:
        std::vector<std::vector<Delta2::State>> states;
    };
}
