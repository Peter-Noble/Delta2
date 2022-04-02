#pragma once

#include "../model/particle.h"
#include <functional>
#include "../common/cli_options.h"

namespace Delta2 {
    namespace timestepping {
        class TimestepScheme
        {
        public:
            TimestepScheme(std::vector<Delta2::Particle>* p, common::Options& opt);
            virtual ~TimestepScheme();
            virtual void stepSetup() = 0;
            virtual void broadPhase() = 0;
            virtual void narrowPhase() = 0;
            virtual void solveContacts() = 0;
            virtual void integrate(u_int32_t p_i, bool apply) = 0;
            virtual void friction() = 0;
            virtual void step();
            void setStepSize(double step);
            //const std::vector<Delta2::Particle>& getParticles();
            TimestepScheme(const TimestepScheme&) = delete;
            void operator=(const TimestepScheme&) = delete;
            std::function<Eigen::Vector3d(const Particle&)> external_force;
        protected:
            common::Options _opt;
            std::vector<Delta2::Particle>* _particles;
            double _step_size;
            std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& _view_draws;
        };

        std::function<Eigen::Vector3d(const Particle&)> default_gravity_callback = [](const Delta2::Particle& p) {
            Eigen::Vector3d f = p.getMass() * Eigen::Vector3d({0, 0, -9.8});
            return f; 
        };
    }
}
