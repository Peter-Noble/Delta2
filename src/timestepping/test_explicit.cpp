#include "../testing/unit_test.h"
#include "../common/utils.h"
#include "../common/cli_options.h"
#include "../model/particle.h"
#include "../common/primitive_geo.h"
#include "explicit.h"
#include "../common/viewer.h"

TEST_CASE("Explicit free fall", "[grp_timestepping][explicit]") {
    REQUIRE( false );
}

TEST_CASE("Explicit cube on plane", "[grp_timestepping][explicit]") {
    REQUIRE( false );
}

TEST_CASE("Explicit cube stack", "[grp_timestepping][explicit]") {
    REQUIRE( false );
}

TEST_CASE("Explicit shouldn't segfault", "[grp_timestepping][explicit]") {
    Delta2::common::Options opt;
    std::vector<Delta2::Particle> particles;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    {
        Delta2::common::sphere(1.0, 5, V, F);

        std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, opt));
        auto& p = particles.emplace_back(M, 1.0, 1.0, 0.1);
        p.current_state.setTime(0.24850000000000019);
        p.current_state.setRotation({-5.2671048623388167e-17, -1.2652746144012484e-29, -7.4270434507300396e-33, 1.0});
        p.current_state.setAngular({-1.7289024151270123e-13, -4.3638715565621033e-26, -1.3451779585193522e-29});
        p.current_state.setTranslation({7.9031482637351628e-30, -3.3448847847070995e-17, 1.196805150000003});
        p.current_state.setVelocity({1.5806296527467047e-26, -6.3788570826964973e-14, -2.4352999999982718});
        p.future_state = p.current_state;
    }

    {
        Delta2::common::plane(2.0, V, F);
        std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, opt, true));
        auto& p = particles.emplace_back(M, 1.0, 1.0, 0.1);
        p.current_state.setTime(0.24850000000000019);
        p.current_state.setRotation({0.043619387365336, 0.0, 0.0, 0.9990482215818578});
        p.current_state.setAngular({0.0, 0.0, 0.0});
        p.current_state.setTranslation({0.0, 0.0, 0.0});
        p.current_state.setVelocity({0.0, 0.0, 0.0});
        p.future_state = p.current_state;
        p.is_static = true;
    }

    Delta2::timestepping::ExplicitScheme ts(&particles, [](const Delta2::Particle& p) {
        Eigen::Vector3d f = p.getMass() * Eigen::Vector3d({0, 0, -9.8});
        return f; 
    }, opt);

    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> view_draws;
    ts.step(0.0005, view_draws);
}


TEST_CASE("Explicit shouldn't segfault 2", "[grp_timestepping][explicit]") {
    Delta2::common::Options opt;
    std::vector<Delta2::Particle> particles;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    {
        Delta2::common::sphere(1.0, 5, V, F);

        std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, opt));
        auto& p = particles.emplace_back(M, 1.0, 1.0, 0.1);
        p.current_state.setTranslation({0.0, 0.0, 1.5});
    }

    {
        Delta2::common::plane(2.0, V, F);
        std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, opt, true));
        auto& p = particles.emplace_back(M, 1.0, 1.0, 0.1);
        p.is_static = true;
        p.current_state.setRotation(Delta2::common::eulerAnglesToQuaternion(Eigen::Vector3d({5.0/360.0*igl::PI*2.0, 0, 0})));
    }

    //Delta2::common::AnimationViewer view(&particles);
    Delta2::timestepping::ExplicitScheme ts(&particles, [](const Delta2::Particle& p) {
        Eigen::Vector3d f = p.getMass() * Eigen::Vector3d({0, 0, -9.8});
        return f; 
    }, opt);

    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> empty;
    //view.recordFrame(empty);

    for (int step = 0; step < 400; step++) {
        std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> view_draws;
        printf("Step: %i\n", step);
        for (int sub_step = 0; sub_step < 5; sub_step++) {
            ts.step(0.0005, view_draws);
        }
        //view.recordFrame(view_draws);
    }

    //view.show();
}
