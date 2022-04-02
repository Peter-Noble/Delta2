#include "broad_phase_embree.h"
#include "continuous_comparison.h"

using namespace Delta2;

std::mutex lock;
std::mutex broad_phase_lock;
model::ParticleHandler* global_particles;
double time_step_size;


embree::Vec3f cast(Eigen::Vector3f v) {
    return embree::Vec3f(v[0], v[1], v[2]);
}

Eigen::Vector3f cast(embree::Vec3f v) {
    return Eigen::Vector3f(v[0], v[1], v[2]);
}

/*
 * We will register this error handler with the device in initializeDevice(),
 * so that we are automatically informed on errors.
 * This is extremely helpful for finding bugs in your code, prevents you
 * from having to add explicit error checking to each Embree API call.
 */
void errorFunction(void *userPtr, enum RTCError error, const char *str)
{
    printf("error %d: %s\n", error, str);
}

/*
 * Embree has a notion of devices, which are entities that can run 
 * raytracing kernels.
 * We initialize our device here, and then register the error handler so that 
 * we don't miss any errors.
 *
 * rtcNewDevice() takes a configuration string as an argument. See the API docs
 * for more information.
 *
 * Note that RTCDevice is reference-counted.
 */
RTCDevice initializeDevice()
{
    RTCDevice device = rtcNewDevice(NULL);

    if (!device)
        printf("error %d: cannot create device\n", rtcGetDeviceError(NULL));

    rtcSetDeviceErrorFunction(device, errorFunction, NULL);
    return device;
}

embree::BBox3fa translate(embree::BBox3fa box, Eigen::Vector3f translation) {
    embree::BBox3fa result = box;
    result.lower += cast(translation);
    result.upper += cast(translation);
    return result;
}

void boundsFunc(const struct RTCBoundsFunctionArguments *args)
{
    void *ptr = args->geometryUserPtr;
    unsigned geomID = (unsigned)(size_t)ptr;
    Delta2::Particle &p = (*global_particles)[geomID];
    // triangle_bvh &p_t = globals::particle_bvhs[geomID];
    
    // embree::BBox3fa bounds = embree::empty;  // TODO why doesn't this work with -O0?
    embree::BBox3fa bounds = embree::BBox3fa(embree::EmptyTy());
    embree::BBox3fa mesh_aabb;
    if (p.is_static) {
        mesh_aabb = p.mesh->getAABB(p.current_state.getTransformation());
        bounds.extend(mesh_aabb);
    } else {
        mesh_aabb = p.mesh->getSphereBounds();
        Eigen::Vector3f start_loc = p.current_state.getTranslation().cast<float>();
        bounds.extend(translate(mesh_aabb, start_loc));
        Eigen::Vector3f end_loc = p.current_state.getTranslation().cast<float>() + p.current_state.getVelocity().cast<float>() * time_step_size;
        bounds.extend(translate(mesh_aabb, end_loc));
    }
    bounds = bounds.enlarge_by(p.geo_eps);

    // args->primID;  // At the moment there is only one prim per object

    *(embree::BBox3fa *)args->bounds_o = bounds;
}

void boundsFuncIndividualTime(const struct RTCBoundsFunctionArguments *args)
{
    void *ptr = args->geometryUserPtr;
    unsigned geomID = (unsigned)(size_t)ptr;
    Delta2::Particle &p = (*global_particles)[geomID];
    // triangle_bvh &p_t = globals::particle_bvhs[geomID];
    
    // embree::BBox3fa bounds = embree::empty;  // TODO why doesn't this work with -O0?
    embree::BBox3fa bounds = embree::BBox3fa(embree::EmptyTy());
    embree::BBox3fa mesh_aabb;
    if (p.is_static) {
        mesh_aabb = p.mesh->getAABB(p.current_state.getTransformation());
        bounds.extend(mesh_aabb);
    } else {
        mesh_aabb = p.mesh->getSphereBounds();
        // Eigen::Vector3f start_loc = p.current_state.getTranslation().cast<float>();
        Eigen::Vector3f start_loc = p.last_state.getTranslation().cast<float>();
        bounds.extend(translate(mesh_aabb, start_loc));
        // Eigen::Vector3f end_loc = p.current_state.getTranslation().cast<float>() + p.current_state.getVelocity().cast<float>() * p.last_time_step_size;
        Eigen::Vector3f end_loc = p.future_state.getTranslation().cast<float>();
        bounds.extend(translate(mesh_aabb, end_loc));
    }
    bounds = bounds.enlarge_by(p.geo_eps);

    // args->primID;  // At the moment there is only one prim per object

    *(embree::BBox3fa *)args->bounds_o = bounds;
}

unsigned int createUserGeom(RTCDevice device, RTCScene scene, Delta2::Particle &p, int ind)
{
    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
    unsigned int geomID = rtcAttachGeometry(scene, geom);
    rtcSetGeometryUserPrimitiveCount(geom, 1);
    rtcSetGeometryUserData(geom, (void *)(size_t)ind); // put the particle index in here
    rtcSetGeometryBoundsFunction(geom, boundsFunc, nullptr);
    rtcSetGeometryIntersectFunction(geom, nullptr);
    rtcCommitGeometry(geom);
    rtcReleaseGeometry(geom);
    return geomID;
}

unsigned int createUserGeomIndividualTime(RTCDevice device, RTCScene scene, Delta2::Particle &p, int ind)
{
    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
    unsigned int geomID = rtcAttachGeometry(scene, geom);
    rtcSetGeometryUserPrimitiveCount(geom, 1);
    rtcSetGeometryUserData(geom, (void *)(size_t)ind); // put the particle index in here
    rtcSetGeometryBoundsFunction(geom, boundsFuncIndividualTime, nullptr);
    rtcSetGeometryIntersectFunction(geom, nullptr);
    rtcCommitGeometry(geom);
    rtcReleaseGeometry(geom);
    return geomID;
}

void collideFunc(void *userPtr, RTCCollision *collisions, unsigned int num_collisions)
{
    for (size_t i = 0; i < num_collisions;)
    {
        bool intersect;

        Delta2::Particle &a = (*global_particles)[collisions[i].geomID0];
        Delta2::Particle &b = (*global_particles)[collisions[i].geomID1];

        if (collisions[i].geomID0 < collisions[i].geomID1 && (!a.is_static || !b.is_static)) {
            double t = Delta2::collision::sphereSphereFirstIntersect(a.current_state.getTranslation(), a.current_state.getVelocity(), a.mesh->getBoundingRadius() + a.geo_eps,
                                                  b.current_state.getTranslation(), b.current_state.getVelocity(), b.mesh->getBoundingRadius() + b.geo_eps);

            intersect = t >= 0.0 && t < time_step_size;
            // Eigen::Vector3d diff = (a.translation + a.linear_velocity * globals::opt.time_step_size) - (b.translation + b.linear_velocity * globals::opt.time_step_size);
            // intersect = diff.norm() < a.getBoundingSphereRadius() + b.getBoundingSphereRadius() + 2 * globals::opt.geo_eps;
        } else {
            intersect = false;
        }

        if (intersect) {
            i++;
        } else {
            collisions[i] = collisions[--num_collisions];
        }
    }

    if (num_collisions == 0)
        return;

    std::lock_guard<std::mutex> guard(lock);
    for (size_t i = 0; i < num_collisions; i++)
    {
        const unsigned geomID0 = collisions[i].geomID0;
        const unsigned primID0 = collisions[i].primID0;
        const unsigned geomID1 = collisions[i].geomID1;
        const unsigned primID1 = collisions[i].primID1;

        static_cast<Delta2::collision::BroadPhaseCollisions *>(userPtr)->push_back(std::make_pair(std::make_pair(geomID0, primID0), std::make_pair(geomID1, primID1)));
    }
}

void collideFuncIndividualTime(void *userPtr, RTCCollision *collisions, unsigned int num_collisions)
{
    for (size_t i = 0; i < num_collisions;)
    {
        bool intersect;

        Delta2::Particle &a = (*global_particles)[collisions[i].geomID0];
        Delta2::Particle &b = (*global_particles)[collisions[i].geomID1];

        if (collisions[i].geomID0 < collisions[i].geomID1 && (!a.is_static || !b.is_static)) {
            Delta2::State a_state = a.current_state;
            Delta2::State b_state = b.current_state;
            
            if (!a.is_static && !b.is_static) {
                if (a.current_state.getTime() < b.current_state.getTime()) {
                    b_state = b.current_state.interpolate(b.last_state, a.current_state.getTime());
                }
                else {
                    a_state = a.current_state.interpolate(a.last_state, b.current_state.getTime());
                }
            }

            double t = Delta2::collision::sphereSphereFirstIntersect(a_state.getTranslation(), a_state.getVelocity(), a.mesh->getBoundingRadius() + a.geo_eps,
                                                  b_state.getTranslation(), b_state.getVelocity(), b.mesh->getBoundingRadius() + b.geo_eps);

            intersect = t >= 0.0 && t < std::max(a.last_time_step_size, b.last_time_step_size);
            // Eigen::Vector3d diff = (a.translation + a.linear_velocity * globals::opt.time_step_size) - (b.translation + b.linear_velocity * globals::opt.time_step_size);
            // intersect = diff.norm() < a.getBoundingSphereRadius() + b.getBoundingSphereRadius() + 2 * globals::opt.geo_eps;
        } else {
            intersect = false;
        }

        if (intersect) {
            i++;
        } else {
            collisions[i] = collisions[--num_collisions];
        }
    }

    if (num_collisions == 0)
        return;

    std::lock_guard<std::mutex> guard(lock);
    for (size_t i = 0; i < num_collisions; i++)
    {
        const unsigned geomID0 = collisions[i].geomID0;
        const unsigned primID0 = collisions[i].primID0;
        const unsigned geomID1 = collisions[i].geomID1;
        const unsigned primID1 = collisions[i].primID1;

        static_cast<Delta2::collision::BroadPhaseCollisions *>(userPtr)->push_back(std::make_pair(std::make_pair(geomID0, primID0), std::make_pair(geomID1, primID1)));
    }
}

Delta2::collision::BroadPhaseCollisions Delta2::collision::broadPhaseEmbree(model::ParticleHandler& particles, double t)
{
    // This is needed because of the global particle store and time_step_size
    std::lock_guard<std::mutex> guard(broad_phase_lock);
    time_step_size = t;

    global_particles = &particles;
    BroadPhaseCollisions sim_collisions;

    RTCDevice device = initializeDevice();
    RTCScene scene = rtcNewScene(device);

    for (int p = 0; p < global_particles->size(); p++) {
        createUserGeom(device, scene, (*global_particles)[p], p);
    }

    rtcCommitScene(scene);

    rtcCollide(scene, scene, collideFunc, &sim_collisions);

    rtcReleaseScene(scene);
    rtcReleaseDevice(device);

    return sim_collisions;
}

Delta2::collision::BroadPhaseCollisions Delta2::collision::broadPhaseEmbree(model::ParticleHandler& particles)
{
    // This is needed because of the global particle store and time_step_size
    std::lock_guard<std::mutex> guard(broad_phase_lock);

    global_particles = &particles;
    BroadPhaseCollisions sim_collisions;

    RTCDevice device = initializeDevice();
    RTCScene scene = rtcNewScene(device);

    for (int p = 0; p < global_particles->size(); p++) {
        createUserGeomIndividualTime(device, scene, (*global_particles)[p], p);
    }

    rtcCommitScene(scene);

    rtcCollide(scene, scene, collideFuncIndividualTime, &sim_collisions);

    rtcReleaseScene(scene);
    rtcReleaseDevice(device);

    return sim_collisions;
}
