//#include "surrogate_components.h"
#include "../common/cli_options.h"
#include "fit_surrogate_triangle.h"

namespace Delta2 {
    namespace model {

        template<typename real>
        void divideAndFit(const std::vector<common::Triangle<real>>& tris, const std::vector<int>& tri_ids, int branching, bool force_branching_factor, common::Options& opt, std::vector<common::Triangle<real>>& tris_out, std::vector<std::vector<int>>& tri_ids_out) {
            if (branching > tri_ids.size()) {
                throw std::invalid_argument("There must be more triangles than branches");
            }
            tris_out.clear();
            for (int i = 0; i < branching; i++) {
                common::Triangle<real> t;
                tris_out.push_back(t);
            }

            tri_ids_out.clear();
            tri_ids_out.resize(branching);

            // Initialise the group triangles from individual triangles in the input
            for (int tri_i = 0; tri_i < branching; tri_i++) {
                tris_out[tri_i] = tris[tri_i];
            }

            // Iterations like k-mean algorithm
            for (int iteration = 0; iteration < 20; iteration++) {
                for (int grp_i = 0; grp_i < branching; grp_i++) {
                    tri_ids_out[grp_i].clear();
                }

                // Sort each triangle of the mesh into one of the groups based on how far the farthest vertex of that triangle is from the group triangle
                for (int tri_i : tri_ids) {
                    real best_dist = std::numeric_limits<real>::infinity();
                    int best_group;

                    for (int grp_i = 0; grp_i < branching; grp_i++) {
                        common::Triangle<real> grp_tri = tris_out[grp_i];
                        real da = grp_tri.distToPoint(tris[tri_i].A);
                        real db = grp_tri.distToPoint(tris[tri_i].B);
                        real dc = grp_tri.distToPoint(tris[tri_i].C);
                        real max = std::max(std::max(da, db), dc);
                        if (max < best_dist) {
                            best_group = grp_i;
                            best_dist = max;
                        }
                    }
                    tri_ids_out[best_group].push_back(tri_i);
                }
                // For each group compute a fitting triangle of all the members of it's group
                for (int grp_i = 0; grp_i < branching; grp_i++) {
                    int grp_size = tri_ids_out[grp_i].size();
                    if (grp_size > 0) {
                        std::vector<Eigen::Vector<real, 3>> pts;
                        for (int g = 0; g < grp_size; g++) {
                            common::Triangle<real> tri = tris[tri_ids_out[grp_i][g]];

                            pts.push_back(tri.A);
                            pts.push_back(tri.B);
                            pts.push_back(tri.C);
                        }
                        
                        common::Triangle<real> selected_tri;
                        real selected_eps;
                        selectSurrogateTriangle(pts, opt, selected_tri, selected_eps);
                        tris_out[grp_i] = selected_tri;
                    }
                }

                // There is a chance groups have no members.  If force_branching_factor is true then these groups will be populated
                if (force_branching_factor) {
                    bool improving = true;
                    bool branching_met = false;
                    while (improving && !branching_met) {
                        int largest_group_ind;
                        int largest_group_size = 0;
                        std::vector<int> zero_groups;
                        // Identify the largest group (to split) and the groups that have no members
                        for (int grp_i = 0; grp_i < branching; grp_i++) {
                            if (tri_ids_out[grp_i].size() > largest_group_size) {
                                largest_group_size = tri_ids_out[grp_i].size();
                                largest_group_ind = grp_i;
                            }
                            if (tri_ids_out[grp_i].size() == 0) {
                                zero_groups.push_back(grp_i);
                            }
                        }
                        if (zero_groups.size() > 0) {
                            // Perform the divide and fit operation but just on the triangle of the largest group.
                            std::vector<common::Triangle<real>> sub_clustered;
                            std::vector<std::vector<int>> sub_group_inds;
                            int num_grps = std::min(tri_ids_out[largest_group_ind].size(), zero_groups.size()+1);
                            divideAndFit(tris, tri_ids_out[largest_group_ind], num_grps, false, opt, sub_clustered, sub_group_inds);

                            // Copy the new groups up to the original group and the empty groups
                            int sub_zero_groups = 0;
                            for (int sub_grp_i = 0; sub_grp_i < sub_clustered.size(); sub_grp_i++) {
                                int grp_ind;
                                if (sub_grp_i == 0) {
                                    grp_ind = largest_group_ind;
                                } else {
                                    grp_ind = zero_groups[sub_grp_i-1];
                                }
                                if (sub_group_inds[sub_grp_i].size() > 0) {
                                    sub_zero_groups++;
                                    tri_ids_out[grp_ind].clear();
                                    for (int tri_ind : sub_group_inds[sub_grp_i]) {
                                        tri_ids_out[grp_ind].push_back(tri_ind);
                                    }
                                    tris_out[grp_ind] = sub_clustered[sub_grp_i];
                                }
                            }
                            improving = sub_zero_groups > 1; // ie this split operation helped
                        } else {
                            branching_met = true;  // We have enough groups
                        }
                    }

                    if (!branching_met) {
                        // Smarter k-mean like clustering methods didn't work so just create one triangle groups
                        std::vector<int> zero_groups;
                        // Repeat the process until all groups have at least one triangle in them
                        do {
                            int largest_group_ind;
                            int largest_group_size = 0;
                            zero_groups.clear();

                            // Find the largest group and any groups that are empty
                            for (int grp_i = 0; grp_i < branching; grp_i++) {
                                if (tri_ids_out[grp_i].size() > largest_group_size) {
                                    largest_group_size = tri_ids_out[grp_i].size();
                                    largest_group_ind = grp_i;
                                }
                                if (tri_ids_out[grp_i].size() == 0) {
                                    zero_groups.push_back(grp_i);
                                }
                            }
                            // Move individual triangles from the biggest group to new groups
                            for (int zero_i = 0; zero_i < std::min(zero_groups.size(), tri_ids_out[largest_group_ind].size()-1); zero_i++) {
                                int zero_ind = zero_groups[zero_i];
                                int tri_i = tri_ids_out[largest_group_ind].back();
                                tri_ids_out[zero_ind].push_back(tri_i);
                                tri_ids_out[largest_group_ind].pop_back();
                                tris_out[zero_ind] = tris[tri_i];
                            }
                            
                            // Recompute the surrogate triangle for the group that we took triangles from
                            int grp_size = tri_ids_out[largest_group_ind].size();
                            
                            std::vector<Eigen::Vector<real, 3>> pts;
                            for (int g = 0; g < grp_size; g++) {
                                common::Triangle<real> tri = tris[tri_ids_out[largest_group_ind][g]];

                                pts.push_back(tri.A);
                                pts.push_back(tri.B);
                                pts.push_back(tri.C);
                            }
                            
                            common::Triangle<real> selected_tri;
                            real selected_eps;
                            selectSurrogateTriangle(pts, opt, selected_tri, selected_eps);
                            tris_out[largest_group_ind] = selected_tri;
                        } while (zero_groups.size() > 0);
                    }
                }
            }
        }
        template void divideAndFit(const std::vector<common::Triangle<float>>&, const std::vector<int>&, int, bool, common::Options&, std::vector<common::Triangle<float>>&, std::vector<std::vector<int>>&);    
        template void divideAndFit(const std::vector<common::Triangle<double>>&, const std::vector<int>&, int, bool, common::Options&, std::vector<common::Triangle<double>>&, std::vector<std::vector<int>>&);    

        template<typename real>
        void calcEpss(const common::Triangle<real>& surrogate, const std::vector<common::Triangle<real>>& tris, const std::vector<int>& tri_ids, real& inner_eps_out, real& eps_out) {
            inner_eps_out = std::numeric_limits<real>::infinity();
            eps_out = 0.0;
            for (int id : tri_ids) {
                const common::Triangle<real>& tri = tris[id];
                // TODO this isn't strictly correct.  Imagine a small surrogate with a large child.  The closest point is in the middle of the child so gets ignored here.
                real a = surrogate.distToPoint(tri.A);
                real b = surrogate.distToPoint(tri.B);
                real c = surrogate.distToPoint(tri.C);
                inner_eps_out = std::min(inner_eps_out, std::min(a, std::min(b, c)));
                eps_out = std::max(eps_out, std::max(a, std::max(b, c)));
            }
        }
        template void calcEpss(const common::Triangle<float>&, const std::vector<common::Triangle<float>>&, const std::vector<int>&, float&, float&);
        template void calcEpss(const common::Triangle<double>&, const std::vector<common::Triangle<double>>&, const std::vector<int>&, double&, double&);

    }
}
