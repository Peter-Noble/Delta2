#pragma once

#include "contact.h"
#include "../model/collision_state.h"
#include "comparison_dist.h"
#include "continuous_comparison.h"
// #include <ittnotify.h>
#include "contact_state.h"
#include "../strategies/contact_detection/contact_detection_strategy.h"

namespace Delta2 {
    namespace collision {
        template<typename real, int branching, int bucket_size>
        std::vector<Contact<real>> compareTreesFull(Particle& a, Particle& b, strategy::ContactDetectionStrategy& bucket_comparison) {
            std::vector<Contact<real>> result;
            
            std::vector<DeferredCompare> dc;
            std::vector<DeferredCompare> dc_next;
            dc.emplace_back(a.mesh->getSurrogateTree().getNode(0), b.mesh->getSurrogateTree().getNode(0));
            std::vector<int> num_hits;
            std::vector<Contact<real>> hits;
            while (!dc.empty()) {
                hits.clear();
                // findContactsBucket<branching, real>(dc, a, b, hits, num_hits);
                bucket_comparison.findContactsBucket(dc, a, b, hits, num_hits);
                int last_hit_num = 0;
                for (int i = 0; i < dc.size(); i++) {
                    if (num_hits[i] > last_hit_num) {
                        bool advance_a = false;
                        bool advance_b = false;
                        if (!dc[i].a.is_inner && !dc[i].b.is_inner) {
                            result.insert(result.end(), hits.begin() + last_hit_num, hits.begin() + num_hits[i]);
                        }
                        else if (!dc[i].a.is_inner) {
                            advance_b = true;
                        }
                        else if (!dc[i].b.is_inner) {
                            advance_a = true;
                        }
                        else {
                            // TODO this decision of which to advance first is a little arbitrary.
                            //    Would probably be better to advance the one with the largest area.
                            if (dc[i].a.depth < dc[i].b.depth) {
                                advance_a = true;
                            }
                            else {
                                advance_b = true;
                            }
                        }
                        if (advance_a) {
                            for (int c = 0; c < dc[i].a.num_children; c++) {
                                dc_next.emplace_back(a.mesh->getSurrogateTree().getNode(c + dc[i].a.child_id_first), dc[i].b);
                            }
                        }
                        else if (advance_b) {
                            for (int c = 0; c < dc[i].b.num_children; c++) {
                                dc_next.emplace_back(dc[i].a, b.mesh->getSurrogateTree().getNode(c + dc[i].b.child_id_first));
                            }
                        }
                    }
                    last_hit_num = num_hits[i];
                }
                dc.swap(dc_next);
                dc_next.clear();
            }
            return result;
        };

        template<typename real, int branching, int bucket_size>
        std::vector<Contact<real>> compareTreesFullCurrent(Particle& a, Particle& b, strategy::ContactDetectionStrategy& bucket_comparison) {
            std::vector<Contact<real>> result;
            
            std::vector<DeferredCompare> dc;
            std::vector<DeferredCompare> dc_next;
            dc.emplace_back(a.mesh->getSurrogateTree().getNode(0), b.mesh->getSurrogateTree().getNode(0));
            std::vector<int> num_hits;
            std::vector<Contact<real>> hits;
            while (!dc.empty()) {
                hits.clear();
                // findContactsBucket<branching, real>(dc, a, b, hits, num_hits);
                bucket_comparison.findContactsBucketCurrent(dc, a, b, hits, num_hits);
                int last_hit_num = 0;
                for (int i = 0; i < dc.size(); i++) {
                    if (num_hits[i] > last_hit_num) {
                        bool advance_a = false;
                        bool advance_b = false;
                        if (!dc[i].a.is_inner && !dc[i].b.is_inner) {
                            result.insert(result.end(), hits.begin() + last_hit_num, hits.begin() + num_hits[i]);
                        }
                        else if (!dc[i].a.is_inner) {
                            advance_b = true;
                        }
                        else if (!dc[i].b.is_inner) {
                            advance_a = true;
                        }
                        else {
                            // TODO this decision of which to advance first is a little arbitrary.
                            //    Would probably be better to advance the one with the largest area.
                            if (dc[i].a.depth < dc[i].b.depth) {
                                advance_a = true;
                            }
                            else {
                                advance_b = true;
                            }
                        }
                        if (advance_a) {
                            for (int c = 0; c < dc[i].a.num_children; c++) {
                                dc_next.emplace_back(a.mesh->getSurrogateTree().getNode(c + dc[i].a.child_id_first), dc[i].b);
                            }
                        }
                        else if (advance_b) {
                            for (int c = 0; c < dc[i].b.num_children; c++) {
                                dc_next.emplace_back(dc[i].a, b.mesh->getSurrogateTree().getNode(c + dc[i].b.child_id_first));
                            }
                        }
                    }
                    last_hit_num = num_hits[i];
                }
                dc.swap(dc_next);
                dc_next.clear();
            }
            return result;
        };

        template<typename real, int branching, int bucket_size>
        std::vector<Contact<real>> compareTreesFullLast(Particle& a, Particle& b, strategy::ContactDetectionStrategy& bucket_comparison) {
            std::vector<Contact<real>> result;
            
            std::vector<DeferredCompare> dc;
            std::vector<DeferredCompare> dc_next;
            dc.emplace_back(a.mesh->getSurrogateTree().getNode(0), b.mesh->getSurrogateTree().getNode(0));
            std::vector<int> num_hits;
            std::vector<Contact<real>> hits;
            while (!dc.empty()) {
                hits.clear();
                // findContactsBucket<branching, real>(dc, a, b, hits, num_hits);
                bucket_comparison.findContactsBucketLast(dc, a, b, hits, num_hits);
                int last_hit_num = 0;
                for (int i = 0; i < dc.size(); i++) {
                    if (num_hits[i] > last_hit_num) {
                        bool advance_a = false;
                        bool advance_b = false;
                        if (!dc[i].a.is_inner && !dc[i].b.is_inner) {
                            result.insert(result.end(), hits.begin() + last_hit_num, hits.begin() + num_hits[i]);
                        }
                        else if (!dc[i].a.is_inner) {
                            advance_b = true;
                        }
                        else if (!dc[i].b.is_inner) {
                            advance_a = true;
                        }
                        else {
                            // TODO this decision of which to advance first is a little arbitrary.
                            //    Would probably be better to advance the one with the largest area.
                            if (dc[i].a.depth < dc[i].b.depth) {
                                advance_a = true;
                            }
                            else {
                                advance_b = true;
                            }
                        }
                        if (advance_a) {
                            for (int c = 0; c < dc[i].a.num_children; c++) {
                                dc_next.emplace_back(a.mesh->getSurrogateTree().getNode(c + dc[i].a.child_id_first), dc[i].b);
                            }
                        }
                        else if (advance_b) {
                            for (int c = 0; c < dc[i].b.num_children; c++) {
                                dc_next.emplace_back(dc[i].a, b.mesh->getSurrogateTree().getNode(c + dc[i].b.child_id_first));
                            }
                        }
                    }
                    last_hit_num = num_hits[i];
                }
                dc.swap(dc_next);
                dc_next.clear();
            }
            return result;
        };

        template<typename real, int branching, int bucket_size>
        std::vector<ContinuousContact<real>> compareTreesFullContinuous(Particle& a, Particle& b, real max_time, real& min_time_out) {
            // __itt_string_handle* bucket_contact_soup_task = __itt_string_handle_create("findContactsBucketContinuousComparison");
            // __itt_string_handle* bucket_contact_connected_task = __itt_string_handle_create("findContactsBucketConnectedContinuousComparison");
            
            std::vector<ContinuousContact<real>> result;
            
            min_time_out = max_time;
            a.projectFutureState(max_time);
            b.projectFutureState(max_time);

            //std::vector<DeferredCompare> dc;
            //std::vector<DeferredCompare> dc_next;

            std::vector<DeferredCompare> dc_soup;
            std::vector<DeferredCompare> dc_connected;
            std::vector<DeferredCompare> dc_soup_next;
            std::vector<DeferredCompare> dc_connected_next;
            dc_soup.emplace_back(a.mesh->getSurrogateTree().getNode(0), b.mesh->getSurrogateTree().getNode(0));
            std::vector<int> num_hits;
            std::vector<ContinuousContact<real>> hits;
            while (!dc_soup.empty() || !dc_connected.empty()) {
                hits.clear();
                num_hits.clear();

                // __itt_task_begin(domain, __itt_null, __itt_null, bucket_contact_soup_task);
                findContactsBucketContinuousComparison<branching, real>(dc_soup, a, b, hits, num_hits);
                // __itt_task_end(domain);

                int last_hit_num = 0;
                for (int i = 0; i < dc_soup.size(); i++) {
                    if (num_hits[i] > last_hit_num) {
                        bool advance_a = false;
                        bool advance_b = false;
                        if (!dc_soup[i].a.is_inner && !dc_soup[i].b.is_inner) {
                            result.insert(result.end(), hits.begin() + last_hit_num, hits.begin() + num_hits[i]);
                            for (ContinuousContact<real>& c : hits) {
                                min_time_out = std::min(min_time_out, c.toc * max_time);
                            }
                        }
                        else if (!dc_soup[i].a.is_inner) {
                            advance_b = true;
                        }
                        else if (!dc_soup[i].b.is_inner) {
                            advance_a = true;
                        }
                        else {
                            // TODO this decision of which to advance first is a little arbitrary.
                            //    Would probably be better to advance the one with the largest area.
                            if (dc_soup[i].a.depth < dc_soup[i].b.depth) {
                                advance_a = true;
                            }
                            else {
                                advance_b = true;
                            }
                        }
                        if (advance_a) {
                            for (int c = 0; c < dc_soup[i].a.num_children; c++) {
                                const model::Node& a_node = a.mesh->getSurrogateTree().getNode(c + dc_soup[i].a.child_id_first);
                                const model::Node& b_node = dc_soup[i].b;
                                if (a_node.is_inner || b_node.is_inner) {
                                    dc_soup_next.emplace_back(a_node, b_node);
                                } else {
                                    dc_connected_next.emplace_back(a_node, b_node);
                                }
                            }
                        }
                        else if (advance_b) {
                            for (int c = 0; c < dc_soup[i].b.num_children; c++) {
                                const model::Node& a_node = dc_soup[i].a;
                                const model::Node& b_node = b.mesh->getSurrogateTree().getNode(c + dc_soup[i].b.child_id_first);
                                if (a_node.is_inner || b_node.is_inner) {
                                    dc_soup_next.emplace_back(a_node, b_node);
                                } else {
                                    dc_connected_next.emplace_back(a_node, b_node);
                                }
                            }
                        }
                    }
                    last_hit_num = num_hits[i];
                }
                

                hits.clear();
                num_hits.clear();
                
                // __itt_task_begin(domain, __itt_null, __itt_null, bucket_contact_connected_task);
                findContactsBucketConnectedContinuousComparison<branching, real>(dc_connected, a, b, hits, num_hits);
                // __itt_task_end(domain);

                last_hit_num = 0;
                for (int i = 0; i < dc_connected.size(); i++) {
                    if (num_hits[i] > last_hit_num) {
                        bool advance_a = false;
                        bool advance_b = false;
                        if (!dc_connected[i].a.is_inner && !dc_connected[i].b.is_inner) {
                            result.insert(result.end(), hits.begin() + last_hit_num, hits.begin() + num_hits[i]);
                            for (ContinuousContact<real>& c : hits) {
                                min_time_out = std::min(min_time_out, c.toc * max_time);
                            }
                        }
                        else if (!dc_connected[i].a.is_inner) {
                            advance_b = true;
                        }
                        else if (!dc_connected[i].b.is_inner) {
                            advance_a = true;
                        }
                        else {
                            // TODO this decision of which to advance first is a little arbitrary.
                            //    Would probably be better to advance the one with the largest area.
                            if (dc_connected[i].a.depth < dc_connected[i].b.depth) {
                                advance_a = true;
                            }
                            else {
                                advance_b = true;
                            }
                        }
                        if (advance_a) {
                            for (int c = 0; c < dc_connected[i].a.num_children; c++) {
                                const model::Node& a_node = a.mesh->getSurrogateTree().getNode(c + dc_connected[i].a.child_id_first);
                                const model::Node& b_node = dc_connected[i].b;
                                if (a_node.is_inner || b_node.is_inner) {
                                    dc_soup_next.emplace_back(a_node, b_node);
                                } else {
                                    dc_connected_next.emplace_back(a_node, b_node);
                                }
                            }
                        }
                        else if (advance_b) {
                            for (int c = 0; c < dc_connected[i].b.num_children; c++) {
                                const model::Node& a_node = dc_connected[i].a;
                                const model::Node& b_node = b.mesh->getSurrogateTree().getNode(c + dc_connected[i].b.child_id_first);
                                if (a_node.is_inner || b_node.is_inner) {
                                    dc_soup_next.emplace_back(a_node, b_node);
                                } else {
                                    dc_connected_next.emplace_back(a_node, b_node);
                                }
                            }
                        }
                    }
                    last_hit_num = num_hits[i];
                }

                dc_soup.swap(dc_soup_next);
                dc_soup_next.clear();
                dc_connected.swap(dc_connected_next);
                dc_connected_next.clear();
            }
            return result;
        };

        std::vector<ContinuousContact<double>> compareTreesFullContinuousWitnesses(Particle& a, Particle& b, ContactStateCache& cache, double max_time, double& min_time_out);
    }
}
