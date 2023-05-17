#include "full_tree_comparison.h"
#include "../globals.h"

using namespace Delta2;

namespace Delta2 {
    namespace collision {
        std::vector<collision::ContinuousContact<double>> compareTreesFullContinuousWitnesses(Particle& a, Particle& b, collision::ContactStateCache& cache, double max_time, double& min_time_out) {
            const int branching = 8;
            
            // __itt_string_handle* bucket_contact_soup_task = __itt_string_handle_create("findContactsBucketContinuousComparison");
            // __itt_string_handle* bucket_contact_connected_task = __itt_string_handle_create("findContactsBucketConnectedContinuousComparison");
            
            std::vector<ContinuousContact<double>> result;
            
            min_time_out = max_time;
            a.projectFutureState(max_time);
            b.projectFutureState(max_time);

            //std::vector<DeferredCompare> dc;
            //std::vector<DeferredCompare> dc_next;

            collision::ContactState& contact_state = cache.getContactState(a, b);

            // Check witnesses and populate DeferredCompare
            // Expand all the DeferredCompares
            // Create witnesses where non-contacts were identified
            // Append noWitnesses with contact

            collision::ContactState next_contact_state;

            for (collision::SeparationWitness& w : contact_state.witnesses) {
                bool separated = w.check(a, b);
                if (separated) {
                    next_contact_state.witnesses.push_back(w);
                    // TODO Needs a heuristic here for when to step back up the tree and a way to implement stepping up all children together.
                }
                else {
                    contact_state.noWitnesses.push_back(w.toDeferredCompare(a, b));
                }
            }

            std::vector<DeferredCompare> dc_soup;
            std::vector<DeferredCompare> dc_connected;
            std::vector<DeferredCompare> dc_soup_next;
            std::vector<DeferredCompare> dc_connected_next;

            for (DeferredCompare& dc : contact_state.noWitnesses) {
                std::shared_ptr<model::Bucket> a_bucket = a.mesh->getSurrogateTree().getBucket(dc.a.bucket_id);
                std::shared_ptr<model::Bucket> b_bucket = b.mesh->getSurrogateTree().getBucket(dc.b.bucket_id);

                if (a_bucket->getConfig() == model::BucketConfig::CONNECTED && b_bucket->getConfig() == model::BucketConfig::CONNECTED) {
                    dc_connected.push_back(dc);
                }
                else {
                    dc_soup.push_back(dc);
                }
            }

            std::vector<int> num_hits;
            std::vector<ContinuousContact<double>> hits;
            while (!dc_soup.empty() || !dc_connected.empty()) {
                hits.clear();
                num_hits.clear();

                // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, bucket_contact_soup_task);
                findContactsBucketContinuousComparison<branching, double>(dc_soup, a, b, hits, min_time_out / max_time, num_hits);
                // __itt_task_end(globals::itt_handles.detailed_domain);

                int last_hit_num = 0;
                for (int i = 0; i < dc_soup.size(); i++) {
                    if (num_hits[i] > last_hit_num) {
                        bool advance_a = false;
                        bool advance_b = false;
                        if (!dc_soup[i].a.is_inner && !dc_soup[i].b.is_inner) {
                            result.insert(result.end(), hits.begin() + last_hit_num, hits.begin() + num_hits[i]);
                            for (ContinuousContact<double>& c : hits) {
                                min_time_out = std::min(min_time_out, c.toc * max_time);
                            }
                            next_contact_state.noWitnesses.push_back(dc_soup[i]);
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
                    else {
                        bool success;
                        next_contact_state.attemptWitness(a, b, dc_soup[i]);
                    }
                    last_hit_num = num_hits[i];
                }
                
                hits.clear();
                num_hits.clear();
                
                // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, bucket_contact_connected_task);
                findContactsBucketConnectedContinuousComparison<branching, double>(dc_connected, a, b, hits, min_time_out / max_time, num_hits);
                // __itt_task_end(globals::itt_handles.detailed_domain);

                last_hit_num = 0;
                for (int i = 0; i < dc_connected.size(); i++) {
                    if (num_hits[i] > last_hit_num) {
                        bool advance_a = false;
                        bool advance_b = false;
                        if (!dc_connected[i].a.is_inner && !dc_connected[i].b.is_inner) {
                            result.insert(result.end(), hits.begin() + last_hit_num, hits.begin() + num_hits[i]);
                            for (ContinuousContact<double>& c : hits) {
                                min_time_out = std::min(min_time_out, c.toc * max_time);
                            }
                            next_contact_state.noWitnesses.push_back(dc_connected[i]);
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
                    else {
                        next_contact_state.attemptWitness(a, b, dc_connected[i]);
                    }
                    last_hit_num = num_hits[i];
                }

                dc_soup.swap(dc_soup_next);
                dc_soup_next.clear();
                dc_soup.swap(dc_soup_next);
                dc_soup_next.clear();
            }
            return result;
        }
    }
}
