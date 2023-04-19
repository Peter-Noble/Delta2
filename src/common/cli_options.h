#pragma once

#include "pcg_basic.h"
#include <string>
#include <stdexcept>
#include "CLI11.h"

namespace Delta2 {
	namespace common {
		const int PRINT_STARTUP = 0;
		const int PRINT_OUTER_STEPS = 1;
		const int PRINT_DEBUG = 2;
		
		struct Options {
			float time_step_size;
			float adaptive_time_step_factor;
			int num_time_steps;
			double final_time;
			bool export_result;
			int scenario;
			bool print_iteration_count;
			bool suppress_all_output;
			int export_skip;
			int scenario_lod;
			int scenario_size;
			std::string mesh_metric;
			float geo_eps;
			bool print_tree;
			bool print_debug;
			bool gui;
			bool view_contacts;
			int threads;
			int print_priority;
			bool limit_clusters_per_step;
			bool check_post_solve;
			bool local_ts;

			int sequential_impulse_total_iterations;
			int sequential_impulse_inner_iterations;
			int sequential_impulse_grain_size;

			int sequential_parallel_threshold;
			int sequential_parallel_individual_colour_threshold;
			int sequential_parallel_grain_size;

			float cluster_separation_factor;

			pcg32_random_t rng;

			Options() {
				time_step_size = 0.01;
				adaptive_time_step_factor = 0.001;
				num_time_steps = 100;
				final_time = -1.0;
				export_result = false;
				scenario = 0;
				print_iteration_count = false;
				suppress_all_output = false;
				export_skip = 0;
				scenario_lod = 0;
				scenario_size = -1;
				mesh_metric = "";
				geo_eps = 0.01;
				print_tree = false;
				print_debug = false;
				gui = false;
				view_contacts = false;
				threads = -1;
				limit_clusters_per_step = false;
				check_post_solve = false;
				local_ts = false;

				print_priority = 0;

				sequential_impulse_total_iterations = 1000;
				sequential_impulse_inner_iterations = 2;
				sequential_impulse_grain_size = 1000;

				sequential_parallel_threshold = 1000;
				sequential_parallel_individual_colour_threshold = 200;
				sequential_parallel_grain_size = 100;

				cluster_separation_factor = 0.75;

				pcg32_srandom_r(&rng, 42u, 54u);
			}
			int rand(int range) {
				return (int)pcg32_boundedrand_r(&rng, range);
			}
			float rand_float(float range) {
				float M = 1073741824;
				return range * (float)pcg32_boundedrand_r(&rng, M) / M;
			}
			void seed_rand(int state, int start) {
				pcg32_srandom_r(&rng, state, start);
			}
			void print() {
				printf("================ Mode ================\n");
				printf("Mesh metric:           %s\n", mesh_metric.c_str());
				printf("================ Settings ================\n");
				printf("Scenario:              %i\n", scenario);
				printf("Scenario LOD:          %i\n", scenario_lod);
				printf("Export skip:           %i\n", export_skip);
				printf("Geo eps:               %f\n", geo_eps);
				if (print_tree) {
					printf("Print tree:            true\n");
				}
				else {
					printf("Print tree:            false\n");
				}
				if (print_debug) {
					printf("Print debug:           true\n");
				}
				else {
					printf("Print debug:           false\n");
				}
				if (gui) {
					printf("GUI:                   true\n");
				}
				else {
					printf("GUI:                   false\n");
				}
				if (view_contacts) {
					printf("View contacts:         true\n");
				}
				else {
					printf("View contacts:         false\n");
				}
				if (limit_clusters_per_step) {
					printf("Limit clusters /step   true\n");
				}
				else {
					printf("Limit cluster /step    false\n");
				}
				if (check_post_solve) {
					printf("Check post solve       true\n");
				}
				else {
					printf("Check post solve       false\n");
				}
				if (local_ts) {
					printf("Local time step        true\n");
				}
				else {
					printf("Local time step        false\n");
				}
				printf("Threads:               %i\n", threads);
				printf("Print priority:        %i\n", print_priority);
				printf("======== Time stepping ========\n");
				printf("Num steps:             %i\n", num_time_steps);
				printf("Final time:            %f\n", final_time);
				printf("Time step size:        %f\n", time_step_size);
				printf("Adaptive time factor:  %f\n", adaptive_time_step_factor);
				if (export_result) {
					printf("Export result:         true\n");
				}
				else {
					printf("Export result:         false\n");
				}
				if (print_iteration_count) {
					printf("Print iteration count: true\n");
				}
				else {
					printf("Print iteration count: false\n");
				}
				if (suppress_all_output) {
					printf("Suppress all output:   true\n");
				}
				else {
					printf("Suppress all output:   false\n");
				}
				printf("======== Sequential Impulses ========\n");
				printf("Seq total iterations:  %i\n", sequential_impulse_total_iterations);
				printf("Seq inner iterations:  %i\n", sequential_impulse_inner_iterations);
				printf("Seq grain size:        %i\n", sequential_impulse_grain_size);
				printf("Seq par threshold:     %i\n", sequential_parallel_threshold);
				printf("Seq par indiv thresh:  %i\n", sequential_parallel_individual_colour_threshold);
				printf("Seq par grain size:    %i\n", sequential_parallel_grain_size);

				printf("==========================================\n");
				printf("Cluster sep factor:    %f\n", cluster_separation_factor);
				printf("==========================================\n");
			}
			void validate() {
				if (scenario_lod < -1) {
					throw std::invalid_argument("Scenario lod must be a positive integer 0-> or -1");
				}
			}
			// Returns exit codes like main
			int fromArgs(int argc, char *argv[]) {
				CLI::App app{"App description"};

				app.add_option("-s,--scenario", scenario, "Scenario");
				app.add_option("-l,--scenario_lod", scenario_lod, "Scenario LOD");
				app.add_option("--scenario_size", scenario_size, "Scenario size");
				app.add_option("-t,--time_step", time_step_size, "Time step size");
				app.add_option("-a,--adaptive_time_step_factor", adaptive_time_step_factor, "How small compared to the target time step the adaptive step sizes are allowed to go");
				app.add_option("-n,--num_steps", num_time_steps, "Number of time steps");
				app.add_option("-f,--final_time", final_time, "Simulate up until this time");
				app.add_flag("-e,--export", export_result, "Export to d2 (custom) file");
				app.add_flag("-p,--print_tree", print_tree, "Print the surrogate tree as a Python compatible string (Not implememted)");
				app.add_flag("-d,--print_debug", print_debug, "Print iteration debug information");
				app.add_flag("-i,--print_iteration_count", print_iteration_count, "Print the number of iterations at each depth");
				app.add_flag("--suppress_all_output", suppress_all_output, "Prevent anything being printed to the console each iteration");
				app.add_option("--mesh_metrics", mesh_metric, "Takes a path to a mesh and computes some metrics about it");
				app.add_option("--geo_eps", geo_eps, "Geometry eps");
				app.add_flag("-g,--gui", gui, "Toggle GUI on");
				app.add_flag("--gui_view_contacts", view_contacts, "Display dots on the contacts between objects");
				app.add_flag("--limit_per_step", limit_clusters_per_step, "Limit the number of cluster handled per step to try and help load balancing");
				app.add_flag("--check_post_solve", check_post_solve, "Redo contact detection after contact solve to see if it was successful");
				app.add_flag("--local_ts", local_ts, "Local time stepping");
				app.add_option("--threads", threads, "Maximum threads");
				app.add_option("--print_priority", print_priority, "Print priority (0-Startup, 1-Outer steps, 2-Debug");

				app.add_option("-c,--cluster_sep_factor", cluster_separation_factor, "Multiply the min time-of-contact by this factor before separating clusters");

				app.add_option("--seq_total_iter", sequential_impulse_total_iterations, "Total number of iterations to use for sequential impulses");
				app.add_option("--seq_inner_iter", sequential_impulse_inner_iterations, "Number of inner iterations to use for sequential impulses");
				app.add_option("--seq_grain_size", sequential_impulse_grain_size, "Grain size to use for sequential impulses");
				
				app.add_option("--seq_par_thresh", sequential_parallel_threshold, "Number of hits in a cluster to activate graph colouring and parallel evaluation");
				app.add_option("--seq_par_indiv_thresh", sequential_parallel_individual_colour_threshold, "Minimum number of hits of one colour to evaluate in parallel");
				app.add_option("--seq_par_grain_size", sequential_parallel_grain_size, "Grain size to use for parallel hit evaluation");

				sequential_parallel_threshold = 200;
				sequential_parallel_individual_colour_threshold = 100;
				sequential_parallel_grain_size = 10;
				CLI11_PARSE(app, argc, argv);

				print();
				validate();

				return 0;
			}
		};
	}
}
