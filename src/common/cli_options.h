#pragma once

#include "pcg_basic.h"
#include <string>
#include <stdexcept>
#include "CLI11.h"

namespace Delta2 {
	namespace common {
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
			std::string mesh_metric;
			float geo_eps;
			bool print_tree;
			bool print_debug;
			bool gui;

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
				scenario_lod = -1;  // -1 for full
				mesh_metric = "";
				geo_eps = 0.01;
				print_tree = false;
				print_debug = false;
				gui = false;

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

				app.add_option("-s,--scenario", scenario, "Scenario (Not implemented)");
				app.add_option("-t,--time_step", time_step_size, "Time step size");
				app.add_option("-a,--adaptive_time_step_factor", adaptive_time_step_factor, "How small compared to the target time step the adaptive step sizes are allowed to go");
				app.add_option("-n,--num_steps", num_time_steps, "Number of time steps");
				app.add_option("-f,--final_time", final_time, "Simulate up until this time");
				app.add_flag("-e,--export", export_result, "Export to VTK (Not implemented)");
				app.add_flag("-p,--print_tree", print_tree, "Print the surrogate tree as a Python compatible string (Not implememted)");
				app.add_flag("-d,--print_debug", print_debug, "Print iteration debug information");
				app.add_flag("-i,--print_iteration_count", print_iteration_count, "Print the number of iterations at each depth");
				app.add_flag("--suppress_all_output", suppress_all_output, "Prevent anything being printed to the console each iteration");
				app.add_option("--mesh_metrics", mesh_metric, "Takes a path to a mesh and computes some metrics about it");
				app.add_option("--geo_eps", geo_eps, "Geometry eps");
				app.add_flag("-g,--gui", gui, "Toggle GUI on");

				CLI11_PARSE(app, argc, argv);

				print();
				validate();

				return 0;
			}
		};
	}
}
