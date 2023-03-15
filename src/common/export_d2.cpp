#include "export_d2.h"

#include "../globals.h"

#include <iostream>
#include <fstream>
using namespace std;

Delta2::D2Writer::D2Writer(std::vector<Delta2::Particle>& particles) {
    for (Delta2::Particle& p: particles) {
        states.push_back({});
    }
}

void Delta2::D2Writer::capture(std::vector<Delta2::Particle>& particles) {
    for (int p_i = 0; p_i < particles.size(); p_i++) {
        Delta2::Particle& p = particles[p_i];

        if (states[p_i].size() > 0) {
            if (p.last_state.getTime() > states[p_i][states[p_i].size()-1].getTime()) {
                states[p_i].push_back(p.last_state);
            }
        } else {
            states[p_i].push_back(p.last_state);
        }
    }
}

void Delta2::D2Writer::write(std::vector<Delta2::Particle>& particles, int frame_rate, std::string file_name) {
    Delta2::globals::logger.printf(1, "Writing d2 file\n");
    for (int p_i = 0; p_i < particles.size(); p_i++) {
        Delta2::Particle& p = particles[p_i];
        states[p_i].push_back(p.current_state);
    }

    std::vector<std::shared_ptr<Delta2::MeshData>> mesh_types;
    // Collect together a list of all the mesh objects used
    for (Delta2::Particle& p: particles) {
        if (std::find(mesh_types.begin(), mesh_types.end(), p.mesh) == mesh_types.end()) {
            mesh_types.push_back(p.mesh);
        }
    }

    ofstream file;
    file.open(file_name);
    if (!file.is_open()) {
        Delta2::globals::logger.printf(1, "Error opening file\n");
        return;
    }

    file << "Meshes:\n";
    for (int m_i = 0; m_i < mesh_types.size(); m_i++) {
        std::shared_ptr<Delta2::MeshData> m = mesh_types[m_i];
        file << "m " << m_i << "\n";
        file << m->serialise();
    }
    file << "\nParticles:\n";
    // For each particle write it's mesh reference 
    for (int p_i = 0; p_i < particles.size(); p_i++) {
        Delta2::Particle& p = particles[p_i];
        int mesh_ref = -1;
        for (int m_i = 0; m_i < mesh_types.size(); m_i++) {
            std::shared_ptr<Delta2::MeshData> m = mesh_types[m_i];
            if (m == p.mesh) {
                mesh_ref = m_i;
                break;
            }
        }
        file << "pm " << p_i << " ";
        file << mesh_ref << "\n";
    }

    Delta2::globals::logger.printf(2, "Writing frames\n");

    float time = 0.0;
    bool reached_end = false;
    std::vector<int> current_capture;  // points to the stored state one after the current export time
    for (Delta2::Particle& p: particles) {
        current_capture.push_back(1);
    }
    while (!reached_end) {
        for (int p_i = 0; p_i < particles.size(); p_i++) {
            if (!particles[p_i].is_static) {
                while (states[p_i][current_capture[p_i]].getTime() <= time) {
                    current_capture[p_i]++;
                    if (current_capture[p_i] >= states[p_i].size()) {
                        reached_end = true;
                        break;
                    }
                }
            }
        }

        if (!reached_end) {
            file << "t " << time << "\n";
            for (int p_i = 0; p_i < particles.size(); p_i++) {
                file << "p " << p_i << "\n";
                file << states[p_i][current_capture[p_i]].interpolate(states[p_i][current_capture[p_i]-1], time).serialise();
                file << "\n";
            }
        }

        time += 1.0 / (double)frame_rate;
    }

    file.flush();
    file.close();
}