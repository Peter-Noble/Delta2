#include "export_d2.h"

#include <iostream>
#include <fstream>
using namespace std;

Delta2::D2Writer::D2Writer(std::vector<Delta2::Particle>& particles) {
    
}

void Delta2::D2Writer::capture(std::vector<Delta2::Particle>& particles) {
    for (int p_i = 0; p_i < particles.size(); p_i++) {
        Delta2::Particle& p = particles[p_i];
        states[p_i].push_back(p.last_state);
    }
}

void Delta2::D2Writer::write(std::vector<Delta2::Particle>& particles, int frame_rate, std::string file_name) {
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
    file << "Meshes:\n";
    for (int m_i = 0; m_i < mesh_types.size(); m_i++) {
        std::shared_ptr<Delta2::MeshData> m = mesh_types[m_i];
        file << "m " << m_i << "\n";
        file << m->serialise();
    }
    file << "Particles:\n";
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

    float time = 0.0;
    bool reached_end = false;
    std::vector<int> current_capture;  // points to the stored state one after the current export time
    for (Delta2::Particle& p: particles) {
        current_capture.push_back(0);
    }
    while (!reached_end) {
        file << "t " << time << "\n";
        for (int p_i = 0; p_i < particles.size(); p_i++) {
            while (states[p_i][current_capture[p_i]].getTime() <= time) {
                current_capture[p_i]++;
                if (current_capture[p_i] >= states[p_i].size()) {
                    reached_end = true;
                    break;
                }
            }
        }

        if (!reached_end) {
            for (int p_i = 0; p_i < particles.size(); p_i++) {
                file << "p " << p_i << "\n";
                file << states[p_i][current_capture[p_i]].interpolate(states[p_i][current_capture[p_i]-1], time).serialise();
                file << "\n";
            }
        }

        time += 1.0 / (double)frame_rate;
    }
}
