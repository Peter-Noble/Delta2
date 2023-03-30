import bpy

import re
from dataclasses import dataclass
from typing import List

wm = bpy.context.window_manager

path = "/home/peter/Documents/particles/export/export.d2"

f = open(path, "r")
lines = [line.strip() for line in f.readlines()]

wm.progress_begin(0, 2)

i = 0

assert(lines[i].strip() == "Meshes:")
i += 1


@dataclass
class Mesh:
    verts: List
    faces: List
    
    
@dataclass
class State:
    time: float
    rotation: [float, float, float, float]
    angular: [float, float, float]
    translation: [float, float, float]
    velocity: [float, float, float]
    
@dataclass
class Particle:
    mesh: int
    states: List
    
meshes = []

all_meshes = False
m_i = 0
while not all_meshes:
    expected = "m {}$".format(m_i)
    r = re.match(expected, lines[i])
    if r:
        verts = []
        faces = []
        i += 1
        
        assert(lines[i] == "V")
        i += 1
        while True:
            print(repr(lines[i]))
            vr = re.match("([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?) *([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?) *([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?)$", lines[i])
            print(vr)
            if vr:
                print("Match v")
                i += 1
                verts.append([float(vr.group(1)), float(vr.group(2)), float(vr.group(3))])
            else:
                break
        
        print("Pre F check", repr(lines[i]))
        assert(lines[i] == "F")
        i += 1
        while True:
            vf = re.match("([0-9]+) *([0-9]+) *([0-9]+) *$", lines[i])
            print("vf", vf)
            if vf:
                i += 1
                faces.append([int(vf.group(1)), int(vf.group(2)), int(vf.group(3))])
            else:
                print("Leaving face search for:", lines[i])
                break
        
        meshes.append(Mesh(verts, faces))
        m_i += 1
    else:
        all_meshes = True
    wm.progress_update(i / len(lines))

assert(len(meshes) > 0)
i += 1
print(lines[i])
assert(lines[i] == "Particles:")
i += 1

particles = []

p_i = 0
while True:
    rpm = re.match("pm {} ([0-9]+)$".format(p_i), lines[i])
    if rpm:
        i += 1
        particles.append(Particle(int(rpm.group(1)), []))
        
        p_i += 1
    else:
        break

wm.progress_update(i / len(lines))

max_time = 0;

while i < len(lines):
    print("i", i, lines[i])
    rt = re.match("t ([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?)$", lines[i])  # target time.  The actual time on states might differ depending on how exporter ends up being implemented
    max_time = max(max_time, float(rt.group(1)))
    print("rt", rt)
    if rt:
        i += 1
        for p_i in range(len(particles)):
            print(i, "p: ", lines[i])
            assert(re.match("p {}$".format(p_i), lines[i]))
            i += 1
            
            rtimep = re.match("([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?)$", lines[i])
            assert(rtimep)
            time = float(rtimep.group(1))
            i += 1

            print(lines[i])
            rrotation = re.match("([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?) ([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?) ([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?) ([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?)$", lines[i])
            assert(rrotation)
            rotation = [float(rrotation.group(1)), float(rrotation.group(2)), float(rrotation.group(3)), float(rrotation.group(4))]
            i += 1
            
            # print(lines[i])
            # rangular = re.match("([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?) ([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?) ([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?)$", lines[i])
            # print("rangular", rangular)
            # assert(rangular)
            # angular = [float(rangular.group(1)), float(rangular.group(2)), float(rangular.group(3))]
            # i += 1
            angular = [0, 0, 0]
            
            rtranslation = re.match("([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?) ([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?) ([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?)$", lines[i])
            assert(rtranslation)
            translation = [float(rtranslation.group(1)), float(rtranslation.group(2)), float(rtranslation.group(3))]
            i += 1
            
            # rvelocity = re.match("([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?) ([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?) ([ |\-]?[0-9]+(?:\.[0-9]+)?(?:e[-|+]?[0-9]+)?)$", lines[i])
            # assert(rvelocity)
            # velocity = [float(rvelocity.group(1)), float(rvelocity.group(2)), float(rvelocity.group(3))]
            # i += 1
            velocity = [0, 0, 0]
            
            particles[p_i].states.append(State(time, rotation, angular, translation, velocity))
            p_i += 1
            i += 1
            wm.progress_update(i / len(lines))
    else:
        break

print("Meshes {} particles {}".format(len(meshes), len(particles)))

b_meshes = []
for i, m in enumerate(meshes):
    b_meshes.append(bpy.data.meshes.new('mesh_{}'.format(i)))
    b_meshes[-1].from_pydata(m.verts, [], m.faces)


b_particles = []
for i, p in enumerate(particles):
    b_particles.append(bpy.data.objects.new("Object {}".format(i), b_meshes[p.mesh]))
    bpy.context.collection.objects.link(b_particles[-1])
    
    for s in p.states:
        b_particles[-1].location = s.translation
        """b_particles[-1].rotation_mode = 'ZYX'"""
        b_particles[-1].rotation_mode = 'QUATERNION'
        b_particles[-1].rotation_quaternion = s.rotation
        frame = round(s.time * bpy.data.scenes["Scene"].render.fps) + 1
        b_particles[-1].keyframe_insert(data_path="location", frame=frame)
        b_particles[-1].keyframe_insert(data_path="rotation_euler", frame=frame)
    
    wm.progress_update(1 + i/len(particles))

wm.progress_end()
