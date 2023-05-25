import re
import os

towers = [1,2,4,6,8,12,16]

path = "/nobackup/wxmv95/problem_size_scaling/hoppers/stdout"
prefix = "stdout"
extension = "log"

files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and f.split('.')[-1] == extension and f.split('.')[0].startswith(prefix)]

files.sort(key=lambda x: x.split('.')[1])
most_recent = files[-len(towers):]

print(most_recent)

data = {t: () for t in towers}
print(data)
for file in most_recent:
    print(f"Opening file {file}")
    opened = open(os.path.join(path, file), 'r')
    lines = opened.readlines()
    opened.close()
    check_towers = int(lines[0].split(" ")[-1])
    timing = lines[1].split()
    user = float(timing[0][:-1])
    system = float(timing[1][:-1])
    elapsed = 60 * float(timing[2].split(":")[0]) + float(timing[2].split(":")[1])
    percent = float(timing[3][:-1])
    page_faults = int(timing[6].split("pf")[0])
    data[check_towers] = ((check_towers, user, system, elapsed, percent, page_faults))

print([data[t] for t in towers])
