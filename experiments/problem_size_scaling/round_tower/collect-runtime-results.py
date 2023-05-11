import re
import os

rows = [8, 16, 24, 32, 40, 48, 56, 64]
towers = [2, 4, 6, 8, 10, 12]

path = "/nobackup/wxmv95/problem_size_scaling/round_tower/stdout"
prefix = "stdout"
extension = "log"

files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and f.split('.')[-1] == extension and f.split('.')[0].startswith(prefix)]

files.sort(key=lambda x: x.split('.')[1])
most_recent = files[-(len(rows) * len(towers)):]

print(most_recent)

data = {t: {r: () for r in rows} for t in towers}
print(data)
for file in most_recent:
    print(f"Opening file {file}")
    opened = open(os.path.join(path, file), 'r')
    lines = opened.readlines()
    opened.close()
    check_rows = int(lines[0].split(" ")[-1])
    check_towers = int(lines[1].split(" ")[-1])
    timing = lines[2].split()
    user = float(timing[0][:-1])
    system = float(timing[1][:-1])
    elapsed = 60 * float(timing[2].split(":")[0]) + float(timing[2].split(":")[1])
    percent = float(timing[3][:-1])
    page_faults = int(timing[6].split("pf")[0])
    print(check_towers, check_rows)
    data[check_towers][check_rows] = ((check_rows, check_towers, user, system, elapsed, percent, page_faults))

print([[data[t][r] for r in rows] for t in towers])
