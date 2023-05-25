import re
import os

rows = [8, 16, 24, 32, 40, 48, 56, 64]
threads = [1, 2, 4, 8, 16, 32, 48, 64]

path = "/nobackup/wxmv95/thread_scaling/round_tower/stdout"
prefix = "stdout"
extension = "log"

files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and f.split('.')[-1] == extension and f.split('.')[0].startswith(prefix)]

files.sort(key=lambda x: x.split('.')[1])
most_recent = files[-(len(rows) * len(threads)):]

print(most_recent)

data = {t: {r: () for r in rows} for t in threads}
print(data)
for file in most_recent:
    print(f"Opening file {file}")
    opened = open(os.path.join(path, file), 'r')
    lines = opened.readlines()
    opened.close()
    check_rows = int(lines[0].split(" ")[-1])
    check_towers = int(lines[1].split(" ")[-1])
    if len(lines) > 3:
        timing = lines[3].split()
        user = float(timing[0][:-1])
        system = float(timing[1][:-1])
        elapsed = 60 * float(timing[2].split(":")[0]) + float(timing[2].split(":")[1])
        percent = float(timing[3][:-1])
        page_faults = int(timing[6].split("pf")[0])
    else:
        timing = None
        user = None
        system = None
        elapsed = None
        percent = None
        page_faults = None
    print(check_towers, check_rows)
    data[check_towers][check_rows] = ((check_rows, check_towers, user, system, elapsed, percent, page_faults))

print([[data[t][r] for r in rows] for t in threads])
