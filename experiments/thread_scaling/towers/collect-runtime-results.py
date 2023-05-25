import re
import os

pairs = [1, 2, 4, 8, 12]
threads = [1, 2, 4, 5, 16, 17, 32]

path = "/nobackup/wxmv95/thread_scaling/towers/stdout"
prefix = "stdout"
extension = "log"

files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and f.split('.')[-1] == extension and f.split('.')[0].startswith(prefix)]

files.sort(key=lambda x: x.split('.')[1])
most_recent = files[-(len(pairs) * len(threads)):]

print(most_recent)

data = {t: {r: () for r in pairs} for t in threads}
print(data)
for file in most_recent:
    print(f"Opening file {file}")
    opened = open(os.path.join(path, file), 'r')
    lines = opened.readlines()
    opened.close()
    check_pairs = int(lines[0].split(" ")[-1])
    check_threads = int(lines[1].split(" ")[-1])
    timing = lines[2].split()
    user = float(timing[0][:-1])
    system = float(timing[1][:-1])
    elapsed = 60 * float(timing[2].split(":")[0]) + float(timing[2].split(":")[1])
    percent = float(timing[3][:-1])
    page_faults = int(timing[6].split("pf")[0])
    print(check_threads, check_pairs)
    data[check_threads][check_pairs] = ((check_pairs, check_threads, user, system, elapsed, percent, page_faults))

print([[data[t][r] for r in pairs] for t in threads])
