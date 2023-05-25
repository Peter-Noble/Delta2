import re
import os

pairs = [64, 128, 512, 2048, 8192, 16384]
threads = [1, 2, 4, 5, 16, 17, 32, 33, 64, 65, 128]

path = "/nobackup/wxmv95/thread_scaling/pairs/stdout"
prefix = "stdout"
extension = "log"

files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and f.split('.')[-1] == extension and f.split('.')[0].startswith(prefix)]

files.sort(key=lambda x: x.split('.')[1])
most_recent = files[-(len(pairs) * len(threads)):]

print(files)
print(len(files))
print(most_recent)
print(len(most_recent))

data = {t: {r: () for r in pairs} for t in threads}
print(data)
for file in most_recent:
    print(f"Opening file {file}")
    opened = open(os.path.join(path, file), 'r')
    lines = opened.readlines()
    opened.close()
    check_pairs = int(lines[0].split(" ")[-1])
    check_threads = int(lines[1].split(" ")[-1])
    if len(lines) > 2:
        timing = lines[2].split()
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
    print(check_threads, check_pairs)
    data[check_threads][check_pairs] = ((check_pairs, check_threads, user, system, elapsed, percent, page_faults))

print([[data[t][r] for r in pairs] for t in threads])
