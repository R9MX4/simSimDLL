import os,time

def count_lines_in_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return len(file.readlines())

current_directory = os.getcwd()
file_paths = [os.path.join(current_directory, f) for f in os.listdir(current_directory) if os.path.isfile(os.path.join(current_directory, f))]

file_cnt = 0
line_cnt = 0
for file_path in file_paths:
    if not file_path.endswith('.h') and not file_path.endswith('.cpp'): continue
    line_count = count_lines_in_file(file_path)
    file_cnt += 1
    line_cnt += line_count
    print(f"File: %-50s\tLines: {line_count}" % file_path)

print(f"Total {file_cnt} Files, {line_cnt} lines")
time.sleep(10)