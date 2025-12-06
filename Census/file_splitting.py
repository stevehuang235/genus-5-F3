# File splitting 
import os
import sys

def split_and_print(file_path, lines_per_file):
    file_number = 1
    with open(file_path, 'r') as original_file:
        while True:
            lines = []
            for _ in range(lines_per_file):
                line = original_file.readline()
                if not line:
                    break
                lines.append(line)
            if not lines:
                break  # Exit the loop if there are no more lines to process

            output_file_path = f"./generic/data_unfiltered/generic_unfiltered_{file_number}.txt"
            with open(output_file_path, 'w') as output_file:
                output_file.writelines(lines)
            file_number += 1


# Replace 'input.txt' with the path to your text file
#for example, file_path = './trigonal/split_node_unfiltered.txt'
file_path = sys.argv[1]
lines_per_file = int(sys.argv[2])
split_and_print(file_path, lines_per_file)
