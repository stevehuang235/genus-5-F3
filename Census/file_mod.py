import sys
import ast
import os 
from collections import defaultdict 

def read_files_into_dictionary(directory_path):
    curves = defaultdict(list)

    # Check if the given path is a directory
    if not os.path.isdir(directory_path):
        print("Error: Input is not a directory.")
        return None

    # Iterate over each file in the directory
    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename)

        # Check if the current item is a file
        if os.path.isfile(file_path) and '.txt' in file_path:
            with open(file_path, 'r') as file:
                # Read all lines from the file and store them in the dictionary
                lines = file.read().split('\n')
                for ele in lines[:-1]:
                    ele = ele.replace('*','')
                    tmp = ast.literal_eval(ele)
                    curves[tuple(tmp[0])].append(tmp[1])
    
    print('There are {} point-count bins'.format(len(curves)))

    return curves

directory_path = sys.argv[1] + 'data_unfiltered_genus/'
curves = read_files_into_dictionary(directory_path)

def write_dict_to_files(input_dict, output_directory, chunk_size):
    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Split the dictionary into chunks of specified size
    chunks = [list(input_dict.keys())[i:i + chunk_size] for i in range(0, len(input_dict), chunk_size)]

    # Write each chunk to a separate file
    for i, chunk in enumerate(chunks):
        output_filename = os.path.join(output_directory, 'with_genus' + f'_{i + 1}.txt')

        with open(output_filename, 'w') as output_file:
            for key in chunk:
                for elem in curves[key]:
                    tmp = str([list(key), elem])
                    # Find the index of the first '[' character
                    first_bracket_index = tmp.find('[')
                    # Find the index of the last ']' character
                    last_bracket_index = tmp.rfind(']')
                    # Add '*' after the first ']' and in front of the last ']'
                    tmp_modified = tmp[:first_bracket_index+1] + '*' + tmp[first_bracket_index+1:last_bracket_index] + '*' + tmp[last_bracket_index:] + '\n'
                    output_file.write(tmp_modified)

output_directory = sys.argv[1] + 'data_unfiltered_genus_updated/'
chunk_size = int(sys.argv[2])
write_dict_to_files(curves, output_directory, chunk_size)