import sys

def sort_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Splitting lines into two lists based on lines starting with '@'
    at_lines = [line for line in lines if line.startswith('@')]
    other_lines = [line for line in lines if not line.startswith('@')]

    # Sorting other lines based on the integer value found in column 3
    sorted_other_lines = sorted(other_lines, key=lambda line: int(line.split()[3]))

    # Reassembling lines in the desired order
    sorted_lines = at_lines + sorted_other_lines

    # Writing the sorted lines back to the file
    with open(file_path, 'w') as file:
        file.writelines(sorted_lines)

# Example usage:
file_path = sys.argv[1]  # Path to your file
sort_file(file_path)