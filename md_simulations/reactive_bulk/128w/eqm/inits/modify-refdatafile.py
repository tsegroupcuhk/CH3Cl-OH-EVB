input_file = "refdatafile.txt"  # Replace with the actual path to your input file
output_file = "new.txt"  # Replace with the desired path for the output file

with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
    for line in f_in:
        line = line.strip()
        if line:
            parts = line.split()
            file_path = parts[1]
            index = int(file_path.split("/")[-1].split(".")[0])
            force_content = "./refdata/force-{}.xyz".format(index)
            modified_line = "{} {}".format(line, force_content)
            f_out.write(modified_line + "\n")

print("Output file generated: {}".format(output_file))
