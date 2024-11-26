import os
import sys
import gzip

__author__ = 'Rob Edwards'

def join(files, raw_output, norm_output):    
    raw = {}
    norm = {}
    alltax = set()
    for f in files:
        raw[f] = {}
        norm[f] = {}
        opener = gzip.open if f.endswith('.gz') else open
        with opener(f, 'rt') as fn:
            for l in fn:
                p = l.strip().split("\t")
                #print (p[2], p[3], p[4])
                if len(p) >= 5:  # Check for valid lines
                    raw[f][p[2].strip()] = p[3]
                    norm[f][p[2].strip()] = p[4]
                    alltax.add(p[2].strip())

    with gzip.open(raw_output, 'wt') as out:
        print("\t".join(["taxonomy"] + files), file=out)  # Write the header
        for t in sorted(alltax):
            # Start the line with the taxonomy ID
            line = [t]
            for f in files:
                if t not in raw[f]:
                    line.append("0")  # Append 0 if the taxonomy is not found
                else:
                    line.append(str(raw[f][t]))  # Append the raw count for the taxonomy
            # Write the complete line and add a newline at the end
            print("\t".join(line), file=out)  # This will automatically add a newline after each line

    with gzip.open(norm_output, 'wt') as out:
        # Print the header
        print("\t".join(["taxonomy"] + files), file=out)
        
        # Iterate over sorted taxonomies
        for t in sorted(alltax):
            # Print the taxonomy and start the line with its name
            line = [t]
            
            # Append the data for each file
            for f in files:
                if t not in norm[f]:
                    line.append("0")  # Add 0 if the taxonomy is not found in the file
                else:
                    line.append(str(norm[f][t]))  # Add the value for the taxonomy in the file

            # Print the entire line and ensure a newline at the end (this is automatic with print)
            print("\t".join(line), file=out)



#join(snakemake.input.files, snakemake.output.raw_output, snakemake.output.norm_output)   

if __name__ == "__main__":
    # Read Nextflow-supplied environment variable for input files
    input_files = sys.argv[1:-2]  # All arguments except the last two
    output_raw = sys.argv[-2]    # Second-to-last argument
    output_norm = sys.argv[-1]   # Last argument

    join(input_files, output_raw, output_norm)

