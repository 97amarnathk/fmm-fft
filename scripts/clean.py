outfile = open("../plotting/cleaned/cluster_out_fmmfft_32_15.txt", "w")
infile =  open("../plotting/raw_data/cluster_out_fmmfft_32_15.txt", "r")

for l in infile:
    if l[0] is not "[" and l[0].isdigit():
        outfile.write(l)