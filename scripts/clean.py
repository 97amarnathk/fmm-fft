outfile = open("out_trial_clean.csv", "w")
infile =  open("out_trial.csv", "r")

for l in infile:
    if l[0] is not "[" and l[0].isdigit():
        outfile.write(l)