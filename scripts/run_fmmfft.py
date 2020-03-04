import subprocess

procs = [2, 4, 8, 16, 32, 64]
boxes = 32
terms = 15

outfile = open("out_fmmfft_32_15.csv", "w")

for i in range(15, 27):
    for p in procs:
        print "2^"+str(i), "processors ", p
        subprocess.call(["mpiexec", "-np", str(p), 
            "--machinefile", "hostfile/cluster.hosts",
            "--map-by", "node",
            "fmmfft.out", str(2**i), str(boxes), str(terms)], stdout=outfile)