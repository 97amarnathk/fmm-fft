import subprocess

procs = [2, 4, 8, 16, 32, 64]

outfile = open("out_serial.csv", "w")

for i in range(15, 27):
    print "2^"+str(i)
    subprocess.call(["mpiexec", "-np", str(1), 
        "--machinefile", "hostfile/cluster.hosts",
        "fftw_serial.out", str(2**i)], stdout=outfile)