import sys

inputname = sys.argv[1]

partitions = {}
with open(inputname) as inputhandle:
	for line in inputhandle:
		if line[:7] == "charset":
			ls = line.strip().split(" ")[-1].rstrip(";").split("=")
			partitions[ls[0]] = ls[1]

with open("mrbblock.txt", "w") as outputhandle:
	print("begin mrbayes;", file=outputhandle)
	print("\tset autoclose=yes nowarn=yes;", file=outputhandle)
	print("\tpartition prt =", len(partitions), ":", ", ".join(partitions.keys()) ,";", file=outputhandle)
	print("\tset partition = prt;", file=outputhandle)
	print("\tlset applyto=(all) nst=6 rates=invgamma;", file=outputhandle)
	print("\tunlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all);", file=outputhandle)
	print("\tprset applyto=(all) ratepr=variable;", file=outputhandle)
	print("\tmcmcp ngen= 10000000 relburnin=yes burninfrac=0.25 printfreq=1000 samplefreq=1000 nchains=4 savebrlens=yes;", file=outputhandle)
	print("\tmcmc;", file=outputhandle)
	print("\tsump;", file=outputhandle)
	print("\tsumt;", file=outputhandle)
	print("end;", file=outputhandle)