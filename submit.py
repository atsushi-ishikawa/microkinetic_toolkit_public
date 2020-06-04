import os,sys
from reaction_tools import get_number_of_reaction

# usage: submit.py inputfile [number of jobs]

submit_sh = "run_kyushu.sh"
argvs = sys.argv
inp   = argvs[1]
rxn_num = get_number_of_reaction(inp)

if len(argvs)==3:
	# -- specify number of jobs
	num_jobs = int(argvs[2])
	each = rxn_num // num_jobs
else:
	# --  specify reactions per jobs
	each = 5

print("rxn_num = %d, rxns per jobs = %d" % (rxn_num, each))

# remove old files
os.system("rm -rf beef-vdw*")
os.system("rm -rf tsdir")
os.system("rm tmp.db")
os.system("rm std*")

st = 0
ed = 1

last = False
while True:
	ed = st + each
	if ed >= rxn_num:
		ed = rxn_num
		last = True

	command = "pjsub -x \"INP={0}\" -x \"ST={1}\" -x \"ED={2}\" {3}".format(inp,st,ed,submit_sh)
	print(command)
	os.system(command)

	if last:
		break
	st = ed + 1

