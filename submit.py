import os,sys
from reaction_tools import get_number_of_reaction

each = 5
inp = "ads.txt"
rxn_num = get_number_of_reaction(inp)
print("rxn_num = %d" % rxn_num)

# remove old files
os.system("rm -rf pbe*")
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

	command = "pjsub -x \"INP={0}\" -x \"ST={1}\" -x \"ED={2}\" run_rxn.sh".format(inp,st,ed)
	print(command)
	os.system(command)

	if last:
		break
	st = ed + 1

