import sys, pickle, argparse
from reaction_tools import *
#
# count species
#
parser = argparse.ArgumentParser()
parser.add_argument("--reactionfile", required=True, help="file with elementary reactions")
argvs = parser.parse_args()
infile  = argvs.reactionfile

speciesfile  = "species.txt"
speciesfile_json = "species.json"

fspecies = open(speciesfile, "w")

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(infile)

species = []
rxnnum = len(r_ads)
for irxn in range(rxnnum):
	#
	# reactant
	#
	for imol, mol in enumerate(r_ads[irxn]):
		mol = mol[0]
		mol = remove_side_and_flip(mol)

		site = r_site[irxn][imol]
		site = site[0]
		try:
			site, site_pos = site.split(".")
		except:
			site_pos = "x1y1"

		if site != 'gas':
			mol += "_surf"  # Not distinguish different sites. Reconsideration may be needed.

		species.append(mol)

	#
	# product
	#
	for imol, mol in enumerate(p_ads[irxn]):
		mol = mol[0]
		mol = remove_side_and_flip(mol)

		site = p_site[irxn][imol]
		site = site[0]
		try:
			site, site_pos = site.split(".")
		except:
			site_pos = "x1y1"

		if site != 'gas':
			mol += "_surf"

		species.append(mol)

species = list(set(species))  # remove duplication
species.sort()

surf = [i for i in species if 'surf' in i]
gas  = [i for i in species if 'surf' not in i]

species = gas + surf

print("Number of compounds: %d, Number of surface species: %d" % (len(species), len(surf)))

# put 'surf' to the last
for item in species:
	if item == 'surf':
		species.remove(item)
		species.append(item)

fspecies.write(str(species))
fspecies.close()
pickle.dump(species, open("species.pickle", "wb"))
