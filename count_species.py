from reaction_tools import *
#
# count species
#
reactionfile = "reaction.txt"
speciesfile  = "species.txt"

fspecies = open(speciesfile, "w")

(r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)

species = []
for imol in r_ads:
	species.append(imol)

for imol in p_ads:
	species.append(imol)

species = [item for sublist in species for item in sublist]

fspecies.write(str(species))
fspecies.close()

