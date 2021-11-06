from reactions import Reaction, Reactions
from tinydb import TinyDB
import os
import glob
for file in glob.glob("test.json"):
	os.remove(file)

rxns = Reactions.from_csv("test.csv")
print(rxns[1]._reaction_str)
