## to open Python in R
library(reticulate)
use_condaenv("bioenv", conda = "/opt/conda/condabin/conda", required = TRUE)
# py_list_packages()
reticulate::repl_python()



# test libraries
import matplotlib
import tqdm
import msprime
import tskit
import demesdraw
import csv

############################ Anopheles gambiae simulation ############################
import numpy
from matplotlib import pyplot as plt

with open('Results1', 'w') as output_file:
	results = []
	for i in [50_000, 75_000, 100_000, 125_000, 150_000, 175_000, 200_000]:

		NpopN1=i
		NpopN2=NpopN1/30
	#NpopN1=30*NpopN2

	# add demography
		demography = msprime.Demography()
		tmp=demography.add_population(name="N1", initial_size=NpopN1)
		tmp=demography.add_population(name="N2", initial_size=NpopN2)

		tmp=demography.add_population(name="ANC", initial_size=7_000_000)
		tmp=demography.add_population_split(time=1000, derived=["N1", "N2"], ancestral="ANC") # split 1k generations ago

		ts = msprime.sim_ancestry(
			{"N1": 10, "N2": 10}, 
			demography=demography, 
			recombination_rate=8.4e-9, # as in humans
			sequence_length=1_000,
			random_seed=1234)
		#print(ts)
		#k=k+1
	 
	 
	# we can add mutations
		mts = msprime.sim_mutations(ts, rate=3.5e-9, random_seed=1234)
		#print(mts.tables.sites)
		
		# Define the samples between which Fst will be calculated
		pop_id = {p.metadata["name"]: p.id for p in mts.populations()}
		sample_sets=[mts.samples(pop_id["N1"]), mts.samples(pop_id["N2"])]
		
		print(mts.Fst(sample_sets))
		thisResults = []
		thisResults.append(mts.Fst(sample_sets))
		thisResults.append(mts.divergence(sample_sets))
		thisResults.extend(mts.Tajimas_D(sample_sets))
		thisResults.extend(mts.segregating_sites(sample_sets))
		results.append(thisResults)
		print(thisResults)
		numpy.array(results)
		output_file.write('\t'.join(map(str, thisResults)) + '\n')


Results_of_msprime=read.table("Results1")
View(Results_of_msprime)
colnames(Results_of_msprime) <- c("Fst", "dxy", "tajima1", "tajima2", "segsites1", "segsites2")
rownames(Results_of_msprime) <- c("50_000", "75_000", "100_000", "125_000", "150_000", "175_000", "200_000")

Fst, dxy, segsites1, segsites2, pi1, pi2, tajima1, tajima2
0.2134, 0.0978, 0.3797, 0.1013, 0.0914, 0.0355, 0.2847, 2.0788




















 results = mts


# show SNPs
for variant in mts.variants():
    print(variant)

# visualise the haplotypes
samples = mts.samples()
for sample_id, h in zip(samples, mts.haplotypes(samples=samples)):
    pop = ts.node(sample_id).population
    print(f"Sample {sample_id:<2} ({ts.population(pop).metadata['name']:^5}): {h}")


# visualise the site frequency spectrum
plt.clf()
afs = mts.allele_frequency_spectrum()
plt.bar(range(mts.num_samples + 1), afs)
plt.title("Allele frequency spectrum")
plt.show()


# Define the samples between which Fst will be calculated
pop_id = {p.metadata["name"]: p.id for p in mts.populations()}
sample_sets=[mts.samples(pop_id["N1"]), mts.samples(pop_id["N2"])]

print(mts.Fst(sample_sets))

