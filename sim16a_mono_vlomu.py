
#### ARTIFICIAL SELECTION

# #
# sim15b_v_low_mu_mono.py
# 
# # Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.  Equal events spacing

# 2 Initial subpopns, expanding to 6
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# phylip, fasta format
# 
# #
# Author: Richard Stephens
# Created: January 22, 2014 08:30:54 AM 
# Modified: January 22, 2014 08:39:39 AM   	  
# #


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
# Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean recombination intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)

def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_a_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_a_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_a_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_vlomu_a_000000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_a_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_a_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_a_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_a_050000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_a_002100_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_a_002100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_002100_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_a_002100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_002100_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# # 		sim.DiscardIf(True, subPops=[
# # 	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_a_002100
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_a_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_a_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_a_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# # 		sim.DiscardIf(True, subPops=[
# # 	        (0, 'Affected')], at = j1),
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_a_000037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_a_052100_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_a_052100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_052100_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_a_052100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_052100_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_a_052100

###################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_a_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_a_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_a_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_a_052137
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_a_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_a_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_a_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_a_050037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_a_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_a_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_a_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_a_210500

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_a_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_a_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_a_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_a_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 100   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_a_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_a_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_a_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_a_372105

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_a_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_a_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_a_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_a_053721

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_a_210537_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_a_210537_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_210537_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_a_210537_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_210537_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_a_210537

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_a_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_a_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_a_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_a_370521

#### ARTIFICIAL SELECTION

# #
# sim15b_v_low_mu_mono.py
# 
# # Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.  Equal events spacing

# 2 Initial subpopns, expanding to 6
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# phylip, fasta format
# 
# #
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM 
# Modified: July 30, 2013 09:04:47 AM 	  
# #


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)

def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_b_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_b_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_b_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_b_000000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_b_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_b_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_b_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_b_050000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_b_002100_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_b_002100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_002100_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_b_002100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_002100_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# # 		sim.DiscardIf(True, subPops=[
# # 	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_b_002100
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_b_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_b_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_b_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# # 		sim.DiscardIf(True, subPops=[
# # 	        (0, 'Affected')], at = j1),
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_b_000037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_b_052100_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_b_052100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_052100_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_b_052100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_052100_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_b_052100

###################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_b_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_b_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_b_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_b_052137
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_b_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_b_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_b_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_b_050037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_b_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_b_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_b_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_b_210500

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_b_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_b_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_b_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_b_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 100   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_b_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_b_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_b_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_b_372105

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_b_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_b_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_b_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_b_053721

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_b_210537_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_b_210537_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_210537_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_b_210537_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_210537_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_b_210537

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_b_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_b_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_b_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_b_370521



#### ARTIFICIAL SELECTION

# #
# sim15b_v_low_mu_mono.py
# 
# # Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.  Equal events spacing

# 2 Initial subpopns, expanding to 6
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# phylip, fasta format
# 
# #
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM 
# Modified: July 30, 2013 09:04:47 AM 	  
# #


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)

def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_c_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_c_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_c_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_c_000000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_c_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_c_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_c_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_c_050000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_c_002100_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_c_002100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_002100_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_c_002100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_002100_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# # 		sim.DiscardIf(True, subPops=[
# # 	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_c_002100
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_c_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_c_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_c_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# # 		sim.DiscardIf(True, subPops=[
# # 	        (0, 'Affected')], at = j1),
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_c_000037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_c_052100_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_c_052100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_052100_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_c_052100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_052100_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_c_052100

###################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_c_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_c_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_c_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_c_052137
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_c_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_c_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_c_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_c_050037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_c_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_c_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_c_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_c_210500

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_c_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_c_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_c_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_c_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 100   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_c_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_c_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_c_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_c_372105

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_c_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_c_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_c_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_c_053721

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_c_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_c_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_c_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_c_370521

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_c_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_c_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_c_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_c_370521

#### ARTIFICIAL SELECTION

# #
# sim15b_v_low_mu_mono.py
# 
# # Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.  Equal events spacing

# 2 Initial subpopns, expanding to 6
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# phylip, fasta format
# 
# #
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM 
# Modified: July 30, 2013 09:04:47 AM 	  
# #


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)

def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_d_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_d_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_d_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_d_000000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_d_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_d_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_d_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_d_050000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_d_002100_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_d_002100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_002100_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_d_002100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_002100_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# # 		sim.DiscardIf(True, subPops=[
# # 	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_d_002100
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_d_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_d_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_d_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# # 		sim.DiscardIf(True, subPops=[
# # 	        (0, 'Affected')], at = j1),
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_d_000037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_d_052100_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_d_052100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_052100_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_d_052100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_052100_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_d_052100

###################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_d_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_d_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_d_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_d_052137
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_d_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_d_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_d_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_d_050037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_d_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_d_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_d_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_d_210500

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_d_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_d_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_d_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_d_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 100   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_d_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_d_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_d_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_d_372105

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_d_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_d_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_d_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_d_053721

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_d_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_d_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_d_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_d_370521

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_d_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_d_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_d_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_d_370521


#### ARTIFICIAL SELECTION

# #
# sim15b_v_low_mu_mono.py
# 
# # Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.  Equal events spacing

# 2 Initial subpopns, expanding to 6
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# phylip, fasta format
# 
# #
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM 
# Modified: July 30, 2013 09:04:47 AM 	  
# #


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)

def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_e_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_e_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_e_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_e_000000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_e_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_e_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_e_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_e_050000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_e_002100_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_e_002100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_002100_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_e_002100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_002100_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# # 		sim.DiscardIf(True, subPops=[
# # 	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_e_002100
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_e_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_e_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_e_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# # 		sim.DiscardIf(True, subPops=[
# # 	        (0, 'Affected')], at = j1),
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_vlomu_e_000037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_e_052100_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_e_052100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_052100_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_e_052100_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_052100_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_e_052100

###################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_e_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_e_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_e_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_e_052137
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_e_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_e_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_e_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_e_050037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_e_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_e_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_e_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_e_210500

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_e_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_e_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_e_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_e_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 100   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_e_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_e_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_e_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_e_372105

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_e_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_e_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_e_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_e_053721

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_e_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_e_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_e_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_e_370521

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_e_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_e_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_e_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_e_370521

import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
# Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)
# Settings ###################################################################################

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 100   # Number of Steps (Generations) before selection event 2
j3 = 250  # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_a_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_a_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_a_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_a_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_a_213705

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 100   # Number of Steps (Generations) before selection event 2
j3 = 250  # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_b_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_b_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_b_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_b_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_b_213705

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 100   # Number of Steps (Generations) before selection event 2
j3 = 250  # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='+_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_c_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_c_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_c_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_vlomu_c_213705

#################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 100   # Number of Steps (Generations) before selection event 2
j3 = 250  # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

#################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_d_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_d_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_d_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_d_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_vlomu_d_213705

#################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 100   # Number of Steps (Generations) before selection event 2
j3 = 250  # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

#################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_vlomu_e_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_vlomu_e_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_vlomu_e_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_vlomu_e_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_vlomu_e_213705

##NATURAL SELECTION

##
# sim15b_v_low_mu_mono_1of2_nat.py
#
# # Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.  Equal events spacing
# 2 initial subpopns, expanding to 6
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) in 1 sub pop before splitting
# If other selection events before that at locus 4, this occurs before the split also

# phylip, fasta format
# 
#
##
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM
# Modified: July 30, 2013 09:04:33 AM	   
##


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)

def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_000000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_000000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_b_000000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_050000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_050000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_b_050000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_001200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_001200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_001200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_001200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_001200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_b_002100
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_000037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_000037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_b_000037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_051200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_051200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_051200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_051200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_051200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_b_052100

###################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_052137_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_052137_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_b_052137
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_050037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_050037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_b_050037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_210500_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_210500_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_b_210500

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_370005_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_370005_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_b_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 100   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_372105_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_372105_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_b_372105

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_053721_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_053721_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_1of2_nat_vlomu_b_053721

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_210537_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_210537_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_210537_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_210537_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_210537_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_1of2_nat_vlomu_b_210537

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_370521_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_370521_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_1of2_nat_vlomu_b_370521


##NATURAL SELECTION

##
# sim15b_v_low_mu_mono_1of2_nat.py
#
# # Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.  Equal events spacing
# 2 initial subpopns, expanding to 6
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) in 1 sub pop before splitting
# If other selection events before that at locus 4, this occurs before the split also

# phylip, fasta format
# 
#
##
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM
# Modified: July 30, 2013 09:04:33 AM	   
##


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)

def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_000000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_000000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_c_000000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_050000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_050000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_c_050000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_001200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_001200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_001200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_001200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_001200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_c_002100
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_000037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_000037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_c_000037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_051200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_051200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_051200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_051200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_051200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_c_052100

###################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_052137_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_052137_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_c_052137
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_050037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_050037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_c_050037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_210500_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_210500_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_c_210500

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_370005_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_370005_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_c_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 100   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_372105_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_372105_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_c_372105

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_053721_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_053721_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_1of2_nat_vlomu_c_053721

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_210537_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_210537_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_210537_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_210537_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_210537_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_1of2_nat_vlomu_c_210537

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_370521_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_370521_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_1of2_nat_vlomu_c_370521


##NATURAL SELECTION

##
# sim15b_v_low_mu_mono_1of2_nat.py
#
# # Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.  Equal events spacing
# 2 initial subpopns, expanding to 6
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) in 1 sub pop before splitting
# If other selection events before that at locus 4, this occurs before the split also

# phylip, fasta format
# 
#
##
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM
# Modified: July 30, 2013 09:04:33 AM	   
##


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)

def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_000000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_000000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_1of2_nat_vlomu_d_000000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_050000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_050000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_d_050000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_001200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_001200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_001200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_001200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_001200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_d_002100
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_000037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_000037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_d_000037
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_051200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_051200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_051200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_051200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_051200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_d_052100
# 
# ###################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_052137_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_052137_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_d_052137
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_050037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_050037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_d_050037
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 250   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_210500_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_210500_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_d_210500
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 250    # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_370005_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_370005_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_d_370005
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 250   # Number of Steps (Generations) before selection event 2
# j3 = 100   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_372105_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_372105_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_d_372105
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 300   # Number of Steps (Generations) before selection event 2
# j3 = 200   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_053721_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_053721_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# #sim15b_v_low_mu_mono_1of2_nat_vlomu_d_053721
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 250   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_210537_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_210537_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_210537_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_210537_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_210537_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# #sim15b_v_low_mu_mono_1of2_nat_vlomu_d_210537
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 250   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_370521_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_370521_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# #sim15b_v_low_mu_mono_1of2_nat_vlomu_d_370521
# 
# 
# ##NATURAL SELECTION
# 
# ##
# # sim15b_v_low_mu_mono_1of2_nat.py
# #
# # # Script to simulate a population with n subpops of i individuals, with population rebound, 
# # natural selection and sampling.  Selection occuring at different periods.  Equal events spacing
# # 2 initial subpopns, expanding to 6
# # output as genotypes for Powermarker
# # Selection at locus 4 by natural selection (bottleneck) in 1 sub pop before splitting
# # If other selection events before that at locus 4, this occurs before the split also
# 
# # phylip, fasta format
# # 
# #
# ##
# # Author: Richard Stephens
# # Created: July 30, 2013 09:04:33 AM
# # Modified: July 30, 2013 09:04:33 AM	   
# ##
# 
# 
# import simuOpt
# import simuPOP as sim
# import math, os
# from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
# from simuPOP.sampling import drawRandomSample
# 
#        
# ## Settings ###################################################################################
# n = 2   # Initial Number of Subpopns
# d = 3 # Divisor for subpopn splitting
# i = 1e5 # Number of Indivs/Subpopn
# l = 30    # Number of Loci per Chromosome
# c = 2 # Number of Chromosomes
# g = 10 # Number of Steps (Generations) before expansion
# t = 2500 # Total Number of Steps (Generations)
# u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
# v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
# m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
# e = 98 # Proportion of selfing
# s1 = 0.008 # selection coefficient
# q = 1e5 # Maximum Population Size
# r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# 
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_000000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_000000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_e_000000
# # 
# # ##################################################################################################
# # 
# # z1 = 500   # Number of Steps (Generations) before population splitting
# # j1 = 400   # Number of Steps (Generations) before selection event 1
# # j2 = 550   # Number of Steps (Generations) before selection event 2
# # j3 = 700   # Number of Steps (Generations) before selection event 3
# # 
# # d1 = 50 # duration of selection event 1 (in generations)
# # d2 = 50 # duration of selection event 2 (in generations)
# # d3 = 50 # duration of selection event 3 (in generations)
# # 
# # ##################################################################################################
# # 
# # pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# # pop.setVirtualSplitter(sim.ProductSplitter([
# #     sim.AffectionSplitter()
# #     ])
# # )
# # z=20
# # def demo(pop,gen):
# # 	global q, g
# # 	rate=2**15
# # 	# if .... stop growing
# # 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# # 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# # 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# # 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# # 	else:
# # 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# # 		
# # def sampleAndExport(pop):
# #     sz = pop.subPopSizes()
# #     new_sz = [x//2000 for x in sz]
# #     sample = drawRandomSample(pop, new_sz)
# #     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
# #     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_050000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
# #     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_050000_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
# #     return True
# # 
# # simu = sim.Simulator(pop, rep=1)
# # simu.evolve(
# #     initOps=[
# # 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# # 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# # 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# # 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# # 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# # 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# # 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# # 	],
# # 	preOps=[
# # 		sim.SNPMutator(u = u, v = v),
# # 		sim.SNPMutator(u = u, v = v),
# # 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# # 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# # 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# # #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# # #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
# #     ],
# # 	
# #     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
# #     matingScheme=sim.HeteroMating([
# #         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
# #         sim.HomoMating(
# #             chooser=sim.CombinedParentsChooser(
# #                 sim.RandomParentChooser(),
# #                 sim.RandomParentChooser()),
# #             generator=sim.OffspringGenerator(
# #                 ops= [
# # 	                sim.MendelianGenoTransmitter(),
# # 	                sim.Recombinator(intensity=r1)
# #                 ],
# #             ),
# #             weight=100-e)
# #         ],
# #     subPopSize=demo
# #     ),
# # 	postOps=[
# # 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
# #     ],
# # 	gen=(t+1),
# # )
# # print 'all done'
# # ##sim15b_v_low_mu_mono_1of2_nat_vlomu_e_050000
# # 
# # ##################################################################################################
# # 
# # z1 = 500   # Number of Steps (Generations) before population splitting
# # j1 = 400   # Number of Steps (Generations) before selection event 1
# # j2 = 550   # Number of Steps (Generations) before selection event 2
# # j3 = 700   # Number of Steps (Generations) before selection event 3
# # 
# # d1 = 50 # duration of selection event 1 (in generations)
# # d2 = 50 # duration of selection event 2 (in generations)
# # d3 = 50 # duration of selection event 3 (in generations)
# # 
# # ##################################################################################################
# # 
# # pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# # pop.setVirtualSplitter(sim.ProductSplitter([
# #     sim.AffectionSplitter()
# #     ])
# # )
# # z=20
# # def demo(pop,gen):
# # 	global q, g
# # 	rate=2**15
# # 	# if .... stop growing
# # 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# # 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# # 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# # 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# # 	else:
# # 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# # 		
# # def sampleAndExport(pop):
# #     sz = pop.subPopSizes()
# #     new_sz = [x//2000 for x in sz]
# #     sample = drawRandomSample(pop, new_sz)
# #     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_001200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
# #     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_001200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_001200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
# #     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_001200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_001200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
# #     return True
# # 
# # simu = sim.Simulator(pop, rep=1)
# # simu.evolve(
# #     initOps=[
# # 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# # 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# # 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# # 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# # 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# # 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# # 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# # 	],
# # 	preOps=[
# # 		sim.SNPMutator(u = u, v = v),
# # 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# # 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# # 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# # #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
# #     ],
# # 	
# #     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
# #     matingScheme=sim.HeteroMating([
# #         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
# #         sim.HomoMating(
# #             chooser=sim.CombinedParentsChooser(
# #                 sim.RandomParentChooser(),
# #                 sim.RandomParentChooser()),
# #             generator=sim.OffspringGenerator(
# #                 ops= [
# # 	                sim.MendelianGenoTransmitter(),
# # 	                sim.Recombinator(intensity=r1)
# #                 ],
# #             ),
# #             weight=100-e)
# #         ],
# #     subPopSize=demo
# #     ),
# # 	postOps=[
# # 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
# #     ],
# # 	gen=(t+1),
# # )
# # print 'all done'
# # ##sim15b_v_low_mu_mono_1of2_nat_vlomu_e_002100
# # 
# # ##################################################################################################
# # 
# # z1 = 500   # Number of Steps (Generations) before population splitting
# # j1 = 400   # Number of Steps (Generations) before selection event 1
# # j2 = 550   # Number of Steps (Generations) before selection event 2
# # j3 = 700   # Number of Steps (Generations) before selection event 3
# # 
# # d1 = 50 # duration of selection event 1 (in generations)
# # d2 = 50 # duration of selection event 2 (in generations)
# # d3 = 50 # duration of selection event 3 (in generations)
# # 
# # ##################################################################################################
# # 
# # pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# # pop.setVirtualSplitter(sim.ProductSplitter([
# #     sim.AffectionSplitter()
# #     ])
# # )
# # z=20
# # def demo(pop,gen):
# # 	global q, g
# # 	rate=2**15
# # 	# if .... stop growing
# # 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# # 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# # 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# # 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# # 	else:
# # 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# # 		
# # def sampleAndExport(pop):
# #     sz = pop.subPopSizes()
# #     new_sz = [x//2000 for x in sz]
# #     sample = drawRandomSample(pop, new_sz)
# #     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
# #     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_000037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
# #     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_000037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
# #     return True
# # 
# # simu = sim.Simulator(pop, rep=1)
# # simu.evolve(
# #     initOps=[
# # 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# # 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# # 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# # 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# # 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# # 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# # 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# # 	],
# # 	preOps=[
# # 		sim.SNPMutator(u = u, v = v),
# # 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# # 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# # 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# # #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
# #     ],
# # 	
# #     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
# #     matingScheme=sim.HeteroMating([
# #         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
# #         sim.HomoMating(
# #             chooser=sim.CombinedParentsChooser(
# #                 sim.RandomParentChooser(),
# #                 sim.RandomParentChooser()),
# #             generator=sim.OffspringGenerator(
# #                 ops= [
# # 	                sim.MendelianGenoTransmitter(),
# # 	                sim.Recombinator(intensity=r1)
# #                 ],
# #             ),
# #             weight=100-e)
# #         ],
# #     subPopSize=demo
# #     ),
# # 	postOps=[
# # 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
# #     ],
# # 	gen=(t+1),
# # )
# # print 'all done'
# # ##sim15b_v_low_mu_mono_1of2_nat_vlomu_e_000037
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_051200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_051200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_051200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_051200_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_051200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_e_052100
# 
# ###################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_052137_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_052137_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_e_052137
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_050037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_050037_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_e_050037
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 250   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_210500_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_210500_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_e_210500
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 250    # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_370005_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_370005_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_e_370005
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 250   # Number of Steps (Generations) before selection event 2
# j3 = 100   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_372105_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_372105_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_1of2_nat_vlomu_e_372105
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 300   # Number of Steps (Generations) before selection event 2
# j3 = 200   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_053721_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_053721_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# #sim15b_v_low_mu_mono_1of2_nat_vlomu_e_053721
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 250   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_210537_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_210537_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_210537_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_210537_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_210537_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# #sim15b_v_low_mu_mono_1of2_nat_vlomu_e_210537
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 250   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_370521_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_370521_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# #sim15b_v_low_mu_mono_1of2_nat_vlomu_e_370521
# 
# import simuOpt
# import simuPOP as sim
# import math, os
# from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
# from simuPOP.sampling import drawRandomSample
# 
#        
# ## Settings ###################################################################################
# n = 2   # Initial Number of Subpopns
# d = 3 # Divisor for subpopn splitting
# i = 1e5 # Number of Indivs/Subpopn
# l = 30    # Number of Loci per Chromosome
# c = 2 # Number of Chromosomes
# g = 10 # Number of Steps (Generations) before expansion
# t = 2500 # Total Number of Steps (Generations)
# u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
# v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
# m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
# e = 98 # Proportion of selfing
# s1 = 0.008 # selection coefficient
# q = 1e5 # Maximum Population Size
# r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 100   # Number of Steps (Generations) before selection event 2
# j3 = 250   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_a_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_a_213705_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_a_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_a_213705_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_a_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# #sim15b_v_low_mu_mono_1of2_nat_vlomu_a_213705
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 100   # Number of Steps (Generations) before selection event 2
# j3 = 250   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_213705_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_213705_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_b_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# #sim15b_v_low_mu_mono_1of2_nat_vlomu_b_213705
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 100   # Number of Steps (Generations) before selection event 2
# j3 = 250   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_213705_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_213705_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_c_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# #sim15b_v_low_mu_mono_1of2_nat_vlomu_c_213705
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 100   # Number of Steps (Generations) before selection event 2
# j3 = 250   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_213705_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_213705_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_d_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# #sim15b_v_low_mu_mono_1of2_nat_vlomu_d_213705
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 100   # Number of Steps (Generations) before selection event 2
# j3 = 250   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_213705_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_213705_sample_%d.phy    sim15b_v_low_mu_mono_1of2_seq_subset_nat_vlomu_e_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#         sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# #sim15b_v_low_mu_mono_1of2_nat_vlomu_e_213705



##NATURAL SELECTION

##
# sim15b_v_low_mu_mono_nat.py
#
# # Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.  Equal events spacing
# 2 initial subpopns, expanding to 6
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) in 1 sub pop before splitting
# If other selection events before that at locus 4, this occurs before the split also

# phylip, fasta format
# 
#
##
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM
# Modified: July 30, 2013 09:04:33 AM	   
##


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)

def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_a_000000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_nat_vlomu_a_050000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_001200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_001200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_001200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_001200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_001200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_nat_vlomu_a_002100
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_nat_vlomu_a_000037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_051200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_051200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_051200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_051200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_051200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_a_052100

###################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_a_052137
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_a_050037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_a_210500

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_a_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 100   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_a_372105

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_a_053721

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_210537_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_210537_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_210537_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_210537_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_210537_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_a_210537

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_a_370521


##NATURAL SELECTION

##
# sim15b_v_low_mu_mono_nat.py
#
# # Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.  Equal events spacing
# 2 initial subpopns, expanding to 6
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) in 1 sub pop before splitting
# If other selection events before that at locus 4, this occurs before the split also

# phylip, fasta format
# 
#
##
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM
# Modified: July 30, 2013 09:04:33 AM	   
##


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)

def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_b_000000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_nat_vlomu_b_050000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_001200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_001200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_001200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_001200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_001200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_nat_vlomu_b_002100
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_nat_vlomu_b_000037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_051200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_051200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_051200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_051200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_051200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_b_052100

###################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_b_052137
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_b_050037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_b_210500

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_b_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 100   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_b_372105

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_b_053721

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_210537_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_210537_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_210537_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_210537_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_210537_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_b_210537

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_b_370521


##NATURAL SELECTION

##
# sim15b_v_low_mu_mono_nat.py
#
# # Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.  Equal events spacing
# 2 initial subpopns, expanding to 6
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) in 1 sub pop before splitting
# If other selection events before that at locus 4, this occurs before the split also

# phylip, fasta format
# 
#
##
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM
# Modified: July 30, 2013 09:04:33 AM	   
##


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)

def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_c_000000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_nat_vlomu_c_050000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_001200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_001200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_001200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_001200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_001200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_nat_vlomu_c_002100
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_nat_vlomu_c_000037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_051200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_051200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_051200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_051200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_051200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_c_052100

###################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_c_052137
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_c_050037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_c_210500

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_c_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 100   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_c_372105

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_c_053721

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_210537_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_210537_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_210537_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_210537_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_210537_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_c_210537

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_c_370521


##NATURAL SELECTION

##
# sim15b_v_low_mu_mono_nat.py
#
# # Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.  Equal events spacing
# 2 initial subpopns, expanding to 6
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) in 1 sub pop before splitting
# If other selection events before that at locus 4, this occurs before the split also

# phylip, fasta format
# 
#
##
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM
# Modified: July 30, 2013 09:04:33 AM	   
##


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)

def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_d_000000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_nat_vlomu_d_050000

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_001200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_001200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_001200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_001200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_001200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_d_002100

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
      sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_d_000037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_051200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_051200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_051200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_051200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_051200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_d_052100

###################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_d_052137
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_d_050037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_d_210500

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_d_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 100   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_d_372105

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_d_053721

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_210537_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_210537_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_210537_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_210537_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_210537_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_d_210537

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_d_370521


##NATURAL SELECTION

##
# sim15b_v_low_mu_mono_nat.py
#
# # Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.  Equal events spacing
# 2 initial subpopns, expanding to 6
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) in 1 sub pop before splitting
# If other selection events before that at locus 4, this occurs before the split also

# phylip, fasta format
# 
#
##
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM
# Modified: July 30, 2013 09:04:33 AM	   
##


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)

def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_000000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_e_000000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_050000_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
# 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_nat_vlomu_e_050000
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_001200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_001200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_001200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_001200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_001200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_nat_vlomu_e_002100
# 
# ##################################################################################################
# 
# z1 = 500   # Number of Steps (Generations) before population splitting
# j1 = 400   # Number of Steps (Generations) before selection event 1
# j2 = 550   # Number of Steps (Generations) before selection event 2
# j3 = 700   # Number of Steps (Generations) before selection event 3
# 
# d1 = 50 # duration of selection event 1 (in generations)
# d2 = 50 # duration of selection event 2 (in generations)
# d3 = 50 # duration of selection event 3 (in generations)
# 
# ##################################################################################################
# 
# pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
# pop.setVirtualSplitter(sim.ProductSplitter([
#     sim.AffectionSplitter()
#     ])
# )
# z=20
# def demo(pop,gen):
# 	global q, g
# 	rate=2**15
# 	# if .... stop growing
# 	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
# 		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
# 	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
# 		return ([q for x in pop.subPopSizes(ancGen=-1)])
# 	else:
# 		return ([x for x in pop.subPopSizes(ancGen=-1)])
# 		
# def sampleAndExport(pop):
#     sz = pop.subPopSizes()
#     new_sz = [x//2000 for x in sz]
#     sample = drawRandomSample(pop, new_sz)
#     export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
#     os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
#     os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_000037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
#     return True
# 
# simu = sim.Simulator(pop, rep=1)
# simu.evolve(
#     initOps=[
# 	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
# 	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
# 	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
# 	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
# 	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
# 	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
# 	],
# 	preOps=[
# 		sim.SNPMutator(u = u, v = v),
# 		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
# 		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
# 	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# # 		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
# #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
#     ],
# 	
#     #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
#     matingScheme=sim.HeteroMating([
#         sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
#         sim.HomoMating(
#             chooser=sim.CombinedParentsChooser(
#                 sim.RandomParentChooser(),
#                 sim.RandomParentChooser()),
#             generator=sim.OffspringGenerator(
#                 ops= [
# 	                sim.MendelianGenoTransmitter(),
# 	                sim.Recombinator(intensity=r1)
#                 ],
#             ),
#             weight=100-e)
#         ],
#     subPopSize=demo
#     ),
# 	postOps=[
# 		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
#     ],
# 	gen=(t+1),
# )
# print 'all done'
# ##sim15b_v_low_mu_mono_nat_vlomu_e_000037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_051200_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_051200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_051200_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_051200_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_051200_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_e_052100

###################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_052137_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_e_052137
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_050037_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_e_050037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_210500_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
#       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_e_210500

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_370005_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_e_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 100   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_372105_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15b_v_low_mu_mono_nat_vlomu_e_372105

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_053721_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_e_053721

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_210537_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_210537_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_210537_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_210537_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_210537_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_e_210537

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_370521_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_e_370521

import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 1e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.3e-5 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 5.3e-6 # Backward Mutation Rate # estimate from Kahler et al. 1984
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 1e5 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 100   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_a_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_a_213705

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 100   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_b_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_b_213705

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 100   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_c_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_c_213705

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 100   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_d_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_d_213705

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 100   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_213705_sample_%d.phy    sim15b_v_low_mu_mono_seq_subset_nat_vlomu_e_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
        sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = [0,600,800,t])
    ],
	gen=(t+1),
)
print 'all done'
#sim15b_v_low_mu_mono_nat_vlomu_e_213705