#NATURAL SELECTION

# #
# sim15d_mono_1of2.py
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
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_nat_lomu_a_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_nat_lomu_a_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_nat_lomu_a_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_a_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_nat_lomu_a_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_a_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_nat_lomu_a_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_nat_lomu_b_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_nat_lomu_b_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_nat_lomu_b_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_b_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_nat_lomu_b_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_b_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_nat_lomu_b_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_nat_lomu_c_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_nat_lomu_c_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_nat_lomu_c_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_c_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_nat_lomu_c_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_c_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_nat_lomu_c_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_nat_lomu_d_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_nat_lomu_d_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_nat_lomu_d_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_d_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_nat_lomu_d_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_d_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_nat_lomu_d_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_nat_lomu_e_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_nat_lomu_e_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_nat_lomu_e_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_e_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_nat_lomu_e_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_e_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_nat_lomu_e_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_nat_lomu_f_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_nat_lomu_f_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_nat_lomu_f_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_f_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_nat_lomu_f_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_f_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_nat_lomu_f_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_nat_lomu_g_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_nat_lomu_g_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_nat_lomu_g_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_g_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_nat_lomu_g_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_g_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_nat_lomu_g_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_nat_lomu_h_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_nat_lomu_h_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_nat_lomu_h_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_h_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_nat_lomu_h_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_h_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_nat_lomu_h_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_nat_lomu_i_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_nat_lomu_i_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_nat_lomu_i_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_i_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_nat_lomu_i_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_i_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_nat_lomu_i_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_nat_lomu_j_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_nat_lomu_j_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_nat_lomu_j_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_j_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_nat_lomu_j_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_nat_lomu_j_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_nat_lomu_j_370005

#NATURAL SELECTION

# #
# sim15d_mono.py
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
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_nat_lomu_a_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_nat_lomu_a_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_nat_lomu_a_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_a_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_nat_lomu_a_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_a_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_nat_lomu_a_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_nat_lomu_b_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_nat_lomu_b_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_nat_lomu_b_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_b_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_nat_lomu_b_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_b_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_nat_lomu_b_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_nat_lomu_c_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_nat_lomu_c_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_nat_lomu_c_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_c_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_nat_lomu_c_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_c_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_nat_lomu_c_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_nat_lomu_d_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_nat_lomu_d_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_nat_lomu_d_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_d_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_nat_lomu_d_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_d_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_nat_lomu_d_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_nat_lomu_e_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_nat_lomu_e_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_nat_lomu_e_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_e_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_nat_lomu_e_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_e_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_nat_lomu_e_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_nat_lomu_f_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_nat_lomu_f_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_nat_lomu_f_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_f_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_nat_lomu_f_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_f_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_nat_lomu_f_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_nat_lomu_g_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_nat_lomu_g_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_nat_lomu_g_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_g_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_nat_lomu_g_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_g_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_nat_lomu_g_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_nat_lomu_h_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_nat_lomu_h_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_nat_lomu_h_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_h_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_nat_lomu_h_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_h_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_nat_lomu_h_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_nat_lomu_i_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_nat_lomu_i_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_nat_lomu_i_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_i_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_nat_lomu_i_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_i_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_nat_lomu_i_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_nat_lomu_j_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_nat_lomu_j_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_nat_lomu_j_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_j_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_nat_lomu_j_370005_sample_%d.phy    sim15d_mono_seq_subset_nat_lomu_j_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),


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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_nat_lomu_j_370005

#NATURAL SELECTION

# #
# sim15d_diphyletic.py
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
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_nat_lomu_a_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_nat_lomu_a_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_nat_lomu_a_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_a_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_nat_lomu_a_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_a_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),



	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_nat_lomu_a_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_nat_lomu_b_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_nat_lomu_b_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_nat_lomu_b_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_b_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_nat_lomu_b_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_b_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),



	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_nat_lomu_b_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_nat_lomu_c_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_nat_lomu_c_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_nat_lomu_c_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_c_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_nat_lomu_c_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_c_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),



	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_nat_lomu_c_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_nat_lomu_d_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_nat_lomu_d_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_nat_lomu_d_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_d_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_nat_lomu_d_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_d_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),



	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_nat_lomu_d_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_nat_lomu_e_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_nat_lomu_e_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_nat_lomu_e_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_e_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_nat_lomu_e_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_e_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),



	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_nat_lomu_e_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_nat_lomu_f_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_nat_lomu_f_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_nat_lomu_f_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_f_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_nat_lomu_f_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_f_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),



	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_nat_lomu_f_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_nat_lomu_g_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_nat_lomu_g_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_nat_lomu_g_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_g_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_nat_lomu_g_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_g_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),



	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_nat_lomu_g_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_nat_lomu_h_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_nat_lomu_h_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_nat_lomu_h_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_h_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_nat_lomu_h_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_h_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),



	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_nat_lomu_h_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_nat_lomu_i_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_nat_lomu_i_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_nat_lomu_i_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_i_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_nat_lomu_i_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_i_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),



	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_nat_lomu_i_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_nat_lomu_j_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_nat_lomu_j_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_nat_lomu_j_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_j_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_nat_lomu_j_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_nat_lomu_j_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),



	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_nat_lomu_j_370005
#NATURAL SELECTION

# #
# sim15d_diphyletic_sep.py
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
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_a_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_a_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_nat_lomu_a_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_a_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_nat_lomu_a_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_a_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),




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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_nat_lomu_a_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_b_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_b_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_nat_lomu_b_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_b_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_nat_lomu_b_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_b_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),




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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_nat_lomu_b_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_c_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_c_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_nat_lomu_c_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_c_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_nat_lomu_c_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_c_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),




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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_nat_lomu_c_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_d_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_d_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_nat_lomu_d_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_d_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_nat_lomu_d_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_d_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),




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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_nat_lomu_d_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_e_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_e_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_nat_lomu_e_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_e_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_nat_lomu_e_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_e_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),




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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_nat_lomu_e_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_f_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_f_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_nat_lomu_f_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_f_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_nat_lomu_f_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_f_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),




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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_nat_lomu_f_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_g_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_g_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_nat_lomu_g_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_g_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_nat_lomu_g_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_g_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),




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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_nat_lomu_g_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_h_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_h_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_nat_lomu_h_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_h_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_nat_lomu_h_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_h_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),




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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_nat_lomu_h_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_i_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_i_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_nat_lomu_i_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_i_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_nat_lomu_i_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_i_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),




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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_nat_lomu_i_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_j_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_nat_lomu_j_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_nat_lomu_j_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_j_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_nat_lomu_j_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_nat_lomu_j_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),




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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_nat_lomu_j_370005


#ARTIFICIAL SELECTION

# #
# sim15d_mono_1of2.py
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
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_lomu_a_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_lomu_a_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_lomu_a_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_a_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_lomu_a_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_a_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_lomu_a_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_lomu_b_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_lomu_b_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_lomu_b_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_b_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_lomu_b_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_b_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_lomu_b_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_lomu_c_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_lomu_c_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_lomu_c_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_c_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_lomu_c_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_c_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_lomu_c_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_lomu_d_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_lomu_d_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_lomu_d_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_d_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_lomu_d_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_d_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_lomu_d_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_lomu_e_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_lomu_e_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_lomu_e_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_e_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_lomu_e_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_e_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_lomu_e_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_lomu_f_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_lomu_f_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_lomu_f_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_f_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_lomu_f_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_f_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_lomu_f_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_lomu_g_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_lomu_g_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_lomu_g_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_g_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_lomu_g_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_g_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_lomu_g_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_lomu_h_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_lomu_h_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_lomu_h_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_h_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_lomu_h_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_h_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_lomu_h_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_lomu_i_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_lomu_i_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_lomu_i_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_i_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_lomu_i_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_i_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_lomu_i_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_1of2_seq_subset_lomu_j_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_1of2_seq_subset_lomu_j_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_1of2_seq_subset_lomu_j_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_j_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_1of2_seq_subset_lomu_j_370005_sample_%d.phy    sim15d_mono_1of2_seq_subset_lomu_j_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_1of2_lomu_j_370005

#ARTIFICIAL SELECTION

# #
# sim15d_mono.py
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
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_lomu_a_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_lomu_a_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_lomu_a_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_a_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_lomu_a_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_a_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_lomu_a_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_lomu_b_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_lomu_b_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_lomu_b_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_b_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_lomu_b_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_b_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_lomu_b_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_lomu_c_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_lomu_c_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_lomu_c_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_c_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_lomu_c_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_c_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_lomu_c_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_lomu_d_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_lomu_d_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_lomu_d_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_d_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_lomu_d_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_d_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_lomu_d_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_lomu_e_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_lomu_e_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_lomu_e_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_e_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_lomu_e_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_e_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_lomu_e_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_lomu_f_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_lomu_f_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_lomu_f_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_f_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_lomu_f_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_f_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_lomu_f_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_lomu_g_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_lomu_g_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_lomu_g_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_g_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_lomu_g_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_g_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_lomu_g_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_lomu_h_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_lomu_h_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_lomu_h_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_h_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_lomu_h_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_h_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_lomu_h_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_lomu_i_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_lomu_i_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_lomu_i_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_i_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_lomu_i_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_i_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_lomu_i_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_mono_seq_subset_lomu_j_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_mono_seq_subset_lomu_j_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_mono_seq_subset_lomu_j_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_j_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_mono_seq_subset_lomu_j_370005_sample_%d.phy    sim15d_mono_seq_subset_lomu_j_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_mono_lomu_j_370005

#ARTIFICIAL SELECTION

# #
# sim15d_diphyletic.py
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
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_lomu_a_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_lomu_a_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_lomu_a_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_a_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_lomu_a_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_a_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_lomu_a_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_lomu_b_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_lomu_b_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_lomu_b_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_b_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_lomu_b_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_b_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_lomu_b_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_lomu_c_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_lomu_c_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_lomu_c_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_c_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_lomu_c_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_c_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_lomu_c_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_lomu_d_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_lomu_d_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_lomu_d_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_d_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_lomu_d_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_d_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_lomu_d_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_lomu_e_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_lomu_e_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_lomu_e_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_e_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_lomu_e_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_e_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_lomu_e_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_lomu_f_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_lomu_f_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_lomu_f_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_f_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_lomu_f_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_f_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_lomu_f_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_lomu_g_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_lomu_g_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_lomu_g_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_g_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_lomu_g_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_g_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_lomu_g_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_lomu_h_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_lomu_h_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_lomu_h_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_h_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_lomu_h_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_h_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_lomu_h_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_lomu_i_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_lomu_i_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_lomu_i_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_i_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_lomu_i_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_i_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_lomu_i_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_seq_subset_lomu_j_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_seq_subset_lomu_j_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_seq_subset_lomu_j_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_j_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_seq_subset_lomu_j_370005_sample_%d.phy    sim15d_diphyletic_seq_subset_lomu_j_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_lomu_j_370005
#ARTIFICIAL SELECTION

# #
# sim15d_diphyletic_sep.py
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
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_lomu_a_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_lomu_a_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_lomu_a_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_a_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_lomu_a_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_a_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1+20)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_lomu_a_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_lomu_b_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_lomu_b_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_lomu_b_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_b_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_lomu_b_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_b_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1+20)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_lomu_b_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_lomu_c_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_lomu_c_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_lomu_c_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_c_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_lomu_c_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_c_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1+20)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_lomu_c_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_lomu_d_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_lomu_d_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_lomu_d_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_d_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_lomu_d_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_d_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1+20)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_lomu_d_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_lomu_e_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_lomu_e_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_lomu_e_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_e_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_lomu_e_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_e_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1+20)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_lomu_e_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_lomu_f_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_lomu_f_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_lomu_f_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_f_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_lomu_f_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_f_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1+20)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_lomu_f_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_lomu_g_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_lomu_g_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_lomu_g_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_g_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_lomu_g_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_g_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1+20)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_lomu_g_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_lomu_h_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_lomu_h_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_lomu_h_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_h_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_lomu_h_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_h_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1+20)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_lomu_h_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_lomu_i_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_lomu_i_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_lomu_i_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_i_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_lomu_i_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_i_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1+20)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_lomu_i_370005
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

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
    export(sample, format='fstat', output='sim15d_diphyletic_sep_seq_subset_lomu_j_370005_sample_%d.dat' % pop.dvars().gen, gui=False),
    export(sample, format='phylip', output='sim15d_diphyletic_sep_seq_subset_lomu_j_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim15d_diphyletic_sep_seq_subset_lomu_j_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_j_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim15d_diphyletic_sep_seq_subset_lomu_j_370005_sample_%d.phy    sim15d_diphyletic_sep_seq_subset_lomu_j_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1+20)
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##sim15d_diphyletic_sep_lomu_j_370005
