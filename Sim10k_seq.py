##
# sim10k_seq_nat.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 3 Initial subpopns, expanding to 30
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: January 31, 2013 12:15:45 PM   
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Initial Number of Subpopns
d = 10 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.01 # selection coefficient
q = 100000 # Maximum Population Size

a = 1
b = 0

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_nat_a_000000_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_nat_a_000000

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_nat_a_030000_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_nat_a_030000

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_nat_a_000700_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_nat_a_000700

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_nat_a_000012_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_nat_a_000012

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_nat_a_030700_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_nat_a_030700

###################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_nat_a_030712_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_nat_a_030712
##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_nat_a_030012_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_nat_a_030012

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 80   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_nat_a_070300_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_nat_a_070300

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 80    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_nat_a_120003_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_nat_a_120003

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 80   # Number of Steps (Generations) before selection event 2
j3 = 60   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_nat_a_120703_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_nat_a_120703

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_nat_a_031207_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_nat_a_031207
  

##
# sim10k_seq.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 3 Initial subpopns, expanding to 30
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: January 28, 2013 09:21:13 AM  
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Initial Number of Subpopns
d = 10 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.01 # selection coefficient
q = 100000 # Maximum Population Size

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_a_000000_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
	    
	    #Penetrance Matrix#
	    #######################
	    ###   A   C   G   T ###
	    ### A b   a   b   b ###
	    ### C a   a   a   a ###
	    ### G b   a   b   b ###
	    ### T b   a   b   b ###
	    #######################
	    
# 		sim.MaPenetrance(loci=3, penetrance= (1,1,0.0), wildtype = 1),
# 		sim.DiscardIf(True, subPops=[
# 	        (x, 'Affected') for x in list(range(n))], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_a_000000

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_a_030000_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=3, penetrance= (1,1,0.0), wildtype = 1),
		sim.DiscardIf(True, subPops=[
	        (x, 'Affected') for x in list(range(n))], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_a_030000

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_a_000700_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=3, penetrance= (1,1,0.0), wildtype = 1),
# 		sim.DiscardIf(True, subPops=[
# 	        (x, 'Affected') for x in list(range(n))], at = j1),
      sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_a_000700

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_a_000012_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=3, penetrance= (1,1,0.0), wildtype = 1),
# 		sim.DiscardIf(True, subPops=[
# 	        (x, 'Affected') for x in list(range(n))], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_a_000012

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_a_030700_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=3, penetrance= (1,1,0.0), wildtype = 1),
		sim.DiscardIf(True, subPops=[
	        (x, 'Affected') for x in list(range(n))], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_a_030700

###################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_a_030712_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=3, penetrance= (1,1,0.0), wildtype = 1),
		sim.DiscardIf(True, subPops=[
	        (x, 'Affected') for x in list(range(n))], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_a_030712
##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_a_030012_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=3, penetrance= (1,1,0.0), wildtype = 1),
		sim.DiscardIf(True, subPops=[
	        (x, 'Affected') for x in list(range(n))], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_a_030012

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 80   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_a_070300_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=3, penetrance= (1,1,0.0), wildtype = 1),
		sim.DiscardIf(True, subPops=[
	        (x, 'Affected') for x in list(range(n))], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_a_070300

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 80    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_a_120003_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=3, penetrance= (1,1,0.0), wildtype = 1),
		sim.DiscardIf(True, subPops=[
	        (x, 'Affected') for x in list(range(n))], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_a_120003

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 80   # Number of Steps (Generations) before selection event 2
j3 = 60   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_a_120703_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=3, penetrance= (1,1,0.0), wildtype = 1),
		sim.DiscardIf(True, subPops=[
	        (x, 'Affected') for x in list(range(n))], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_a_120703

##################################################################################################

z1 = 120   # Number of Steps (Generations) before population splitting
j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, infoFields=['migrate_to','fitness'])
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
    new_sz = [x//1000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_seq_a_031207_sample_%d.dat' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 3, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 7, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 12, freq=(0.0, 0.2, 0.0, 0.8)),
	    sim.InitGenotype(loci = 16, freq=(0.8, 0.0, 0.2, 0.0)),
	    sim.InitGenotype(loci = [0,15,19], freq=(1.0, 0.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [2,6,8,13], freq=(0.0, 1.0, 0.0, 0.0)),
	    sim.InitGenotype(loci = [14,17,18], freq=(0.0, 0.0, 1.0, 0.0)),
	    sim.InitGenotype(loci = [1,4,5,9,10,11], freq=(0.0, 0.0, 0.0, 1.0)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=3, penetrance= (1,1,0.0), wildtype = 1),
		sim.DiscardIf(True, subPops=[
	        (x, 'Affected') for x in list(range(n))], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(), 
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at=[0,500]),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_seq_a_031207
