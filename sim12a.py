##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################

##
# sim12a.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different times.  Most sites monomorphic
# to give "sequence-like" data output using "phylip" Exporter format.  Non-selected polymorphic sites added
# to simulate selective sweeps.
#
# simplified with fewer sites for a network analysis cf sim12.py
# 
# Switched to MaSelector as well in order to use ACGT allele names, which require that InitGenotype
# takes four arguments, rather than two as with 0,1 alleles.  Set wildtype= to allele which is being 
# selected FOR and fitness of all but the wildtype homozygote to 1-s1 (qv).
#
# With and without recombination, muation and migration, several values of s, 5 reps. 
#
# 
#
##
# Author: Richard Stephens
# Created: December 19, 2012 09:36:55 AM 
# Modified: December 20, 2012 01:35:36 PM      
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##sim12a_s002_a

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

print('check1')
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_a_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_a_000000_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
# 		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# 		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
    	sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_a_000000

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_a_041700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_a_041700_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# 		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),    
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_a_041700

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 500   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_a_170400_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_a_170400_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_a_170400

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_a_041731_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_a_041731_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_a_041731

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 1500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_a_311704_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_a_311704_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_a_311704

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 1500   # Number of Steps (Generations) before selection event 2
j3 = 750   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_a_043117_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_a_043117_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_a_043117

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_a_040031_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_a_040031_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_a_040031

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 1500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_a_310004_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_a_310004_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_a_310004

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 500   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_a_170431_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_a_170431_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_a_170431

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 1500   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_a_310417_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_a_310417_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_a_310417

##
# sim12.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different times.  Most sites monomorphic
# to give "sequence-like" data output using "phylip" Exporter format.  Non-selected polymorphic sites added
# to simulate selective sweeps.
# 
# Switched to MaSelector as well in order to use ACGT allele names, which require that InitGenotype
# takes four arguments, rather than two as with 0,1 alleles.  Set wildtype= to allele which is being 
# selected FOR and fitness of all but the wildtype homozygote to 1-s1 (qv).
#
# With and without recombination, muation and migration, several values of s, 5 reps. 
#
# 
#
##
# Author: Richard Stephens
# Created: December 19, 2012 09:36:55 AM 
# Modified: December 19, 2012 09:37:03 AM     
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##sim12a_s002_a

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_b_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_b_000000_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
# 		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# 		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
    	sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_b_000000

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_b_041700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_b_041700_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# 		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),    
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_b_041700

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 500   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_b_170400_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_b_170400_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_b_170400

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_b_041731_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_b_041731_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_b_041731

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 1500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_b_311704_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_b_311704_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_b_311704

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 1500   # Number of Steps (Generations) before selection event 2
j3 = 750   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_b_043117_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_b_043117_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_b_043117

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_b_040031_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_b_040031_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_b_040031

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 1500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_b_310004_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_b_310004_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_b_310004

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 500   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_b_170431_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_b_170431_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_b_170431

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 1500   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_b_310417_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_b_310417_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_b_310417

##
# sim12.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different times.  Most sites monomorphic
# to give "sequence-like" data output using "phylip" Exporter format.  Non-selected polymorphic sites added
# to simulate selective sweeps.
# 
# Switched to MaSelector as well in order to use ACGT allele names, which require that InitGenotype
# takes four arguments, rather than two as with 0,1 alleles.  Set wildtype= to allele which is being 
# selected FOR and fitness of all but the wildtype homozygote to 1-s1 (qv).
#
# With and without recombination, muation and migration, several values of s, 5 reps. 
#
# 
#
##
# Author: Richard Stephens
# Created: December 19, 2012 09:36:55 AM 
# Modified: December 19, 2012 09:37:03 AM     
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##sim12a_s002_a

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_c_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_c_000000_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
# 		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# 		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
    	sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_c_000000

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_c_041700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_c_041700_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# 		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),    
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_c_041700

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 500   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_c_170400_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_c_170400_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_c_170400

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_c_041731_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_c_041731_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_c_041731

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 1500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_c_311704_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_c_311704_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_c_311704

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 1500   # Number of Steps (Generations) before selection event 2
j3 = 750   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_c_043117_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_c_043117_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_c_043117

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_c_040031_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_c_040031_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_c_040031

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 1500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_c_310004_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_c_310004_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_c_310004

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 500   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_c_170431_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_c_170431_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_c_170431

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 1500   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_c_310417_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_c_310417_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_c_310417

##
# sim12.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different times.  Most sites monomorphic
# to give "sequence-like" data output using "phylip" Exporter format.  Non-selected polymorphic sites added
# to simulate selective sweeps.
# 
# Switched to MaSelector as well in order to use ACGT allele names, which require that InitGenotype
# takes four arguments, rather than two as with 0,1 alleles.  Set wildtype= to allele which is being 
# selected FOR and fitness of all but the wildtype homozygote to 1-s1 (qv).
#
# With and without recombination, muation and migration, several values of s, 5 reps. 
#
# 
#
##
# Author: Richard Stephens
# Created: December 19, 2012 09:36:55 AM 
# Modified: December 19, 2012 09:37:03 AM     
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##sim12a_s002_a

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_d_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_d_000000_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
# 		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# 		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
    	sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_d_000000

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_d_041700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_d_041700_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# 		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),    
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_d_041700

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 500   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_d_170400_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_d_170400_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_d_170400

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_d_041731_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_d_041731_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_d_041731

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 1500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_d_311704_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_d_311704_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_d_311704

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 1500   # Number of Steps (Generations) before selection event 2
j3 = 750   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_d_043117_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_d_043117_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_d_043117

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_d_040031_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_d_040031_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_d_040031

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 1500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_d_310004_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_d_310004_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_d_310004

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 500   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_d_170431_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_d_170431_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_d_170431

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 1500   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_d_310417_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_d_310417_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_d_310417

##
# sim12.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different times.  Most sites monomorphic
# to give "sequence-like" data output using "phylip" Exporter format.  Non-selected polymorphic sites added
# to simulate selective sweeps.
# 
# Switched to MaSelector as well in order to use ACGT allele names, which require that InitGenotype
# takes four arguments, rather than two as with 0,1 alleles.  Set wildtype= to allele which is being 
# selected FOR and fitness of all but the wildtype homozygote to 1-s1 (qv).
#
# With and without recombination, muation and migration, several values of s, 5 reps. 
#
# 
#
##
# Author: Richard Stephens
# Created: December 19, 2012 09:36:55 AM 
# Modified: December 19, 2012 09:37:03 AM     
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##sim12a_s002_a

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_e_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_e_000000_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
# 		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# 		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
    	sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_e_000000

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_e_041700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_e_041700_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
# 		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),    
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_e_041700

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 500   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_e_170400_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_e_170400_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_e_170400

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_e_041731_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_e_041731_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_e_041731

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 1500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_e_311704_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_e_311704_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_e_311704

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 1500   # Number of Steps (Generations) before selection event 2
j3 = 750   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_e_043117_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_e_043117_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_e_043117

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_e_040031_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_e_040031_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_e_040031

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 1500   # Number of Steps (Generations) before selection event 1
j2 = 750   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_e_310004_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_e_310004_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
# 		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_e_310004

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 500   # Number of Steps (Generations) before selection event 2
j3 = 1500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_e_170431_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_e_170431_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_e_170431

import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 20   # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2000 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.02 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-5 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##################################################################################################

j1 = 750   # Number of Steps (Generations) before selection event 1
j2 = 1500   # Number of Steps (Generations) before selection event 2
j3 = 500   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,l,1))*c, infoFields=['migrate_to','fitness'])
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
    export(sample, format='phylip', output='sim12a_s002_e_310417_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    export(sample, format='fstat', output='sim12a_s002_e_310417_sample_%d.dat' % pop.dvars().gen, gui=False),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	     sim.InitGenotype(loci = [17,30], freq=(0.0, 0.80, 0.20, 0.0)), #S
         sim.InitGenotype(loci = [3,5,16,18,31], freq=(0.2, 0.0, 0.8, 0.0)), #R
         sim.InitGenotype(loci = [4,32], freq=(0.0, 0.2, 0.0, 0.8)), #Y
         sim.InitGenotype(loci = [0,10,19,24,25,35,36], freq=(1.0, 0.0, 0.0, 0.0)), #A
         sim.InitGenotype(loci = [1,7,12,14,21,28,37,39], freq=(0.0, 1.0, 0.0, 0.0)), #C
         sim.InitGenotype(loci = [6,11,20,26,27,33,34,38], freq=(0.0, 0.0, 1.0, 0.0)), #G
         sim.InitGenotype(loci = [2,8,9,13,15,22,23,29], freq=(0.0, 0.0, 0.0, 1.0)), #T
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MaSelector(loci=4, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),
		sim.MaSelector(loci=17, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
		sim.MaSelector(loci=31, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
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
        sim.PyOperator(func = sampleAndExport, at = t),
    ],
	gen=(t+1),
)
print 'all done'
##sim12a_s002_e_310417

