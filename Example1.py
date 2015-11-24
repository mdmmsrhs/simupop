#ARTIFICIAL SELECTION

# #
# Example1.py
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.  Equal events spacing
# 
# phylip, fasta, fstat, popgene format
# 
# #
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM 
# Modified: February 25, 2014 02:11:54 PM 	  
# #

# Import requisite modules

import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

# Set up parameter options
options=[
	{'name': 'n',
	 'default': 1,
	 'label': 'Initial subPops',
	 'type': 'integer',
	 'description': 'The initial number of sub-populations',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'd',
	 'default': 6,
	 'label': 'subPop Split Divisor',
	 'type': 'integer',
	 'description': 'The number of pieces to split the initial sub_population into',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'i',
	 'default': 1e5,
	 'label': 'Initial subPopSize',
	 'type': 'integer',
	 'description': 'The initial number of individuals per sub-population',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'l',
	 'default': 30,
	 'label': 'Loci',
	 'type': 'integer',
	 'description': 'The number of loci',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	  {'name': 'c',
	 'default': 2,
	 'label': 'Chromosomes',
	 'type': 'integer',
	 'description': 'The number of chromosomes',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'g',
	 'default': 10,
	 'label': 'Burnin Gens',
	 'type': 'integer',
	 'description': 'The number of generations before expansion allowed',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 't',
	 'default': 2500,
	 'label': 'Total Generations',
	 'type': 'integer',
	 'description': 'The total number of generations',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'j1',
	 'default': 400,
	 'label': 'eventGens1',
	 'type': 'integer',
	 'description': 'Number of generations before event 1',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'j2',
	 'default': 550,
	 'label': 'eventGens2',
	 'type': 'integer',
	 'description': 'Number of generations before event 2',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'j3',
	 'default': 700,
	 'label': 'eventGens3',
	 'type': 'integer',
	 'description': 'Number of generations before event 3',
	 'validator': simuOpt.valueGT(0)
	 },
	 	 
	 {'name': 'd1',
	 'default': 50,
	 'label': 'eventDur1',
	 'type': 'integer',
	 'description': 'Duration of event 1 in Gens',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'd2',
	 'default': 50,
	 'label': 'eventDur2',
	 'type': 'integer',
	 'description': 'Duration of event 2 in Gens',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'd3',
	 'default': 50,
	 'label': 'eventDur3',
	 'type': 'integer',
	 'description': 'Duration of event 3 in Gens',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'q',
	 'default': 1e5,
	 'label': 'Final PopSize',
	 'type': 'integer',
	 'description': 'The maximum or final sub-population size',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'e',
	 'default': 98,
	 'label': 'propSelf',
	 'type': 'float',
	 'description': 'Proportion of mating that is self-mating',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'r1',
	 'default': 2.5e-4,
	 'label': 'recombInt',
	 'type': 'float',
	 'description': 'Recombination Intensity',
	 'validator': simuOpt.valueGT(0)
	 },
	 	 
	 {'name': 'u',
	 'default': 0.0001,
	 'label': 'forwardMutation',
	 'type': 'float',
	 'description': 'Forward Mutation Rate',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'v',
	 'default': 0.00001,
	 'label': 'backMutation',
	 'type': 'float',
	 'description': 'The Backward Mutation Rate',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'm',
	 'default': 0.01,
	 'label': 'migRate',
	 'type': 'float',
	 'description': 'Migration rate between sub populations',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'h',
	 'default': 0.5,
	 'label': 'domCoeff',
	 'type': 'float',
	 'description': 'The dominance coefficient',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 's1',
	 'default': 0.5,
	 'label': 'selCoeff1',
	 'type': 'float',
	 'description': 'The selection coefficient',
	 'validator': simuOpt.valueGT(0)
	 },
	 
	 {'name': 'z1',
	 'default': 25,
	 'label': 'artSelect1',
	 'type': 'integer',
	 'description': 'Number of generations before sub-Population splitting',
	 'validator': simuOpt.valueGT(0)
	 },
]

## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 2.5e5 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 5.35e-5 # Forward Mutation Rate
v = 5.35e-6# Backward Mutation Rate
m = 0.0005  # Migration rate
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

# Initialise the population
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])

#initialise a splitter function to allow penetrance models
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)

# Define a Function to rapidly re-expand the population if it decreases below i
def demo(pop,gen):
	global q, g
	rate=2**15
	# if .... stop growing
	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
		return ([q for x in pop.subPopSizes(ancGen=-1)])
	else:
		return ([x for x in pop.subPopSizes(ancGen=-1)])
		
# define a function to draw a random sample if individuals and export their sequences		
def sampleAndExport(pop):
    sz = pop.subPopSizes()
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='Example1_a_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False),
    export(sample, format='fstat', output='Example1_a_000000_sample_%d.dat' % pop.dvars().gen, gui = False),
    export(sample, format='phylip', output='Example1_a_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   Example1_a_000000_sample_%d.phy    Example1_a_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   Example1_a_000000_sample_%d.phy    Example1_a_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
# 		sim.DiscardIf(True, subPops=[
#	        (1, 'Affected')], at = j1),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
## Example1_a_000000

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

# Initialise the population
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])

#initialise a splitter function to allow penetrance models
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
    export(sample, format='genepop', output='Example1_a_050000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False),
    export(sample, format='fstat', output='Example1_a_050000_sample_%d.dat' % pop.dvars().gen, gui = False),
    export(sample, format='phylip', output='Example1_a_050000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   Example1_a_050000_sample_%d.phy    Example1_a_050000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   Example1_a_050000_sample_%d.phy    Example1_a_050000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
 		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##Example1_a_050000

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

# Initialise the population
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])

#initialise a splitter function to allow penetrance models
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
    export(sample, format='genepop', output='Example1_a_002100_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False),
    export(sample, format='fstat', output='Example1_a_002100_sample_%d.dat' % pop.dvars().gen, gui = False),
    export(sample, format='phylip', output='Example1_a_002100_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   Example1_a_002100_sample_%d.phy    Example1_a_002100_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   Example1_a_002100_sample_%d.phy    Example1_a_002100_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
# 		sim.DiscardIf(True, subPops=[
#	        (1, 'Affected')], at = j1),
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##Example1_a_002100

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

# Initialise the population
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])

#initialise a splitter function to allow penetrance models
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
    export(sample, format='genepop', output='Example1_a_000037_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False),
    export(sample, format='fstat', output='Example1_a_000037_sample_%d.dat' % pop.dvars().gen, gui = False),
    export(sample, format='phylip', output='Example1_a_000037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   Example1_a_000037_sample_%d.phy    Example1_a_000037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   Example1_a_000037_sample_%d.phy    Example1_a_000037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
# 		sim.DiscardIf(True, subPops=[
#	        (1, 'Affected')], at = j1),
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
##Example1_a_000037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

# Initialise the population
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])

#initialise a splitter function to allow penetrance models
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
    export(sample, format='genepop', output='Example1_a_052100_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False),
    export(sample, format='fstat', output='Example1_a_0052100_sample_%d.dat' % pop.dvars().gen, gui = False),
    export(sample, format='phylip', output='Example1_a_052100_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   Example1_a_052100_sample_%d.phy    Example1_a_052100_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   Example1_a_052100_sample_%d.phy    Example1_a_052100_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
	        (1, 'Affected')], at = j1),	  
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##Example1_a_052100

###################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

# Initialise the population
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])

#initialise a splitter function to allow penetrance models
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
    export(sample, format='genepop', output='Example1_a_052137_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False),
    export(sample, format='fstat', output='Example1_a_052137_sample_%d.dat' % pop.dvars().gen, gui = False),
    export(sample, format='phylip', output='Example1_a_052137_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   Example1_a_052137_sample_%d.phy    Example1_a_052137_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   Example1_a_052137_sample_%d.phy    Example1_a_052137_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
	        (1, 'Affected')], at = j1),	        
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##Example1_a_052137
##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

# Initialise the population
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])

#initialise a splitter function to allow penetrance models
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
    export(sample, format='genepop', output='Example1_a_050037_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False),
    export(sample, format='fstat', output='Example1_a_050037_sample_%d.dat' % pop.dvars().gen, gui = False),
    export(sample, format='phylip', output='Example1_a_050037_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   Example1_a_050037_sample_%d.phy    Example1_a_050037_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   Example1_a_050037_sample_%d.phy    Example1_a_050037_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
	        (1, 'Affected')], at = j1),	        
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
##Example1_a_050037

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

# Initialise the population
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])

#initialise a splitter function to allow penetrance models
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
    export(sample, format='genepop', output='Example1_a_210500_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False),
    export(sample, format='fstat', output='Example1_a_210500_sample_%d.dat' % pop.dvars().gen, gui = False),
    export(sample, format='phylip', output='Example1_a_210500_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   Example1_a_210500_sample_%d.phy    Example1_a_210500_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   Example1_a_210500_sample_%d.phy    Example1_a_210500_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
	        (1, 'Affected')], at = j1),	        
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##Example1_a_210500

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 100    # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

# Initialise the population
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])

#initialise a splitter function to allow penetrance models
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
    export(sample, format='genepop', output='Example1_a_370005_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False),
    export(sample, format='fstat', output='Example1_a_370005_sample_%d.dat' % pop.dvars().gen, gui = False),
    export(sample, format='phylip', output='Example1_a_370005_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   Example1_a_370005_sample_%d.phy    Example1_a_370005_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   Example1_a_370005_sample_%d.phy    Example1_a_370005_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
	        (1, 'Affected')], at = j1),	        
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
##Example1_a_370005

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 100   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

# Initialise the population
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])

#initialise a splitter function to allow penetrance models
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
    export(sample, format='genepop', output='Example1_a_372105_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False),
    export(sample, format='fstat', output='Example1_a_372105_sample_%d.dat' % pop.dvars().gen, gui = False),
    export(sample, format='phylip', output='Example1_a_372105_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   Example1_a_372105_sample_%d.phy    Example1_a_372105_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   Example1_a_372105_sample_%d.phy    Example1_a_372105_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
	        (1, 'Affected')], at = j1),	        
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##Example1_a_372105

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 300   # Number of Steps (Generations) before selection event 2
j3 = 200   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

# Initialise the population
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])

#initialise a splitter function to allow penetrance models
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
    export(sample, format='genepop', output='Example1_a_053721_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False),
    export(sample, format='fstat', output='Example1_a_053721_sample_%d.dat' % pop.dvars().gen, gui = False),
    export(sample, format='phylip', output='Example1_a_053721_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   Example1_a_053721_sample_%d.phy    Example1_a_053721_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   Example1_a_053721_sample_%d.phy    Example1_a_053721_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
	        (0, 'Affected')], at = j1),
 		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1),	        
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##Example1_a_053721

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 550   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

# Initialise the population
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])

#initialise a splitter function to allow penetrance models
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
    export(sample, format='genepop', output='Example1_a_370521_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False),
    export(sample, format='fstat', output='Example1_a_370521_sample_%d.dat' % pop.dvars().gen, gui = False),
    export(sample, format='phylip', output='Example1_a_370521_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   Example1_a_370521_sample_%d.phy    Example1_a_370521_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   Example1_a_370521_sample_%d.phy    Example1_a_370521_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
	        (0, 'Affected')], at = j1),
 		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1),	        
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##Example1_a_370521

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 250   # Number of Steps (Generations) before selection event 2
j3 = 700   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

# Initialise the population
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])

#initialise a splitter function to allow penetrance models
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
    export(sample, format='genepop', output='Example1_a_210537_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False),
    export(sample, format='fstat', output='Example1_a_210537_sample_%d.dat' % pop.dvars().gen, gui = False),
    export(sample, format='phylip', output='Example1_a_210537_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   Example1_a_210537_sample_%d.phy    Example1_a_210537_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   Example1_a_210537_sample_%d.phy    Example1_a_210537_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
	        (0, 'Affected')], at = j1),
 		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1),	        
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##Example1_a_210537

##################################################################################################

z1 = 500   # Number of Steps (Generations) before population splitting
j1 = 400   # Number of Steps (Generations) before selection event 1
j2 = 100   # Number of Steps (Generations) before selection event 2
j3 = 250   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

# Initialise the population
pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])

#initialise a splitter function to allow penetrance models
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
    export(sample, format='genepop', output='Example1_a_213705_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False),
    export(sample, format='fstat', output='Example1_a_213705_sample_%d.dat' % pop.dvars().gen, gui = False),
    export(sample, format='phylip', output='Example1_a_213705_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   Example1_a_213705_sample_%d.phy    Example1_a_213705_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   Example1_a_213705_sample_%d.phy    Example1_a_213705_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
	        (0, 'Affected')], at = j1),
 		sim.DiscardIf(True, subPops=[
	        (1, 'Affected')], at = j1),	        
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
		sim.PyOperator(func=sampleAndExport, at = [0,500,t])
    ],
	gen=(t+1),
)
print 'all done'
##Example1_a_213705