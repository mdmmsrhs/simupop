##
# sim10k_mono_1of2_nat.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) in 1 sub pop before splitting
# If other selection events before that at locus 4, this occurs before the split also

# genepop, fstat format
# 
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 08, 2013 10:48:27 AM 	   
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_a_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_a_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_a_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_a_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_a_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_a_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_a_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_a_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_a_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_a_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_a_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_a_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_a_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_a_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_a_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_a_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_a_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_a_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_a_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_a_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_a_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_a_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_a_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_a_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_a_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_a_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_a_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_a_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_a_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_a_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_a_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_a_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_a_043414
  

##
# sim10k_mono.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 08, 2013 10:48:27 AM 	  
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_a_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_a_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_a_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_a_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_a_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_a_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_a_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_a_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_a_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_a_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_a_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_a_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_a_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_a_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_a_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_a_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_a_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_a_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_a_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_a_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_a_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_a_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_a_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_a_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_a_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_a_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_a_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_a_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_a_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_a_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_a_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_a_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_a_043414



##
# sim10k_mono_nat.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 1 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:49:26 AM   
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_a_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_a_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_a_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_a_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_a_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_a_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_a_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_a_00001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_a_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_a_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_a_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_a_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_a_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_a_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_a_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_a_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_a_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_a_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_a_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_a_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_a_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_a_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_a_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_a_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_a_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_a_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_a_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_a_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_a_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_a_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_a_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_a_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_a_043414
  

#
# sim10k_mono.py
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 1 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:49:35 AM  
#


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
# Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

#################################################################################################

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_a_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_a_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_a_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_a_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_a_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_a_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_a_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_a_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_a_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_a_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_a_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_a_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_a_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_a_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_a_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_a_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_a_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_a_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_a_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_a_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_a_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_a_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_a_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_a_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_a_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_a_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_a_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_a_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_a_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_a_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_a_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_a_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_a_043414

##
# sim10k_mono_1of2_nat.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) in 1 sub pop before splitting
# If other selection events before that at locus 4, this occurs before the split also

# genepop, fstat format
# 
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:50:59 AM 	   
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_b_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_b_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_b_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_b_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_b_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_b_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_b_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_b_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_b_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_b_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_b_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_b_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_b_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_b_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_b_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_b_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_b_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_b_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_b_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_b_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_b_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_b_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_b_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_b_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_b_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_b_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_b_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_b_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_b_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_b_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_b_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_b_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_b_043414
  

##
# sim10k_mono.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 08, 2013 10:48:27 AM 	  
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_b_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_b_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_b_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_b_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_b_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_b_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_b_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_b_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_b_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_b_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_b_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_b_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_b_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_b_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_b_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_b_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_b_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_b_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_b_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_b_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_b_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_b_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_b_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_b_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_b_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_b_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_b_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_b_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_b_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_b_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_b_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_b_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_b_043414



##
# sim10k_mono_nat.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 1 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:49:26 AM   
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_b_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_b_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_b_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_b_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_b_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_b_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_b_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_b_00001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_b_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_b_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_b_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_b_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_b_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_b_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_b_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_b_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_b_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_b_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_b_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_b_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_b_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_b_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_b_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_b_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_b_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_b_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_b_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_b_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_b_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_b_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_b_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_b_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_b_043414
  

#
# sim10k_mono.py
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 1 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:49:35 AM  
#


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
# Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

#################################################################################################

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_b_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_b_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_b_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_b_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_b_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_b_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_b_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_b_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_b_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_b_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_b_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_b_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_b_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_b_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_b_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_b_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_b_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_b_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_b_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_b_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_b_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_b_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_b_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_b_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_b_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_b_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_b_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_b_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_b_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_b_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_b_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_b_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_b_043414

##
# sim10k_mono_1of2_nat.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) in 1 sub pop before splitting
# If other selection events before that at locus 4, this occurs before the split also

# genepop, fstat format
# 
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:50:59 AM 	   
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_c_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_c_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_c_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_c_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_c_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_c_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_c_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_c_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_c_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_c_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_c_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_c_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_c_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_c_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_c_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_c_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_c_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_c_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_c_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_c_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_c_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_c_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_c_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_c_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_c_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_c_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_c_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_c_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_c_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_c_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_c_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_c_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_c_043414
  

##
# sim10k_mono.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 08, 2013 10:48:27 AM 	  
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_c_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_c_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_c_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_c_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_c_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_c_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_c_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_c_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_c_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_c_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_c_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_c_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_c_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_c_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_c_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_c_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_c_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_c_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_c_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_c_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_c_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_c_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_c_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_c_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_c_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_c_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_c_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_c_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_c_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_c_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_c_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_c_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_c_043414



##
# sim10k_mono_nat.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 1 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:49:26 AM   
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_c_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_c_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_c_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_c_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_c_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_c_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_c_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_c_00001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_c_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_c_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_c_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_c_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_c_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_c_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_c_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_c_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_c_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_c_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_c_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_c_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_c_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_c_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_c_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_c_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_c_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_c_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_c_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_c_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_c_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_c_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_c_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_c_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_c_043414
  

#
# sim10k_mono.py
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 1 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:49:35 AM  
#


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
# Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

#################################################################################################

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_c_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_c_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_c_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_c_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_c_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_c_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_c_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_c_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_c_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_c_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_c_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_c_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_c_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_c_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_c_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_c_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_c_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_c_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_c_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_c_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_c_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_c_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_c_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_c_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_c_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_c_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_c_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_c_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_c_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_c_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_c_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_c_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_c_043414

##
# sim10k_mono_1of2_nat.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) in 1 sub pop before splitting
# If other selection events before that at locus 4, this occurs before the split also

# genepop, fstat format
# 
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:50:59 AM 	   
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_d_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_d_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_d_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_d_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_d_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_d_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_d_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_d_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_d_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_d_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_d_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_d_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_d_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_d_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_d_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_d_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_d_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_d_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_d_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_d_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_d_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_d_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_d_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_d_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_d_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_d_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_d_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_d_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_d_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_d_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_d_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_d_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_d_043414
  

##
# sim10k_mono.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 08, 2013 10:48:27 AM 	  
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_d_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_d_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_d_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_d_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_d_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_d_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_d_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_d_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_d_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_d_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_d_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_d_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_d_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_d_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_d_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_d_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_d_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_d_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_d_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_d_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_d_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_d_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_d_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_d_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_d_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_d_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_d_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_d_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_d_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_d_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_d_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_d_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_d_043414



##
# sim10k_mono_nat.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 1 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:49:26 AM   
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_d_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_d_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_d_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_d_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_d_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_d_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_d_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_d_00001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_d_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_d_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_d_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_d_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_d_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_d_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_d_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_d_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_d_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_d_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_d_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_d_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_d_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_d_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_d_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_d_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_d_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_d_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_d_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_d_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_d_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_d_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_d_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_d_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_d_043414
  

#
# sim10k_mono.py
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 1 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:49:35 AM  
#


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
# Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

#################################################################################################

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_d_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_d_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_d_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_d_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_d_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_d_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_d_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_d_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_d_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_d_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_d_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_d_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_d_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_d_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_d_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_d_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_d_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_d_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_d_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_d_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_d_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_d_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_d_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_d_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_d_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_d_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_d_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_d_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_d_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_d_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_d_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_d_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_d_043414

##
# sim10k_mono_1of2_nat.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) in 1 sub pop before splitting
# If other selection events before that at locus 4, this occurs before the split also

# genepop, fstat format
# 
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:50:59 AM 	   
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_e_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_e_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_e_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_e_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_e_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_e_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_e_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_e_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_e_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_e_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_e_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_e_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_e_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_e_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_e_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_e_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_e_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_e_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_e_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_e_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_e_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_e_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_e_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_e_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_e_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_e_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_e_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_e_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_e_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_e_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_nat_e_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_nat_e_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_nat_e_043414
  

##
# sim10k_mono.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 08, 2013 10:48:27 AM 	  
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_e_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_e_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_e_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_e_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_e_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_e_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_e_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_e_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_e_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_e_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_e_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_e_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_e_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_e_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_e_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_e_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_e_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_e_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_e_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_e_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_e_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_e_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_e_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_e_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_e_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_e_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_e_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_e_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_e_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_e_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_1of2_e_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_1of2_e_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_e_043414



##
# sim10k_mono_nat.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 1 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by natural selection (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:49:26 AM   
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_e_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_e_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_e_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_e_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_e_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_e_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_e_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_e_00001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_e_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_e_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_e_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_e_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_e_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_e_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_e_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_e_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_e_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_e_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_e_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_e_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_e_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_e_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_e_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_e_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_e_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_e_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_e_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_e_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_e_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_e_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_nat_e_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_nat_e_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_nat_e_043414
  

#
# sim10k_mono.py
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 1 Initial subpopns, expanding to 14
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:49:35 AM  
#


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
# Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

#################################################################################################

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_e_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_e_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_e_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_e_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_e_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_e_040000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_e_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_e_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
      sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_e_001400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_e_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_e_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
# 		sim.DiscardIf(True, subPops=[
# 	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_e_000034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_e_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_e_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_e_041400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_e_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_e_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_e_041434
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_e_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_e_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_e_040034

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_e_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_e_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_e_140400

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_e_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_e_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
#       sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_e_340004

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_e_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_e_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_e_341404

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10k_mono_e_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10k_mono_e_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaPenetrance(loci=4, penetrance=[1, 1, 0.0]),
		sim.DiscardIf(True, subPops=[
	        (0, 'Affected')], at = j1),
        sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_e_043414


##
# sim10k_mono_1of2_seq_nat.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 10
# output as genotypes for Powermarker
# art selection in 1 of 3 original popns

# phylip format
# 
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:49:43 AM   
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)


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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_nat_a_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
# 		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_nat_a_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_nat_a_030000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_nat_a_030000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_nat_a_000700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
# 		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_nat_a_000700

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_nat_a_000012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
# 		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_nat_a_000012

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_nat_a_030700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_nat_a_030700

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_nat_a_030712_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_nat_a_030712
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_nat_a_030012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_nat_a_030012

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_nat_a_070300_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_nat_a_070300

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_nat_a_120003_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_nat_a_120003

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_nat_a_120703_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_nat_a_120703

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_nat_a_031207_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
		sim.MaSelector(loci=3, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_nat_a_031207
  

##
# sim10k_mono_1of2_seq.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 10
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# phylip format
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:49:52 AM 
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_a_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
# 	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_a_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_a_030000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_a_030000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_a_000700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
# 	        (0, 'Affected')], at = j1),
      sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_a_000700

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_a_000012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
# 	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_a_000012

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_a_030700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_a_030700

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_a_030712_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_a_030712
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_a_030012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_a_030012

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_a_070300_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_a_070300

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_a_120003_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_a_120003

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_a_120703_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_a_120703

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_1of2_seq_a_031207_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_1of2_seq_a_031207



##
# sim10k_mono_seq_nat.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 1 Initial subpopns, expanding to 10
# output as genotypes for Powermarker
# Selection at locus 3 by natural selection (bottleneck) before splitting
# If other selection events before that at locus 3, this occurs before the split also

# selectio befiore split in all three base populations
# 
# phylip format
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:50:27 AM    
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)


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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_nat_a_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_nat_a_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_nat_a_030000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_nat_a_030000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_nat_a_000700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_nat_a_000700

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_nat_a_000012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_nat_a_000012

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_nat_a_030700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_nat_a_030700

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_nat_a_030712_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_nat_a_030712
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_nat_a_030012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_nat_a_030012

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_nat_a_070300_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_nat_a_070300

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_nat_a_120003_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_nat_a_120003

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_nat_a_120703_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_nat_a_120703

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_nat_a_031207_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_nat_a_031207
  

##
# sim10k_mono_seq.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 1 Initial subpopns, expanding to 10
# output as genotypes for Powermarker
# Selection at locus 3 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 3, this occurs before the split also
# 
# phylip format
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 06, 2013 08:50:51 AM 
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_a_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
# 	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_a_000000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_a_030000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_a_030000

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_a_000700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
# 	        (0, 'Affected')], at = j1),
      sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_a_000700

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_a_000012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
# 	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_a_000012

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_a_030700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_a_030700

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_a_030712_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_a_030712
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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_a_030012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_a_030012

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_a_070300_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_a_070300

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_a_120003_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
#       sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_a_120003

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_a_120703_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_a_120703

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
    new_sz = [x//2000 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='phylip', output='sim10k_mono_seq_a_031207_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
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
	        (0, 'Affected')], at = j1),
        sim.MaSelector(loci=7, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci=12, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10k_mono_seq_a_031207

