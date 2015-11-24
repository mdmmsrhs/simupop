##
# sim10l_mono_1of2_nat.py
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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_a_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_a_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_a_000000

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_a_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_a_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_a_040000

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_a_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_a_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_a_001400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_a_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_a_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_a_000034

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_a_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_a_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_a_041400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_a_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_a_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_a_041434
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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_a_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_a_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_a_040034

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_a_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_a_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_a_140400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_a_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_a_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_a_340004

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_a_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_a_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_a_341404

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_a_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_a_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
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

#
# sim10l_mono_1of2_seq_nat.py
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 10
# output as genotypes for Powermarker
# art selection in 1 of 3 original popns
# 
# phylip format
# 
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 11, 2013 11:49:36 AM 
# #


import simuOpt
import simuPOP as sim
import math
import os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 3    # Number of Loci per Chromosome
c = 1 # Number of Chromosomes
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_a_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_a_000000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_a_000000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_a_000000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_a_030000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_a_030000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_030000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_a_030000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_030000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_a_030000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_a_000700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_a_000700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_000700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_a_000700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_000700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_a_000700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_a_000012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_a_000012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_000012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_a_000012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_000012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_a_000012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_a_030700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_a_030700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_030700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_a_030700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_030700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_a_030700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_a_030712_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_a_030712_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_030712_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_a_030712_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_030712_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_a_030712
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_a_030012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_a_030012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_030012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_a_030012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_030012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_a_030012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_a_070300_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_a_070300_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_070300_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_a_070300_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_070300_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_a_070300

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_a_120003_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_a_120003_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_120003_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_a_120003_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_120003_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_a_120003

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_a_120703_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_a_120703_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_120703_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_a_120703_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_120703_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_a_120703

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_a_031207_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_a_031707_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_031207_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_a_031207_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_a_031207_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_a_031207
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 10
# output as genotypes for Powermarker
# art selection in 1 of 3 original popns
# 
# phylip format
# 
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 11, 2013 11:49:36 AM 
# #


import simuOpt
import simuPOP as sim
import math
import os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
# Settings ###################################################################################
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_a_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_a_000000_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_a_000000_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_a_000000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_a_030000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_a_030000_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_030000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_a_030000_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_030000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_a_030000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_a_000700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_a_000700_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_000700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_a_000700_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_000700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_a_000700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_a_000012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_a_000012_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_000012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_a_000012_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_000012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_a_000012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_a_030700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_a_030700_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_030700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_a_030700_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_030700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_a_030700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_a_030712_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_a_030712_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_030712_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_a_030712_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_030712_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_a_030712
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_a_030012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_a_030012_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_030012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_a_030012_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_030012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_a_030012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_a_070300_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_a_070300_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_070300_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_a_070300_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_070300_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_a_070300

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_a_120003_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_a_120003_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_120003_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_a_120003_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_120003_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_a_120003

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_a_120703_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_a_120703_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_120703_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_a_120703_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_120703_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_a_120703

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_a_031207_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_a_031207_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_031207_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_a_031207_sample_%d.phy    sim10l_mono_1of2_seq_nat_a_031207_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_a_031207

##
# sim10l_mono_1of2_nat.py
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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_b_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_b_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_b_000000

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_b_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_b_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_b_040000

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_b_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_b_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_b_001400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_b_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_b_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_b_000034

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_b_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_b_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_b_041400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_b_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_b_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_b_041434
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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_b_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_b_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_b_040034

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_b_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_b_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_b_140400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_b_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_b_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_b_340004

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_b_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_b_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_b_341404

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_b_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_b_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
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

#
# sim10l_mono_1of2_seq_nat.py
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 10
# output as genotypes for Powermarker
# art selection in 1 of 3 original popns
# 
# phylip format
# 
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 11, 2013 11:49:36 AM 
# #


import simuOpt
import simuPOP as sim
import math
import os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 3    # Number of Loci per Chromosome
c = 1 # Number of Chromosomes
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_b_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_b_000000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_b_000000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_b_000000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_b_030000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_b_030000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_030000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_b_030000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_030000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_b_030000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_b_000700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_b_000700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_000700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_b_000700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_000700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_b_000700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_b_000012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_b_000012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_000012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_b_000012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_000012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_b_000012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_b_030700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_b_030700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_030700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_b_030700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_030700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_b_030700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_b_030712_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_b_030712_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_030712_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_b_030712_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_030712_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_b_030712
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_b_030012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_b_030012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_030012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_b_030012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_030012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_b_030012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_b_070300_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_b_070300_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_070300_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_b_070300_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_070300_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_b_070300

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_b_120003_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_b_120003_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_120003_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_b_120003_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_120003_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_b_120003

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_b_120703_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_b_120703_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_120703_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_b_120703_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_120703_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_b_120703

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_b_031207_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_b_031707_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_031207_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_b_031207_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_b_031207_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_b_031207
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 10
# output as genotypes for Powermarker
# art selection in 1 of 3 original popns
# 
# phylip format
# 
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 11, 2013 11:49:36 AM 
# #


import simuOpt
import simuPOP as sim
import math
import os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
# Settings ###################################################################################
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_b_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_b_000000_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_b_000000_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_b_000000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_b_030000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_b_030000_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_030000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_b_030000_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_030000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_b_030000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_b_000700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_b_000700_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_000700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_b_000700_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_000700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_b_000700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_b_000012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_b_000012_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_000012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_b_000012_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_000012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_b_000012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_b_030700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_b_030700_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_030700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_b_030700_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_030700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_b_030700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_b_030712_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_b_030712_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_030712_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_b_030712_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_030712_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_b_030712
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_b_030012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_b_030012_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_030012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_b_030012_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_030012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_b_030012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_b_070300_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_b_070300_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_070300_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_b_070300_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_070300_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_b_070300

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_b_120003_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_b_120003_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_120003_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_b_120003_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_120003_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_b_120003

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_b_120703_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_b_120703_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_120703_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_b_120703_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_120703_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_b_120703

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_b_031207_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_b_031207_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_031207_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_b_031207_sample_%d.phy    sim10l_mono_1of2_seq_nat_b_031207_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_b_031207

##
# sim10l_mono_1of2_nat.py
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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_c_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_c_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_c_000000

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_c_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_c_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_c_040000

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_c_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_c_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_c_001400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_c_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_c_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_c_000034

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_c_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_c_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_c_041400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_c_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_c_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_c_041434
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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_c_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_c_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_c_040034

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_c_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_c_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_c_140400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_c_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_c_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_c_340004

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_c_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_c_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_c_341404

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_c_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_c_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
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

#
# sim10l_mono_1of2_seq_nat.py

# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 10
# output as genotypes for Powermarker
# art selection in 1 of 3 original popns
# 
# phylip format
# 
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 11, 2013 11:49:36 AM 
# #


import simuOpt
import simuPOP as sim
import math
import os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 3    # Number of Loci per Chromosome
c = 1 # Number of Chromosomes
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_c_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_c_000000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_c_000000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_c_000000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_c_030000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_c_030000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_030000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_c_030000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_030000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_c_030000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_c_000700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_c_000700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_000700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_c_000700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_000700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_c_000700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_c_000012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_c_000012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_000012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_c_000012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_000012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_c_000012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_c_030700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_c_030700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_030700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_c_030700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_030700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_c_030700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_c_030712_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_c_030712_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_030712_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_c_030712_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_030712_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_c_030712
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_c_030012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_c_030012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_030012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_c_030012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_030012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_c_030012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_c_070300_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_c_070300_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_070300_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_c_070300_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_070300_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_c_070300

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_c_120003_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_c_120003_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_120003_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_c_120003_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_120003_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_c_120003

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_c_120703_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_c_120703_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_120703_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_c_120703_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_120703_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_c_120703

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_c_031207_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_c_031707_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_031207_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_c_031207_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_c_031207_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_c_031207
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 10
# output as genotypes for Powermarker
# art selection in 1 of 3 original popns
# 
# phylip format
# 
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 11, 2013 11:49:36 AM 
# #


import simuOpt
import simuPOP as sim
import math
import os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
# Settings ###################################################################################
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_c_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_c_000000_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_c_000000_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_c_000000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_c_030000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_c_030000_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_030000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_c_030000_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_030000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_c_030000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_c_000700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_c_000700_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_000700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_c_000700_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_000700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_c_000700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_c_000012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_c_000012_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_000012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_c_000012_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_000012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_c_000012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_c_030700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_c_030700_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_030700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_c_030700_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_030700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_c_030700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_c_030712_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_c_030712_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_030712_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_c_030712_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_030712_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_c_030712
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_c_030012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_c_030012_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_030012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_c_030012_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_030012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_c_030012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_c_070300_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_c_070300_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_070300_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_c_070300_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_070300_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_c_070300

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_c_120003_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_c_120003_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_120003_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_c_120003_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_120003_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_c_120003

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_c_120703_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_c_120703_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_120703_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_c_120703_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_120703_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_c_120703

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_c_031207_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_c_031207_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_031207_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_c_031207_sample_%d.phy    sim10l_mono_1of2_seq_nat_c_031207_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_c_031207

##
# sim10l_mono_1of2_nat.py
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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_d_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_d_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_d_000000

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_d_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_d_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_d_040000

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_d_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_d_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_d_001400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_d_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_d_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_d_000034

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_d_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_d_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_d_041400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_d_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_d_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_d_041434
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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_d_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_d_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_d_040034

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_d_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_d_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_d_140400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_d_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_d_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_d_340004

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_d_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_d_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_d_341404

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_d_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_d_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
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

#
# sim10l_mono_1of2_seq_nat.py
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 10
# output as genotypes for Powermarker
# art selection in 1 of 3 original popns
# 
# phylip format
# 
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 11, 2013 11:49:36 AM 
# #


import simuOpt
import simuPOP as sim
import math
import os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 3    # Number of Loci per Chromosome
c = 1 # Number of Chromosomes
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_d_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_d_000000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_d_000000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_d_000000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_d_030000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_d_030000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_030000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_d_030000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_030000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_d_030000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_d_000700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_d_000700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_000700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_d_000700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_000700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_d_000700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_d_000012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_d_000012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_000012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_d_000012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_000012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_d_000012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_d_030700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_d_030700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_030700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_d_030700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_030700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_d_030700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_d_030712_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_d_030712_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_030712_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_d_030712_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_030712_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_d_030712
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_d_030012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_d_030012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_030012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_d_030012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_030012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_d_030012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_d_070300_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_d_070300_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_070300_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_d_070300_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_070300_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_d_070300

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_d_120003_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_d_120003_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_120003_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_d_120003_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_120003_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_d_120003

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_d_120703_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_d_120703_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_120703_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_d_120703_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_120703_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_d_120703

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_d_031207_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_d_031707_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_031207_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_d_031207_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_d_031207_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_d_031207

# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 10
# output as genotypes for Powermarker
# art selection in 1 of 3 original popns
# 
# phylip format
# 
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 11, 2013 11:49:36 AM 
# #


import simuOpt
import simuPOP as sim
import math
import os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
# Settings ###################################################################################
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_d_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_d_000000_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_d_000000_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_d_000000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_d_030000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_d_030000_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_030000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_d_030000_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_030000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_d_030000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_d_000700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_d_000700_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_000700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_d_000700_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_000700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_d_000700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_d_000012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_d_000012_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_000012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_d_000012_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_000012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_d_000012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_d_030700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_d_030700_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_030700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_d_030700_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_030700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_d_030700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_d_030712_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_d_030712_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_030712_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_d_030712_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_030712_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_d_030712
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_d_030012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_d_030012_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_030012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_d_030012_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_030012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_d_030012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_d_070300_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_d_070300_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_070300_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_d_070300_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_070300_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_d_070300

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_d_120003_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_d_120003_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_120003_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_d_120003_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_120003_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_d_120003

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_d_120703_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_d_120703_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_120703_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_d_120703_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_120703_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_d_120703

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_d_031207_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_d_031207_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_031207_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_d_031207_sample_%d.phy    sim10l_mono_1of2_seq_nat_d_031207_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_d_031207

##
# sim10l_mono_1of2_nat.py
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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_e_000000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_e_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_e_000000

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_e_040000_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_e_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_e_040000

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_e_001400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_e_001400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_e_001400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_e_000034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_e_000034_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_e_000034

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_e_041400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_e_041400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_e_041400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_e_041434_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_e_041434_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_e_041434
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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_e_040034_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_e_040034_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_e_040034

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_e_140400_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_e_140400_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_e_140400

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_e_340004_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_e_340004_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_e_340004

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_e_341404_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_e_341404_sample_%d.dat' % pop.dvars().gen, gui = False)
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
##sim10l_mono_1of2_nat_e_341404

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
    new_sz = [x//200 for x in sz]
    sample = drawRandomSample(pop, new_sz)
    export(sample, format='genepop', output='sim10l_mono_1of2_nat_e_043414_sample_%d.gen' % pop.dvars().gen, adjust = 0, gui = False)
    export(sample, format='fstat', output='sim10l_mono_1of2_nat_e_043414_sample_%d.dat' % pop.dvars().gen, gui = False)
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

#
# sim10l_mono_1of2_seq_nat.py
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 10
# output as genotypes for Powermarker
# art selection in 1 of 3 original popns
# 
# phylip format
# 
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 11, 2013 11:49:36 AM 
# #


import simuOpt
import simuPOP as sim
import math
import os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 3    # Number of Loci per Chromosome
c = 1 # Number of Chromosomes
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_e_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_e_000000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_e_000000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_e_000000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_e_030000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_e_030000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_030000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_e_030000_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_030000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_e_030000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_e_000700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_e_000700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_000700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_e_000700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_000700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_e_000700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_e_000012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_e_000012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_000012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_e_000012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_000012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
# 		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
      sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_e_000012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_e_030700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_e_030700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_030700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_e_030700_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_030700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_e_030700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_e_030712_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_e_030712_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_030712_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_e_030712_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_030712_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_e_030712
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_e_030012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_e_030012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_030012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_e_030012_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_030012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_e_030012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_e_070300_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_e_070300_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_070300_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_e_070300_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_070300_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_e_070300

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_e_120003_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_e_120003_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_120003_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_e_120003_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_120003_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),

	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
#       sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_e_120003

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_e_120703_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_e_120703_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_120703_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_e_120703_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_120703_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_e_120703

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_subset_nat_e_031207_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_subset_nat_e_031707_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_031207_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_subset_nat_e_031207_sample_%d.phy    sim10l_mono_1of2_seq_subset_nat_e_031207_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(loci = 0, freq=(0.0, 0.80, 0.20, 0.0)),
	    sim.InitGenotype(loci = 1, freq=(0.2, 0.0, 0.8, 0.0)),
	    sim.InitGenotype(loci = 2, freq=(0.0, 0.2, 0.0, 0.8)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
		sim.MaSelector(loci = 0, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1), subPops = 0),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 1, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
        sim.MaSelector(loci = 2, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.PyOperator(func=sampleAndExport, at = 500),
    ],
	gen=(t+1),
)
print 'all done'
##sim10l_mono_1of2_seq_subset_nat_e_031207
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# 2 Initial subpopns, expanding to 10
# output as genotypes for Powermarker
# art selection in 1 of 3 original popns
# 
# phylip format
# 
# 
# #
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 11, 2013 11:49:36 AM 
# #


import simuOpt
import simuPOP as sim
import math
import os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
# Settings ###################################################################################
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_e_000000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_e_000000_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_000000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_e_000000_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_000000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_e_000000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_e_030000_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_e_030000_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_030000_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_e_030000_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_030000_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_e_030000

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_e_000700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_e_000700_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_000700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_e_000700_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_000700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_e_000700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_e_000012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_e_000012_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_000012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_e_000012_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_000012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_e_000012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_e_030700_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_e_030700_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_030700_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_e_030700_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_030700_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_e_030700

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_e_030712_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_e_030712_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_030712_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_e_030712_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_030712_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_e_030712
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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_e_030012_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_e_030012_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_030012_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_e_030012_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_030012_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_e_030012

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_e_070300_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_e_070300_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_070300_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_e_070300_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_070300_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_e_070300

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_e_120003_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_e_120003_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_120003_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_e_120003_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_120003_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_e_120003

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_e_120703_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_e_120703_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_120703_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_e_120703_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_120703_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_e_120703

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
    export(sample, format='phylip', output='sim10l_mono_1of2_seq_nat_e_031207_sample_%d.phy' % pop.dvars().gen, alleleNames = ('A','C','G','T'), gui=False),
    os.system('perl convert_diploid.pl   N   sim10l_mono_1of2_seq_nat_e_031207_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_031207_merged_sample_%d.phy' % (pop.dvars().gen, pop.dvars().gen)),
    os.system('perl phylip_to_fasta.pl   sim10l_mono_1of2_seq_nat_e_031207_sample_%d.phy    sim10l_mono_1of2_seq_nat_e_031207_sample_%d.fas' % (pop.dvars().gen, pop.dvars().gen)),
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
##sim10l_mono_1of2_seq_nat_e_031207

