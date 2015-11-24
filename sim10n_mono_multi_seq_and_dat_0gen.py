##
# sim10n_mono.py
#
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural selection and sampling.  Selection occuring at different periods.
# # 1 Initial subpopns, expanding to 6
# output as genotypes for Powermarker
# Selection at locus 4 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 4, this occurs before the split also
# 
# genepop, fstat format
#
##
# Author: Richard Stephens
# Created: October 29, 2012  08:56:33 AM
# Modified: February 25, 2013 12:24:31 PM	  
##


import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 6   # Initial Number of Subpopns
d = 1 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 150 # Total Number of Steps (Generations)
u = 0.005 # Forward Mutation Rate
v = 0.0005# Backward Mutation Rate
m = 0.0005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)

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
    export(sample, format='fstat', output='sim10n_mono_a_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
		sim.PyOperator(func=sampleAndExport, at = 0)
    ],
	gen=(t+1),
)
print 'all done'
##sim10n_mono_a_000000

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
    export(sample, format='fstat', output='sim10n_mono_b_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
		sim.PyOperator(func=sampleAndExport, at = 0)
    ],
	gen=(t+1),
)
print 'all done'
##sim10n_mono_b_000000

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
    export(sample, format='fstat', output='sim10n_mono_c_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
		sim.PyOperator(func=sampleAndExport, at = 0)
    ],
	gen=(t+1),
)
print 'all done'
##sim10n_mono_c_000000

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
    export(sample, format='fstat', output='sim10n_mono_d_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
		sim.PyOperator(func=sampleAndExport, at = 0)
    ],
	gen=(t+1),
)
print 'all done'
##sim10n_mono_d_000000

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
    export(sample, format='fstat', output='sim10n_mono_e_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
		sim.PyOperator(func=sampleAndExport, at = 0)
    ],
	gen=(t+1),
)
print 'all done'
##sim10n_mono_e_000000

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
    export(sample, format='fstat', output='sim10n_mono_f_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
		sim.PyOperator(func=sampleAndExport, at = 0)
    ],
	gen=(t+1),
)
print 'all done'
##sim10n_mono_f_000000

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
    export(sample, format='fstat', output='sim10n_mono_g_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
		sim.PyOperator(func=sampleAndExport, at = 0)
    ],
	gen=(t+1),
)
print 'all done'
##sim10n_mono_g_000000

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
    export(sample, format='fstat', output='sim10n_mono_h_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
		sim.PyOperator(func=sampleAndExport, at = 0)
    ],
	gen=(t+1),
)
print 'all done'
##sim10n_mono_h_000000

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
    export(sample, format='fstat', output='sim10n_mono_i_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
		sim.PyOperator(func=sampleAndExport, at = 0)
    ],
	gen=(t+1),
)
print 'all done'
##sim10n_mono_i_000000

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
    export(sample, format='fstat', output='sim10n_mono_j_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
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
		sim.PyOperator(func=sampleAndExport, at = 0)
    ],
	gen=(t+1),
)
print 'all done'
##sim10n_mono_j_000000
