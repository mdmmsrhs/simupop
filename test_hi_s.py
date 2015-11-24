import simuOpt
import simuPOP as sim
import math
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

## Settings ###################################################################################
n = 3   # Number of Subpopns
i = 100000 # Number of Indivs/Subpopn
l = 10    # Number of Loci per Chromosome
c = 5 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 500 # Total Number of Steps (Generations)
u = 0.00005 # Forward Mutation Rate
v = 0.000005# Backward Mutation Rate
m = 0.000005  # Migration rate
e = 98 # Proportion of selfing
s1 = 0.1 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 4Mb/cM) ie.(1/4/100)

##sim10g_hi_s_a

##################################################################################################

j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,20,2))*5, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
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
    export(sample, format='fstat', output='sim10g_hi_s_a_000000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
# 		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#         sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
		sim.Stat(pop, alleleFreq=list(range(9)), step = 50),
		sim.PyEval('gen', step=100, output = '>>sim10g_0_output.txt'),
		
        sim.PyEval("'\t%.3f' % alleleFreq[0][1]", step = 50, output = '>>sim10g_0_output.txt'),
        
        sim.PyEval("'\t%.3f' % alleleFreq[1][1]", step = 50, output = '>>sim10g_0_output.txt'),
        
        sim.PyEval("'\t%.3f' % alleleFreq[2][1]", step = 50, output = '>>sim10g_0_output.txt'),
     
        sim.PyEval("'\t%.3f' % alleleFreq[3][1]", step=50, output = '>>sim10g_0_output.txt'),
        
        sim.PyEval("'\t%.3f' % alleleFreq[4][1]", step=50, output = '>>sim10g_0_output.txt'),
                
        sim.PyOutput('\n', step=50, output = '>>sim10g_0_output.txt'),	
        
	sim.PyOperator(func=sampleAndExport, at=[0,500]),
		
    ],
	gen=(t+1),
)
print 'all done'
##sim10g_hi_s_a_000000    

##################################################################################################

j1 = 100   # Number of Steps (Generations) before selection event 1
j2 = 200   # Number of Steps (Generations) before selection event 2
j3 = 300   # Number of Steps (Generations) before selection event 3

d1 = 50 # duration of selection event 1 (in generations)
d2 = 50 # duration of selection event 2 (in generations)
d3 = 50 # duration of selection event 3 (in generations)

##################################################################################################

pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,20,2))*5, infoFields=['migrate_to','fitness'])
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter()
    ])
)
z=20
def demo(pop,gen):
	global q, g
	rate=2**15
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
    export(sample, format='fstat', output='sim10g_hi_s_a_040000_sample_%d.dat' % pop.dvars().gen, gui = False)
    return True

simu = sim.Simulator(pop, rep=1)
simu.evolve(
    initOps=[
	    sim.InitGenotype(freq=(0.5, 0.5)),
	],
	preOps=[
		sim.SNPMutator(u = u, v = v),
		sim.Migrator(rate=migrIslandRates(m, n)),
		sim.MapSelector(loci=4, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
#         sim.MapSelector(loci=14, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j2, end = (j2+d2)),# Positive Selection for allele 1 at locus given
#         sim.MapSelector(loci=34, fitness = {(0,0):1-s1, (0,1):1-s1, (1,1):1}, begin = j3, end = (j3+d3)),# Positive Selection for allele 1 at locus given
    ],
	
    #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
    matingScheme=sim.HeteroMating([
        sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
        sim.HomoMating(
            chooser=sim.CombinedParentsChooser(
                sim.RandomParentChooser(),
                sim.RandomParentChooser()),
            generator=sim.OffspringGenerator(
                ops= [
	                sim.MendelianGenoTransmitter(),
	                sim.Recombinator(intensity=r1)
                ],
            ),
            weight=100-e)
        ],
    subPopSize=demo
    ),
	postOps=[
	sim.Stat(pop, alleleFreq=list(range(9)), step = 50),
		sim.PyEval('gen', step=100, output = '>>sim10g_4_output.txt'),
		
        sim.PyEval("'\t%.3f' % alleleFreq[0][1]", step = 50, output = '>>sim10g_4_output.txt'),
        
        sim.PyEval("'\t%.3f' % alleleFreq[1][1]", step = 50, output = '>>sim10g_4_output.txt'),
        
        sim.PyEval("'\t%.3f' % alleleFreq[2][1]", step = 50, output = '>>sim10g_4_output.txt'),
     
        sim.PyEval("'\t%.3f' % alleleFreq[3][1]", step=50, output = '>>sim10g_4_output.txt'),
        
        sim.PyEval("'\t%.3f' % alleleFreq[4][1]", step=50, output = '>>sim10g_4_output.txt'),
                
        sim.PyOutput('\n', step=50, output = '>>sim10g_4_output.txt'),	
        
	sim.PyOperator(func=sampleAndExport, at=[0,500]),
		
    ],
	gen=(t+1),
)
print 'all done'
##sim10g_hi_s_a_040000   