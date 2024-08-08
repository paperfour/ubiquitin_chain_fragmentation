from fragmentation import *
from parsing import *
import time
import tracemalloc

s = "substrate, 41(Ub,6(Ub),48(Ub)), 113(Ub,63(Ub,6(Ub)))"
#print(s)


#c = simpleParse(s, SimpleProtein(simpleUbSites(GFP_SEQ)))
#print(simpleUnparse(c))

t = time.time()
for i in range(0, 763000):
    c = getNestedListFromString(s, Protein(GFP_SEQ))
t1 = time.time()
print("In " + str(t1 - t))

print(unparse(getProteinCollection(c)))



print("=================")

raise Exception(":)")

x = []

t = time.time()

tracemalloc.start()


for i in range(0, 10000):
    c = simpleParse(s, SimpleProtein(simpleUbSites(GFP_SEQ)))
    x.append(simpleUnparse(c))

r = tracemalloc.get_traced_memory()
t1 = time.time()

print("In " + str(t1 - t) + " using " + str(r))
print("=================")

x = []

t = time.time()

tracemalloc.reset_peak()

for i in range(0, 10000):
    c = getNestedListFromString(s, Protein(GFP_SEQ))
    x.append(unparse(getProteinCollection(c)))

r = tracemalloc.get_traced_memory()
t1 = time.time()

print("In " + str(t1 - t) + " using " + str(r))
