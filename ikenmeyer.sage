n = 6
g = gap.SymmetricGroup([1..n])


# character of the SL2 invariants in (C2)^n, from schur weyl duality we obtain
# this from the two row and n/2 column representation
tbl = gap.CharacterTable('\"symmetric\"',n)
chiix = tbl.CharacterParameters().Position([1,[n//2,n//2]])
# strictly speaking maybe we should use CompatibleConjugacyClasses to match up
# the columns of tbl with those of the character table of g. I think gap
# computes that for a known symmetric group this way anyway, so there is no
# issue
chi = g.Irr()[chiix]

# character of action on (C2)^n
chifull = g.ClassFunction([
    2**(n-int(r.NrMovedPoints())+ r.CycleStructurePerm().SortedList().Sum())
            for c in g.CharacterTable().ConjugacyClasses() for r in [c.Representative()]])

psis = []
for h in g.ConjugacyClassesSubgroups():
    h = h.Representative()
    print(h)
    for psi in h.LinearCharacters():
        print (psi)
        mult = gap.ScalarProduct(gap.RestrictedClassFunction(chi,h),psi)
        print (mult)
        if mult == 1:
            psis.append(psi)
            print ('found')

import numpy as np
arr = np.zeros((2,2),dtype='object')
arr[0,1] = 1
arr[1,0] = -1

T = reduce(lambda a,b: np.tensordot(a,b,0),[arr]*(n//2))

Ss = []
for psi in psis:
    print (psi)
    # h = psi.UnderlyingGroup()
    S = np.zeros(T.shape,dtype='object')
    for e,c in zip(psi,psi.UnderlyingCharacterTable().ConjugacyClasses()):
        for g in c.List():
            S += e.sage()*T.transpose([int(i^g-1) for i in [1..n]])
    S /= psi.UnderlyingGroup().Size()
    Ss.append(S)





# vim: ft=python
