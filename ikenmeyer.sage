n = 6
g = gap.SymmetricGroup([1..n])
load('all.sage')

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

# psis = []
# for h in g.ConjugacyClassesSubgroups():
#     h = h.Representative()
#     print(h)
#     for psi in h.LinearCharacters():
#         print (psi)
#         mult = gap.ScalarProduct(gap.RestrictedClassFunction(chi,h),psi)
#         print (mult)
#         if mult == 1:
#             psis.append(psi)
#             print ('found')

# lat = g.LatticeSubgroups()
# minsup = lat.MinimalSupergroupsLattice().sage()
# lat_graph = DiGraph()
# for i,above in enumerate(minsup):
#     for j,_ in above:
#         lat_graph.add_edge(i+1,j)

# psis = []
# for i in lat_graph.topological_sort()[::-1]:
#     h = lat.ConjugacyClassesSubgroups()[i].Representative()

psis = []
psisl = []
for h in sorted(g.ConjugacyClassesSubgroups(),key=lambda h: -h.Representative().Size()):
    h = h.Representative()
    print 'hsize',h.Size()
    e_already = {} # just to avoid redundant computation
    for psio in gap.OrbitsDomain(g.Normalizer(h),h.LinearCharacters()):
        psi = psio[1] 
        # check one linear character in each conjugacy class of the normalizer
        print psi
        mult = gap.ScalarProduct(gap.RestrictedClassFunction(chi,h),psi)
        if mult != 1:
            continue
        print mult
        already = False
        for psi2i,psi2 in enumerate(psis):
            k = psi2.UnderlyingGroup()
            if not psi2.IsSubsetSet(psi):
                continue
            if psi2i not in e_already:
                e_already[psi2i] = gap.ContainedConjugates(g,k,h)
            for _,x in e_already[psi2i]:
                for psiother in psio:
                    # pure gap for performance reasons
                    if gap('ForAll(List(ConjugacyClasses(%s),Representative),clr->(clr^%s)^%s = clr^%s)' %\
                                    (h.name(),x.name(),psi2.name(),psiother.name())):
                    # # equivalent to 
                    # if all((clr^x)^psi2 == clr^psiother for cl in
                    #         h.ConjugacyClasses() for clr in [cl.Representative()]):
                        psisl[psi2i] = (psi,x)
                        already = True
                        break
                if already:
                    break
            if already:
                break
        if not already:
            psis.append(psi)
            psisl.append(psi)
            print 'psis',len(psis)


# next two functions should be identical, second inlines some functions and
# works with sparse matrices
def psi_projection(psi):
    M = matrix(UniversalCyclotomicField(),
        catalan_number(n//2),catalan_number(n//2))
    for e,c in zip(psi,psi.UnderlyingCharacterTable().ConjugacyClasses()):
        for g in c.List():
            M += e.sage()*matrix_of_perm(g,n)
    return M

def psi_projection_direct(psi,conj=None):
    M = {}
    for Li,L in enumerate(DyckWords(n//2)):
        L = dyck_word_to_pairing(L)
        for e,c in zip(psi,psi.UnderlyingCharacterTable().ConjugacyClasses()):
            for g in c.List():
                if conj is not None:
                    conj.Inverse()*g*conj
                Lperm = [(int((i+1)^g-1),int((j+1)^g-1)) for i,j in L]
                cs = Basis(n).expand_uncrossing(Lperm)
                for w,f in cs.items():
                    ix = (Li,DyckWords(n//2).rank(pairing_to_dyck_word(w)))
                    M[ix] = M.get(ix,0) + e*f
    return matrix(UniversalCyclotomicField(),
        catalan_number(n//2),catalan_number(n//2), M)


projs = []
for psi in psisl:
    print (psi)
    projs.append(psi_projection_direct(psi))

Ts = [p.row_space().basis_matrix().list() for p in projs]

rep = gap.GroupHomomorphismByImages(g,[matrix_of_perm(e,n)
    for e in g.GeneratorsOfGroup()])

gens = g.GeneratorsOfGroup()
acts = [e^rep for e in g.GeneratorsOfGroup()]

print("removing duplicates (Sn-conjugates)")
Ts_by_stab = {}
for T in Ts:
    stab = gap.Stabilizer(g,T,gens,acts,gap.OnLines)
    Ts_by_stab.setdefault(int(stab.Size()),[]).append((T,stab))
Ts_dedup = []
for _,THs in sorted(Ts_by_stab.items(),key=lambda p: -p[0]):
    for i,(T,stab) in enumerate(THs):
        occ_later = False
        for S,_ in THs[i+1:]:
            # IsBool <-> fail here
            if not g.RepresentativeAction(T,S,gens,acts,gap.OnLines).IsBool():
                occ_later = True
                break
        if not occ_later:
            Ts_dedup.append(T)

# # print("all orbits, very slow n>6")
# # Tsall = gap.Orbits(rep.Image(),Ts,gap.OnLines).sage()
# # stabs = [gap.Stabilizer(g,o[0],gens,acts,gap.OnLines) for o in Tsall]

# # M = matrix([T for o in Tsall for T in o])
# # M = M.change_ring(CyclotomicField(5*4*3))

# vim: ft=python
