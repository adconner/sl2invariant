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

gap.Read('\"sl2inv.gap\"');
psis, psisl = gap.SL2DeterminedByStabilizer(n)
# gap.Read('\"psisn10.gap\"');
# psis, psisl = gap('psisn10')

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

def psi_row_space_direct(psi):
    print 'group size',psi.UnderlyingGroup().Size()
    for L in DyckWords(n//2):
    # for L in sorted(DyckWords(n//2),key=lambda _: RDF.random_element()):
        print L
        L = dyck_word_to_pairing(L)
        vec = vector(UniversalCyclotomicField(),catalan_number(n//2))
        for e,c in zip(psi,psi.UnderlyingCharacterTable().ConjugacyClasses()):
            for g in c.List():
                Lperm = [(int((i+1)^g-1),int((j+1)^g-1)) for i,j in L]
                cs = Basis(n).expand_uncrossing(Lperm)
                for w,f in cs.items():
                    vec[DyckWords(n//2).rank(pairing_to_dyck_word(w))] += e*f
        if not vec.is_zero():
            vec /= next(e for e in vec if e != 0)
            return vec.list()

Ts = []
for i,(psi,x) in enumerate(psisl):
    print (i,psi)
    Ts.append(psi_row_space_direct(psi))



# projs = []
# for psi,x in psisl:
#     print (psi)
#     projs.append(psi_projection_direct(psi))

# Ts = [p.row_space().basis_matrix().list() for p in projs]

# rep = gap.GroupHomomorphismByImages(g,[matrix_of_perm(e,n)
#     for e in g.GeneratorsOfGroup()])

# gens = g.GeneratorsOfGroup()
# acts = [e^rep for e in g.GeneratorsOfGroup()]

# print("removing duplicates (Sn-conjugates)")
# Ts_by_stab = {}
# for T in Ts:
#     stab = gap.Stabilizer(g,T,gens,acts,gap.OnLines)
#     Ts_by_stab.setdefault(int(stab.Size()),[]).append((T,stab))
# Ts_dedup = []
# for _,THs in sorted(Ts_by_stab.items(),key=lambda p: -p[0]):
#     for i,(T,stab) in enumerate(THs):
#         occ_later = False
#         for S,_ in THs[i+1:]:
#             # IsBool <-> fail here
#             if not g.RepresentativeAction(T,S,gens,acts,gap.OnLines).IsBool():
#                 occ_later = True
#                 break
#         if not occ_later:
#             Ts_dedup.append(T)

# # print("all orbits, very slow n>6")
# # Tsall = gap.Orbits(rep.Image(),Ts,gap.OnLines).sage()
# # stabs = [gap.Stabilizer(g,o[0],gens,acts,gap.OnLines) for o in Tsall]

# # M = matrix([T for o in Tsall for T in o])
# # M = M.change_ring(CyclotomicField(5*4*3))

# vim: ft=python
