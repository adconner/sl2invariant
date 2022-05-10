# n Length(hs)
# 2 1
# 3 0
# 4 4
# 5 0
# 6 13
# 7 0
# 8 40
# 9 0
# 10 117
# 12 525

n := 6;
g := SymmetricGroup(n);

chi := ClassFunction(g, List(List(ConjugacyClasses(g),Representative),r->
            2^(n-NrMovedPoints(r)+Sum(SortedList(CycleStructurePerm(r))))));

hs := [];
for h in ConjugacyClassesSubgroups(g) do
    h := Representative(h);
    Print(h,"\n");
    for psi in Irr(h) do
        if psi[1] = 1 then
            mult := ScalarProduct(RestrictedClassFunction(chi,h),psi);
            Print(psi," ",mult,"\n");
            if mult = 1 then
                Add(hs,psi); # psi remembers h: h = UnderlyingGroup(psi)
            fi;
        fi;
    od;
    Print(h,"\n");
od;

Print(Length(hs),"\n");

# vim: ft=python
