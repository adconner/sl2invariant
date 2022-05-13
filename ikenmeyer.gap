SL2DeterminedByStabilizer := function (n)
    local g,hi,h,cc_cache,tbl,chiix,chi,psis,psisl,hs,psio,psi,mult,already,psi2i,psi2,k,p,x;
    g := SymmetricGroup(n);

    tbl := CharacterTable("symmetric",n);
    chiix := Position(CharacterParameters(tbl),[1,[n/2,n/2]]);
    # strictly speaking maybe we should use CompatibleConjugacyClasses to match up
    # the columns of tbl with those of the character table of g. I think gap
    # computes that for a known symmetric group this way anyway, so there is no
    # issue
    chi := Irr(g)[chiix];

    psis := [];
    psisl := [];
    hs := List(ConjugacyClassesSubgroups(g),Representative);
    SortBy(hs, h -> -Size(h));
    for hi in [1..Size(hs)] do
        h := hs[hi];
        Print(hi,"/",Size(hs),": size ",Size(h),", ",Length(psis)," found\n");
        cc_cache := [];
        for psio in OrbitsDomain(Normalizer(g,h),LinearCharacters(h)) do
            psi := psio[1];
            # Print(psi,"\n");
            mult := ScalarProduct(RestrictedClassFunction(chi,h),psi);
            if mult <> 1 then
                continue;
            fi;
            Print("found, checking if duplicate...\c");
            already := false;
            for psi2i in [1..Length(psis)] do
                psi2 := psis[psi2i];
                k := UnderlyingGroup(psi2);
                if not IsSubsetSet(psi2,psi) then
                    continue;
                fi;
                if not IsBound(cc_cache[psi2i]) then
                    Print("*\c");
                    cc_cache[psi2i] := ContainedConjugates(g,k,h);
                    Print(".\c");
                fi;
                for p in cc_cache[psi2i] do
                    x := p[2];
                    if ForAny(psio,psiother->
                            ForAll(List(ConjugacyClasses(h),Representative),
                            clr->(clr^x)^psi2 = clr^psiother)) then
                        psisl[psi2i] := [psi,x];
                        Print(" yes\n");
                        already := true;
                        break;
                    fi;

                od;
                if already then
                    break;
                fi;
            od;
            if not already then
                Add(psis,psi);
                Add(psisl,[psi,()]);
                Print(" no\n");
            fi;
        od;
    od;
    return [psis,psisl];
end;

# vim: ft=python
