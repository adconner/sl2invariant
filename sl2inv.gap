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

CheckUnique := function(psis,n)
    local g,ok,psii,hi,h,cc_cache,tbl,chiix,chi,psisl,hs,psio,psi,mult,already,psi2i,psi2,k,p,x;
    g := SymmetricGroup(n);

    tbl := CharacterTable("symmetric",n);
    chiix := Position(CharacterParameters(tbl),[1,[n/2,n/2]]);
    # strictly speaking maybe we should use CompatibleConjugacyClasses to match up
    # the columns of tbl with those of the character table of g. I think gap
    # computes that for a known symmetric group this way anyway, so there is no
    # issue
    chi := Irr(g)[chiix];

    ok := true;

    for psii in [1..Length(psis)] do
        psi := psis[psii];
        h := UnderlyingGroup(psi);
        if ScalarProduct(RestrictedClassFunction(chi,h),psi) <> 1 then
            Print("fail 1 ",psii,"\n");
            ok := false;
        fi;

        for psi2i in [1..psii-1] do
            psi2 := psis[psi2i];
            k := UnderlyingGroup(psi2);
            for p in ContainedConjugates(g,h,k) do
                x := p[2];
                if ForAny(Orbit(Normalizer(g,h),psi),psiother->
                        ForAll(List(ConjugacyClasses(h),Representative),
                        clr->(clr^x)^psi2 = clr^psiother)) then
                    Print("fail 2 ",psii," ",psi2i,"\n");
                    ok := false;
                fi;

            od;
        od;
    od;
    return ok;
end;


# vim: ft=python
