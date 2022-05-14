
# V = 8 
# L = ((0,4),(1,5),(2,6),(3,7))


#V = 12
#L = ((0,6),(1,7),(2,8),(3,9),(4,10),(5,11))

#####################################################


import sys
import itertools



class Basis():
    def __init__(self,V):
        self.V = V

    def is_intersecting(self,tup0,tup1):
    #this is very slow, but probably good enough for the sizes in our project (V = 18 at most)
      pointer = tup0[0]
      occ = []
      hits = [tup0[0],tup0[1],tup1[0],tup1[1]]
      while True:
        if pointer in hits:
          occ.append(pointer)
        pointer += 1
        if pointer == self.V:
          pointer = 0
        if pointer == occ[0]:
          break
      #there should be 4 entries in "occ" now
      if len(occ)!=4:
        print("ERROR!")
        sys.exit(1)
      if (occ[0] in tup0 and occ[1] in tup0) or (occ[1] in tup0 and occ[2] in tup0) or (occ[2] in tup0 and occ[3] in tup0) or (occ[3] in tup0 and occ[0] in tup0):
        return False
      return True

    def sorted_correctly(self,tup0,tup1):
      if (tup0[0]<tup1[0] and tup0[0]<tup1[1]) or (tup0[1]<tup1[0] and tup0[1]<tup1[1]):
        return True
      return False

    def change_at_pos(self,L,i,t):
      ret = []
      for j in range(len(L)):
        if i != j:
          ret.append(L[j])
        else:
          ret.append(t)
      return tuple(ret)


    def expansion_step(self,L,scalar):
      for i in range(self.V//2):
        if not self.order_correct(L[i]):
          newtup = (L[i][1],L[i][0])
          newL = self.change_at_pos(L,i,newtup)
          return {newL:-scalar}
      for pair in itertools.combinations(range(self.V//2),2):
        if not self.sorted_correctly(L[pair[0]],L[pair[1]]):
          newtup0=L[pair[1]]
          newtup1=L[pair[0]]
          newL = L
          newL = self.change_at_pos(newL,pair[0],newtup0)
          newL = self.change_at_pos(newL,pair[1],newtup1)
          return {newL:scalar}
      tup_index0 = -1
      tup_index1 = -1
      for pair in itertools.combinations(range(self.V//2),2):
        if self.is_intersecting(L[pair[0]],L[pair[1]]):
          tup_index0 = pair[0]
          tup_index1 = pair[1]
      if tup_index0 == -1:
        return {L:scalar}
      else:
        newtup0=(L[tup_index0][0],L[tup_index1][1])
        newtup1=(L[tup_index1][0],L[tup_index0][1])
        newLA = L
        newLA = self.change_at_pos(newLA,tup_index0,newtup0)
        newLA = self.change_at_pos(newLA,tup_index1,newtup1)
        newtup0=(L[tup_index0][0],L[tup_index1][0])
        newtup1=(L[tup_index1][1],L[tup_index0][1])
        newLB = L
        newLB = self.change_at_pos(newLB,tup_index0,newtup0)
        newLB = self.change_at_pos(newLB,tup_index1,newtup1)
        return {newLA:scalar,newLB:-scalar}

    def roman(self,i):
      return ["-","I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII"][i]

    def tikz_single(self,L):
      s = "\\pic"+self.roman(self.V)+"{"
      for i in range(self.V//2):
        s += "\\ar{"
        s += str(L[i][0])+"}{"+str(L[i][1])
        s += "}"
      s += "}"
      return s

    def tikz_lincomb(self,dic):
      s = ""
      for k in dic:
        s+="\\rb{"
        if dic[k]>0: s+="+"
        s += str(dic[k]) +"}" + self.tikz_single(k) + "\n"
      return s[:-1] #remove the last newline

    def sum_up_lincombs(self,dic0,dic1):
      if len(dic0)==0:
        return dic1
      if len(dic1)==0:
        return dic0
      ret = dict()
      for k in set.union(set(dic0.keys()),set(dic1.keys())):
        if k not in dic0:
          ret[k] = dic1[k]
        elif k not in dic1:
          ret[k] = dic0[k]
        else:
          ret[k] = dic0[k] + dic1[k]
        if ret[k]==0:
          ret.pop(k)
      return ret


    def dist(self,x,y):
      return min(abs(x-y), abs(self.V+x-y), abs(x-y-self.V), abs(self.V+y-x),
              abs(y-x-self.V))

    def counterclockwise_half_starting_at(self,x):
      ret = []
      for i in range(self.V//2):
        ret.append((x+i)%self.V)
      return ret

    def order_correct(self,tup):
      if self.dist(tup[0],tup[1])<self.V/2:
        if not tup[1] in self.counterclockwise_half_starting_at(tup[0]):
          return False
      else:
        if not tup[1] in self.counterclockwise_half_starting_at(0):
          return False
      return True
        


    def done_single(self,L):
      for i in range(self.V//2):
        if not self.order_correct(L[i]):
          return False
      for pair in itertools.combinations(range(self.V//2),2):
        if not self.sorted_correctly(L[pair[0]],L[pair[1]]):
          return False
        if self.is_intersecting(L[pair[0]],L[pair[1]]):
          return False
      return True
      
    def expand_uncrossing(self,L):
        L = tuple([tuple(ix) for ix in L])
        cs = {L:1}
        while True:
            done = True
            for L in cs.keys():
                if not self.done_single(L):
                    oldcoeff = cs[L]
                    found = {L:-oldcoeff}
                    cs = self.sum_up_lincombs(cs,found)
                    cs = self.sum_up_lincombs(cs,self.expansion_step(L,oldcoeff))
                    done = False
                    break
            if done:
                break
        cs_normalized = {}
        for P,e in cs.items():
            Pnorm = []
            enorm = e
            for (i,j) in P:
                if i < j:
                    Pnorm.append((i,j))
                else:
                    Pnorm.append((j,i))
                    enorm = -enorm
            Pnorm.sort()
            cs_normalized[tuple(Pnorm)] = enorm
        return cs_normalized


    def print_latex(self,L):
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        print(self.tikz_single(L))
        print("\\rb{=}")

        cs = self.expand_uncrossing(L)

        print(self.tikz_lincomb(cs))
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

# Basis(V).print_latex(L)
