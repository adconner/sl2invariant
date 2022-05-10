
V = 8
L = ((0,4),(1,5),(2,6),(3,7))


#V = 12
#L = ((0,6),(1,7),(2,8),(3,9),(4,10),(5,11))








#####################################################


import sys
import itertools





def is_intersecting(tup0,tup1):
#this is very slow, but probably good enough for the sizes in our project (V = 18 at most)
  pointer = tup0[0]
  occ = []
  hits = [tup0[0],tup0[1],tup1[0],tup1[1]]
  while True:
    if pointer in hits:
      occ.append(pointer)
    pointer += 1
    if pointer == V:
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

def sorted_correctly(tup0,tup1):
  if (tup0[0]<tup1[0] and tup0[0]<tup1[1]) or (tup0[1]<tup1[0] and tup0[1]<tup1[1]):
    return True
  return False

def change_at_pos(L,i,t):
  ret = []
  for j in range(len(L)):
    if i != j:
      ret.append(L[j])
    else:
      ret.append(t)
  return tuple(ret)


def expansion_step(L,scalar):
  for i in range(V//2):
    if not order_correct(L[i]):
      newtup = (L[i][1],L[i][0])
      newL = change_at_pos(L,i,newtup)
      return {newL:-scalar}
  for pair in itertools.combinations(range(V//2),2):
    if not sorted_correctly(L[pair[0]],L[pair[1]]):
      newtup0=L[pair[1]]
      newtup1=L[pair[0]]
      newL = L
      newL = change_at_pos(newL,pair[0],newtup0)
      newL = change_at_pos(newL,pair[1],newtup1)
      return {newL:scalar}
  tup_index0 = -1
  tup_index1 = -1
  for pair in itertools.combinations(range(V//2),2):
    if is_intersecting(L[pair[0]],L[pair[1]]):
      tup_index0 = pair[0]
      tup_index1 = pair[1]
  if tup_index0 == -1:
    return {L:scalar}
  else:
    newtup0=(L[tup_index0][0],L[tup_index1][1])
    newtup1=(L[tup_index1][0],L[tup_index0][1])
    newLA = L
    newLA = change_at_pos(newLA,tup_index0,newtup0)
    newLA = change_at_pos(newLA,tup_index1,newtup1)
    newtup0=(L[tup_index0][0],L[tup_index1][0])
    newtup1=(L[tup_index1][1],L[tup_index0][1])
    newLB = L
    newLB = change_at_pos(newLB,tup_index0,newtup0)
    newLB = change_at_pos(newLB,tup_index1,newtup1)
    return {newLA:scalar,newLB:-scalar}

def roman(i):
  return ["-","I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII"][i]

def tikz_single(L):
  s = "\\pic"+roman(V)+"{"
  for i in range(V//2):
    s += "\\ar{"
    s += str(L[i][0])+"}{"+str(L[i][1])
    s += "}"
  s += "}"
  return s

def tikz_lincomb(dic):
  s = ""
  for k in dic:
    s+="\\rb{"
    if dic[k]>0: s+="+"
    s += str(dic[k]) +"}" + tikz_single(k) + "\n"
  return s[:-1] #remove the last newline

def sum_up_lincombs(dic0,dic1):
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


def dist(x,y):
  return min(abs(x-y), abs(V+x-y), abs(x-y-V), abs(V+y-x), abs(y-x-V))

def counterclockwise_half_starting_at(x):
  ret = []
  for i in range(V//2):
    ret.append((x+i)%V)
  return ret

def order_correct(tup):
  if dist(tup[0],tup[1])<V/2:
    if not tup[1] in counterclockwise_half_starting_at(tup[0]):
      return False
  else:
    if not tup[1] in counterclockwise_half_starting_at(0):
      return False
  return True
    


def done_single(L):
  for i in range(V//2):
    if not order_correct(L[i]):
      return False
  for pair in itertools.combinations(range(V//2),2):
    if not sorted_correctly(L[pair[0]],L[pair[1]]):
      return False
    if is_intersecting(L[pair[0]],L[pair[1]]):
      return False
  return True
  


print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print(tikz_single(L))
print("\\rb{=}")

dic = {L:1}
while True:
  #print(tikz_lincomb(dic))
  #print("\\\\\\rb{=}")
  done = True
  for L in dic.keys():
    if not done_single(L):
      oldcoeff = dic[L]
      found = {L:-oldcoeff}
      dic = sum_up_lincombs(dic,found)
      dic = sum_up_lincombs(dic,expansion_step(L,oldcoeff))
      done = False
      break
  if done:
    break

print(tikz_lincomb(dic))
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

