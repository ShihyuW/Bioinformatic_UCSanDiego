from shihyu.Bioinformatics1 import PatternCount

def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(ExtendedGenome[i:i+(n//2)],symbol)
    return array

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    Extendedgenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2])
    for i in range(1, n):
        array[i] = array [i-1]
        if Extendedgenome[i-1] == symbol:
            array[i]-=1
        if Extendedgenome[i+(n//2)-1] == symbol:
            array[i]+=1                                                  
    return array

def SkewArray(Genome):
    Skew = [0]
    for i in range (len(Genome)):
        if Genome[i] == "A" or Genome[i] =="T":
            Skew.append(Skew[i]+0)
        elif Genome[i] == "C":
            Skew.append(Skew[i]-1)
        elif Genome[i] == "G":
            Skew.append(Skew[i]+1)
    return Skew

def MinimumSkew(Genome):
    positions = [] 
    M = SkewArray(Genome)
    for i in range(len(Genome)):
        if M[i] == min(M):
            positions.append(i)
   
    return positions

def HammingDistance(p, q):
    distance = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            distance +=1
    return distance

def ApproximatePatternMatching(Text, Pattern, D):
    positions = [] 
    n = len(Text)
    m = len(Pattern)
    for i in range(n-m+1):
        q = Text[i:i+m]
        d = HammingDistance(Pattern, D)
        if d <= D:        
            positions.append(i)
    
    return positions

def ApproximatePatternCount(Pattern, Text, D):
    count = 0 
    n = len(Text)
    m = len(Pattern)
    for i in range(n-m+1):
        q = Text[i:i+m]
        d = HammingDistance(Pattern, D)
        if d <= D:        
            count += 1
    return count