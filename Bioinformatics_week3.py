def Count(Motifs):
    Count = {} # initializing the count dictionary
    k = len(Motifs[0])
    for s in "ACGT":
        Count[s] = []
        for j in range(k):
            Count[s].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            indi = Motifs[i][j]
            Count[indi][j]+=1
            
    return Count

def Profile(Motifs):
    t = len(Motifs)#number of seq
    k = len(Motifs[0])#lenth of each seq
    profile = {}
    c = Count(Motifs)#the dic
    for i in "ACGT":
        profile[i] = []
        for j in range(k):
            profile[i].append(c[i][j]/t)                                
    
    return profile

def Consensus(Motifs):
    cs = ""
    k = len(Motifs[0])
    c = Count(Motifs)        
    for i in range(k):
        m = 0
        s = ""
        for sym in "ACGT":
            if c[sym][i] >m:
                    m = c[sym][i]
                    s = sym
        cs +=s
    return cs

def Score(Motifs):
    score = 0
    con = Consensus(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    for i in range(k):
        for j in range(t):
            if Motifs[j][i] != con[i]:
                score +=1
    return score

#GreedyMotifSearch
def Pr(Text, Profile):
    l = len(Text)
    p = 1
    for i in range(l):
        p *= Profile[Text[i]][i]
    return p

def ProfileMostProbableKmer(text, k, profile):
    n = len(text)
    mp = -1
    ms = ""
    for i in range(n-k+1):
        T = text[i:i+k]
        prob = Pr(T, profile)
        if prob > mp:
            mp = prob
            ms = T
    return ms

def GreedyMotifSearch(Dna, k, t):
    bestmotifs = []
    for i in range(0,t):
        bestmotifs.append(Dna[i][0:k])
    n = len(Dna[0])    
    for m in range(n-k+1):
        motifs = []
        motifs.append(Dna[0][m:m+k])
        for j in range(1,t):
            P = Profile(motifs[0:j])
            motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(motifs) < Score(bestmotifs):
            bestmotifs = motifs
    return bestmotifs   