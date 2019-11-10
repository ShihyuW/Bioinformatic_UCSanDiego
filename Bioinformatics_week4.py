from shihyu.Bioinformatics_week3 import *

def CountWithPseudocounts(Motifs):
    Count = {} 
    k = len(Motifs[0])
    for s in "ACGT":
        Count[s] = []
        for j in range(k):
            Count[s].append(1)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            indi = Motifs[i][j]
            Count[indi][j]+=1
            
    return Count

def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {} # output variable
    c = CountWithPseudocounts(Motifs)
    for i in "ACGT":
        profile[i] = []
        for j in range(k):
            profile[i].append(c[i][j]/(t+4))
    return profile

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    bestmotifs = [] # output variable
    for i in range(0,t):
        bestmotifs.append(Dna[i][0:k])
    n = len(Dna[0])    
    for m in range(n-k+1):
        motifs = []
        motifs.append(Dna[0][m:m+k])
        for j in range(1,t):
            P = ProfileWithPseudocounts(motifs[0:j])
            motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(motifs) < Score(bestmotifs):
            bestmotifs = motifs
    return bestmotifs

def Motifs(Profile, Dna):
    pmpk = []
    k = len(Profile['A'])
    for i in range(len(Dna)):
        text = Dna[i]
        p = ProfileMostProbableKmer(text, k, Profile)
        pmpk.append(p)
    return pmpk


#Random Motifs Search
import random
def RandomMotifs(Dna, k, t):
    RandomS = []
    for i in range(t):
        a = random.randint(1,len(Dna[0])-k)
        RandomS.append(Dna[i][a:a+k+1])
    return RandomS  

def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)#kmer array (list
    Bestmotifs = M
        
    while True:
        Profile = ProfileWithPseudocounts(M)#profile of kmer array(Dictionary
        M = Motifs(Profile, Dna)#每一列的most probable kmer(list,array
        if Score(M) < Score(Bestmotifs):
            Bestmotifs = M
        else:
            return Bestmotifs

#Gibbs Sampling
def Pr(Text, Profile):
    l = len(Text)
    p = 1
    for i in range(l):
        p *= Profile[Text[i]][i]
    return p

def Normalize(Probabilities):
    Normalized = {}
    for i in Probabilities.keys():
        Normalized[i] = Probabilities[i]/sum(Probabilities.values())
    return Normalized

def WeightedDie(Probabilities):
    kmer = ""
    r = random.uniform(0,1)
    N = Normalize(Probabilities)
    p = 0
    for i in N:
        p += N[i]
        if r <= p:
            kmer = i
            break
    return kmer

def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    pr = {}
    for i in range(0,n-k+1):
        pr[Text[i:i+k]] = Pr(Text[i:i+k],profile) 
    pr = Normalize(pr)
    return WeightedDie(pr)

def GibbsSampler(Dna, k, t, N):
    BestMotifs = []
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    for j in range(N):
        i = random.randint(0,t-1)
        del M[i]
        profile = ProfileWithPseudocounts(M)
        M.insert(i, ProfileGeneratedString(Dna[i], profile, k))
        if Score(M) < Score(BestMotifs):
            BestMotifs = M   
    return BestMotifs