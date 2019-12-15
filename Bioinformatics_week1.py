class PatternMatching:
	def __init__(self, Text, Pattern, k, Genome):
		self.Text = Text
		self.Pattern = Pattern
		self.k = k
		self.Genome = Genome



def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
        for i in range(n-k+1):
            if Text[i:i+k] == Pattern:
                freq[Pattern] +=1       
    return freq

def FrequentWords(Text, k):
    word = []
    freq = FrequencyMap(Text,k)
    m = max(freq.values())
    for i in freq:
         if freq[i] == m:
                word.append(i)
    return word

def ReverseComplement(Pattern):
    def Reverse(Pattern):  
        P = Pattern[::-1]
        return P

    def Complement(Pattern):        
        rev = ""
        for i in Reverse(Pattern):
            if i == "A":
                rev += "T"
            elif i == "T":
                rev += "A"
            elif i == "C":
                rev += "G"
            elif i == "G":
                rev += "C"   
            else:
                return "not exist"
        return rev
    return Complement(Pattern)

def PatternMatching(Pattern, Genome):
    positions = []
    n = len (Pattern)
    k = len (Genome)
    for i in range(k-n+1):
        if Genome[i:i+n] == Pattern:
            positions.append(i)        
    return positions




