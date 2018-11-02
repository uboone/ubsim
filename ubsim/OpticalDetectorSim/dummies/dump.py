import sys
from ROOT import TChain
ch=TChain("tree")
ch.AddFile(sys.argv[1])
ch.GetEntry(0)
print ch.max_v.size()
for i in range(ch.max_v.size()):
    print ch.ch_v[i],ch.max_v[i],ch.ts_v[i]
