import sys
from ROOT import TChain
ch=TChain("tree")
ch.AddFile(sys.argv[1])
ch.GetEntry(0)

max_v = ch.max_v
tp_v  = ch.tpeak_v
ts_v  = ch.tstart_v
ch_v  = ch.ch_v
print ch.max_v.size()
for i in range(ch.max_v.size()):
    print 'waveform',i,'Ch.',ch_v[i],'Tstart',ts_v[i],'[us] Max ADC',max_v[i],'@ index',tp_v[i]
