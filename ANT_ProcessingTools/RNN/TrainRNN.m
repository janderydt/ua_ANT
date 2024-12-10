function TrainRNN

% we call a separate tensorflow script

pct = [960 970 980 990 995];

for ii=1:numel(pct)
    
        pyrunfile("RNN_network_optimisation_JDR.py "+string(pct(ii)));

end