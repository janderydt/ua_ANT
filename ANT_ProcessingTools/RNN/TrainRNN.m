function TrainRNN

% we call a separate tensorflow script

pct = [95:99];

for ii=1:numel(pct)
    
        pyrunfile("RNN_network_optimisation_JDR.py "+string(pct(ii)));

end