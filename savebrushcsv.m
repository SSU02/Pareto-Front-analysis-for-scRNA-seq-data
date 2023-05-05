function[] = savebrushcsv(brushed1,OutputName)
C = [brushed1];
C = round(C,4);
writematrix(C,['/Users/srisruthi/Documents/',OutputName,'.csv']);
end