function[] = mergebrushedsri2(brushed1,brushed2,OutputName)
C = [brushed1;brushed2];
C = round(C,4);
writematrix(C,['/Users/srisruthi/Documents/',OutputName,'.csv']);
end