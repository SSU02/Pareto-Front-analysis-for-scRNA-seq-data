function[] = mergebrushedsri3(brushed1,brushed2,brushed3,OutputName)
C = [brushed1;brushed2;brushed3];
C = round(C,4);
writematrix(C,['/Users/srisruthi/Documents/',OutputName,'.csv']);
end