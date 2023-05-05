function[percentpoints] = percentpoints(ans,PCAfile)
for i = 1:20
    percentpoints(i) = (ans(i)./length(PCAfile))*100;
end