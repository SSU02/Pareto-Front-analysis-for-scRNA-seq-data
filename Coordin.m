function[] = Coordin(X,Y,Z,X1,Coor, OutputName)
    %first
    p1 = sqrt((Coor(:,1)-X(1)).^2 + (Coor(:,2)-X(2)).^2 + (Coor(:,3)-X(3)).^2)<6;
    Coor_withindistance_r1 = Coor(p1,:);
    C1 = round(Coor_withindistance_r1,4);
    assignin('base','arc1',C1);
    writematrix(C1,['/Users/srisruthi/Documents/ARCfiles/',OutputName,'ARC1','.csv']);
    %second
    p2 = sqrt((Coor(:,1)-Y(1)).^2 + (Coor(:,2)-Y(2)).^2 + (Coor(:,3)-Y(3)).^2)<6;
    Coor_withindistance_r2 = Coor(p2,:);
    C2 = round(Coor_withindistance_r2,4);
    assignin('base','arc2',C2);
    writematrix(C2,['/Users/srisruthi/Documents/ARCfiles/',OutputName,'ARC2','.csv']);
    %third
    p3 = sqrt((Coor(:,1)-Z(1)).^2 + (Coor(:,2)-Z(2)).^2 + (Coor(:,3)-Z(3)).^2)<6;
    Coor_withindistance_r3 = Coor(p3,:);
    C3 = round(Coor_withindistance_r3,4);
    assignin('base','arc3',C3);
    writematrix(C3,['/Users/srisruthi/Documents/ARCfiles/',OutputName,'ARC3','.csv']);
    %fourth
    p4 = sqrt((Coor(:,1)-X1(1)).^2 + (Coor(:,2)-X1(2)).^2 + (Coor(:,3)-X1(3)).^2)<6;
    Coor_withindistance_r4 = Coor(p4,:);
    C4 = round(Coor_withindistance_r4,4);
    assignin('base','arc4',C4);
    writematrix(C4,['/Users/srisruthi/Documents/ARCfiles/',OutputName,'ARC4','.csv']);
    %fifth
    %p5 = sqrt((Coor(:,1)-Y1(1)).^2 + (Coor(:,2)-Y1(2)).^2 + (Coor(:,3)-Y1(3)).^2)<6;
    %Coor_withindistance_r5 = Coor(p5,:);
    %C5 = round(Coor_withindistance_r5,4);
    %assignin('base','arc5',C5);
    %writematrix(C5,['/Users/srisruthi/Documents/ARCfiles/',OutputName,'ARC5','.csv']);
    %sixth
    %p6 = sqrt((Coor(:,1)-Z1(1)).^2 + (Coor(:,2)-Z1(2)).^2 + (Coor(:,3)-Z1(3)).^2)<6;
    %Coor_withindistance_r6 = Coor(p6,:);
    %C6 = round(Coor_withindistance_r6,4);
    %assignin('base','arc6',C6);
    %writematrix(C6,['/Users/srisruthi/Documents/ARCfiles/',OutputName,'ARC6','.csv']);
    %sEVENTH
    %p7 = sqrt((Coor(:,1)-W1(1)).^2 + (Coor(:,2)-W1(2)).^2 + (Coor(:,3)-W1(3)).^2)<6;
    %Coor_withindistance_r7 = Coor(p7,:);
    %C7 = round(Coor_withindistance_r7,4);
    %assignin('base','arc7',C7);
    %writematrix(C7,['/Users/srisruthi/Documents/ARCfiles/',OutputName,'ARC7','.csv']);
end