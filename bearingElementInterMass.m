%% bearingElementInterMass
% generate the mass, stiffness, damping matrix of a intermediate bearing 
% element with mass
%% Syntax
% [Me, Ke, Ce] = bearingElementInterMass(AMBearing)
%% Description
% AMBearing is a struct saving the physical parameters of a bearing element
% with fields: dofOfEachNodes, stiffness, damping, mass, dofOnShaftNode
%
% Me, Ke, Ce are mass, stiffness, damping cell of a intermediate bearing 
% element. 
%
% cell2mat(Me) is n*n matrix;
% 
% Ke, Ce is 1*7 cell saving the non-zero matrix of intermediate bearing:
%
% 1->Inner shaft; 2->Inner shaft with m1 (1,2); 3->Inner shaft with m1
% (2,1); 4->Outer shaft; 5->Outer shaft with mn (2,n); 6->Outer shaft with
% mn (n,2); 7->m1 m2 ... mn
%% Symbols
% m: stiffness of bearing
% 
% k: stiffness of bearing
%
% c: damping of bearing


function [Me, Ke, Ce] = bearingElementInterMass(AMBearing)

%% initial
% constants
m = AMBearing.mass;
k = AMBearing.stiffness;
c = AMBearing.damping;
dofShaft = AMBearing.dofOnShaftNode;
dofBearing = AMBearing.dofOfEachNodes;
m = m(m~=0); % get the non-zero mass
massNum = length(m);
if length(AMBearing.mass)~=massNum
    k = k(1:massNum+1);
    c = c(1:massNum+1); 
    dofBearing = dofBearing(1:massNum);% get first n+1 sitffness and damping, n is the number of the non-zero mass of bearings
end
dofElement = [dofShaft, dofBearing];
dofNum = sum(dofElement);
dofShaftNum = sum(dofShaft);
dofBearingNum = sum(dofBearing);

% initial
K = zeros(dofNum, dofNum);
C = zeros(dofNum, dofNum);


%% add part of the inner shaft (stiffness and damping)
% stiffness matrix
Kin = [ k(1), 0;...
        0, k(1) ];
K = addElementIn(K, Kin, [1,1]);
K = addElementIn(K, -Kin, [1,dofShaftNum+1]);

% damping matrix
Cin = [ c(1), 0;...
        0, c(1) ];
C = addElementIn(C, Cin, [1,1]);
C = addElementIn(C, -Cin, [1,dofShaftNum+1]);


%% add part of the outer shaft (stiffness and damping)
% stiffness matrix
Kin = [k(end), 0;...
       0,      k(end)];
K = addElementIn(K, Kin, [dofShaft(1)+1, dofShaft(1)+1]);
K = addElementIn(K, -Kin, [dofShaft(1)+1, sum(dofElement(1:end-1))+1]);

% damping matrix
Cin = [c(end), 0;...
       0,      c(end)];
C = addElementIn(C, Cin, [dofShaft(1)+1, dofShaft(1)+1]);
C = addElementIn(C, -Cin, [dofShaft(1)+1, sum(dofElement(1:end-1))+1]);

%% add part of the mass bearing (stiffness and damping)

for im = 1:1:massNum
    isLast = im==massNum; % boolean, last element
    isFirst = im==1;
    if isLast&&isFirst
        [K1, C1] = mFirstLast(k(1),k(2),c(1),c(2),dofShaft(1),dofShaft(2),dofBearing);
        K = addElementIn(K, K1, [1,1]);
        C = addElementIn(C, C1, [1,1]);
    elseif isFirst&&~isLast
        [Kj, Cj] = mj(k(1),k(2),c(1),c(2),dofShaft(1),dofBearing(1),dofBearing(2));
        Kje = mat2cell(Kj, [dofShaft(1), dofBearing(1)+dofBearing(2)], [dofShaft(1), dofBearing(1)+dofBearing(2)]);
        Cje = mat2cell(Cj, [dofShaft(1), dofBearing(1)+dofBearing(2)], [dofShaft(1), dofBearing(1)+dofBearing(2)]);
        K = addElementIn(K, Kje{2,1}, [dofShaftNum+1,1]); % Kje{1,1}, kje{2,1} are 0 matrix
        K = addElementIn(K, Kje{2,2}, [dofShaftNum+1,dofShaftNum+1]);
        C = addElementIn(C, Cje{2,1}, [dofShaftNum+1,1]);
        C = addElementIn(C, Cje{2,2}, [dofShaftNum+1,dofShaftNum+1]);
    elseif isLast&&~isFirst
        [Kn, Cn] = mn(k(end-1),k(end),c(end-1),c(end),dofShaft(2),dofBearing(end-1),dofBearing(end));
        Kne = mat2cell(Kn, [dofShaft(2), dofBearing(end-1)+dofBearing(end)], [dofShaft(2), dofBearing(end-1)+dofBearing(end)]);
        Cne = mat2cell(Cn, [dofShaft(2), dofBearing(end-1)+dofBearing(end)], [dofShaft(2), dofBearing(end-1)+dofBearing(end)]);
        dofHere = sum(dofElement(1:end-2))+1;
        K = addElementIn(K, Kne{2,1}, [dofHere,dofShaft(1)+1]); % Kne{1,1}, kne{2,1} are 0 matrix
        K = addElementIn(K, Kne{2,2}, [dofHere,dofHere]);
        C = addElementIn(C, Cne{2,1}, [dofHere,dofShaft(1)+1]); % Kne{1,1}, kne{2,1} are 0 matrix
        C = addElementIn(C, Cne{2,2}, [dofHere,dofHere]);
    else
        [Kj, Cj] = mj(k(im),k(im+1),c(im),c(im+1),dofBearing(im-1),dofBearing(im),dofBearing(im+1));
        dofHere = sum(dofElement(1:im))+1;
        K = addElementIn(K, Kj, [dofHere,dofHere]);
        C = addElementIn(C, Cj, [dofHere,dofHere]);
    end % end if
end % end for im


%% divide into 7 sub-matrix (stiffness and damping)

% there are 7 non-zeros sub-matrix in this intermediate bearing element
% 1->Inner shaft; 2->Inner shaft with m1 (1,2); 3->Inner shaft with m1
% (2,1); 4->Outer shaft; 5->Outer shaft with mn (2,n); 6->Outer shaft with
% mn (n,2); 7->m1 m2 ... mn
switch massNum
    case 1
        KeTemp = mat2cell(K, [dofShaft, dofBearing], [dofShaft, dofBearing]);
        CeTemp = mat2cell(C, [dofShaft, dofBearing], [dofShaft, dofBearing]);
        Ke = {KeTemp{1,1}, KeTemp{1,3}, KeTemp{3,1}, KeTemp{2,2}, KeTemp{2,3}, KeTemp{3,2}, KeTemp{3,3}};
        Ce = {CeTemp{1,1}, CeTemp{1,3}, CeTemp{3,1}, CeTemp{2,2}, CeTemp{2,3}, CeTemp{3,2}, CeTemp{3,3}};
    case 2
        KeTemp = mat2cell(K, [dofShaft, dofBearing], [dofShaft, dofBearing]);
        KeTemp2 = mat2cell(K, [dofShaftNum,dofBearingNum], [dofShaftNum,dofBearingNum]);
        CeTemp = mat2cell(C, [dofShaft, dofBearing], [dofShaft, dofBearing]);
        CeTemp2 = mat2cell(C, [dofShaftNum,dofBearingNum], [dofShaftNum,dofBearingNum]);
        Ke = {KeTemp{1,1}, KeTemp{1,3}, KeTemp{3,1}, KeTemp{2,2}, KeTemp{2,4}, KeTemp{4,2}, KeTemp2{2,2}};
        Ce = {CeTemp{1,1}, CeTemp{1,3}, CeTemp{3,1}, CeTemp{2,2}, CeTemp{2,4}, CeTemp{4,2}, CeTemp2{2,2}};
    otherwise % dofBearingNum >= 3
        divideIndex = [dofShaft, dofBearing(1), sum(dofBearing(2:end-1)), dofBearing(end)];
        divideIndex2 = [dofShaftNum,dofBearingNum];
        KeTemp = mat2cell(K, divideIndex, divideIndex);
        KeTemp2 = mat2cell(K, divideIndex2, divideIndex2);
        CeTemp = mat2cell(C, divideIndex, divideIndex);
        CeTemp2 = mat2cell(C, divideIndex2, divideIndex2);
        Ke = {KeTemp{1,1}, KeTemp{1,3}, KeTemp{3,1}, KeTemp{2,2}, KeTemp{2,5}, KeTemp{5,2}, KeTemp2{2,2}};
        Ce = {CeTemp{1,1}, CeTemp{1,3}, CeTemp{3,1}, CeTemp{2,2}, CeTemp{2,5}, CeTemp{5,2}, CeTemp2{2,2}};
end


%% mass

Min = zeros(sum(dofBearingNum));
for im = 1:1:massNum
    Mi = [ m(im), 0;...
            0,     m(im)];
    dofHere = sum(dofBearing(1:im-1))+1;
    Min = addElementIn(Min, Mi, [dofHere, dofHere]);
end

M11 = zeros(dofShaftNum);      M12 = zeros(dofShaftNum, dofBearingNum);
M21 = M12';                    M22 = zeros(dofBearingNum);

M22 = addElementIn(M22, Min, [1,1]);

% output
Me = {M11, M12;...
      M21, M22 };
  
  
%% sub-function
    
    % sub-function 1
    function [Kn, Cn] = mn(kn,kn1,cn,cn1,dof1,dof2,dof3)
        dofNum1 = dof1+dof2+dof3;
        % initial
        Kn = zeros(dofNum1);
        Cn = zeros(dofNum1);
        % construct 
        Kn1 = [-kn, 0;...
               0,   -kn];
        Kn2 = [kn+kn1, 0;...
               0,      kn+kn1];
        Kn3 = [-kn1, 0;...
               0,   -kn1];
        Cn1 = [-cn, 0;...
               0,   -cn];
        Cn2 = [cn+cn1, 0;...
               0,      cn+cn1];
        Cn3 = [-cn1, 0;...
               0,   -cn1];
        % assembly
        Kn = addElementIn(Kn, Kn3, [dof1+dof2+1, 1]);
        Kn = addElementIn(Kn, Kn1, [dof1+dof2+1, dof1+1]);
        Kn = addElementIn(Kn, Kn2, [dof1+dof2+1, dof1+dof2+1]);
        Cn = addElementIn(Cn, Cn3, [dof1+dof2+1, 1]);
        Cn = addElementIn(Cn, Cn1, [dof1+dof2+1, dof1+1]);
        Cn = addElementIn(Cn, Cn2, [dof1+dof2+1, dof1+dof2+1]);
    end


    % sub-function 2
    function [Kj, Cj] = mj(kj,kj1,cj,cj1,dof1,dof2,dof3)
        dofNum1 = dof1 + dof2 + dof3;
        % initial
        Kj = zeros(dofNum1);
        Cj = zeros(dofNum1);
        % construct
        Kj1 = [-kj, 0;...
               0,   -kj];
        Kj2 = [kj+kj1, 0;...
               0,      kj+kj1];
        Kj3 = [-kj1, 0;...
               0,   -kj1];  
        Cj1 = [-cj, 0;...
               0,   -cj];
        Cj2 = [cj+cj1, 0;...
               0,      cj+cj1];
        Cj3 = [-cj1, 0;...
               0,   -cj1];
        % assembly
        Kj = addElementIn(Kj, Kj1, [dof1+1, 1]);
        Kj = addElementIn(Kj, Kj2, [dof1+1, dof1+1]);
        Kj = addElementIn(Kj, Kj3, [dof1+1, dof1+dof2+1]);
        Cj = addElementIn(Cj, Cj1, [dof1+1, 1]);
        Cj = addElementIn(Cj, Cj2, [dof1+1, dof1+1]);
        Cj = addElementIn(Cj, Cj3, [dof1+1, dof1+dof2+1]);
    end
    
    % sub-function 3
    function [K1, C1] = mFirstLast(k1, k2, c1, c2, dof1, dof2, dof3)
        dofNum1 = dof1 + dof2 + dof3;
        % initial
        K1 = zeros(dofNum1);
        C1 = zeros(dofNum1);
        % construct
        K11 = [-k1, 0;...
               0,   -k1];
        K12 = [-k2, 0;...
               0,   -k2];
        K13 = [k1+k2, 0;...
               0,     k1+k2];
        C11 = [-c1, 0;...
               0,   -c1];
        C12 = [-c2, 0;...
               0,   -c2];
        C13 = [c1+c2, 0;...
               0,     c1+c2];
        % assembly
        K1 = addElementIn(K1, K11, [dof1+dof2+1, 1]);
        K1 = addElementIn(K1, K12, [dof1+dof2+1, dof1+1]);
        K1 = addElementIn(K1, K13, [dof1+dof2+1, dof1+dof2+1]);
        C1 = addElementIn(C1, C11, [dof1+dof2+1, 1]);
        C1 = addElementIn(C1, C12, [dof1+dof2+1, dof1+1]);
        C1 = addElementIn(C1, C13, [dof1+dof2+1, dof1+dof2+1]);
    end
end  