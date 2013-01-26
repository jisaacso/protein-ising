function ProteinIsingWater(Seq,T,tf,intvec)

% Protein ising water uses Monte Carlo simulation of an ising model
% to predict protein structure on a 2D square lattice
%
% INPUT:
%    Seq: vector of 1s and -1s corresponding to AA sequence,
%           -1 = hydrophilic
%           +1 = hydrophobic
%
%      T: Energy of the system in KbT
%
%     Tf: Number of MC time steps
%
% intvec: Vector of interaction energies between
%           . hydrophobic and water (1/0) (-5)
%           . hydrophilic and water (-1/0) (5)
%           . same amino acids (1/1 or -1/-1) (10)
%           . opposite amino acids (-1/1) (-5)
%
%   sample commandline call:
%       ProteinIsingWater([-1,1,-1,1,1,-1,-1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,-1],10,10000,[-5,5,10,-5])
%
% OUTPUT:
%       Visualization of protein fold after tf steps
%
% AUTHOR: Joseph Isaacson
%         Brown University Box 3937
%         Providence, RI 02902
    

    %Construct matrix A0, A and J0
    A0 = sparse(length(Seq)+10,length(Seq)+10);
    A0(floor(length(A0)/2),5:4+length(Seq))=Seq;


    A = A0;
    L = length(A);
    J0 = sparse(L^2,L^2);
    %initialize J0 such that ALL adjacent Amino acids have a HUGE
    %energy penality for switching (aka AA bonds are never broken)

    %find all amino acids
    x = [find(A==1);find(A==-1)];

    %for each amino acid, set J0 to be 100 if your neighbor is another aa
    for(i=x')
        %find all neighbors of n1
        nhbrs = [i-L-1:i-L+1,i-1,i+1,i+L-1:i+L+1];
        clear n2;
        idx=1;
        %for each neighbor of n1
        for(k=nhbrs)
            %if the neighbor is an amino acid
            if(A(k)==1 || A(k)==-1)
                %mark it
                n2(idx) = k;
                idx = idx+1;
            end
        end
        J0(i,n2)=100;
    end



    for(tStep=1:tf)
        %find all amino acids
        x = [find(A==1);find(A==-1)];
        %randomly pick one to be flipped
        i = x(randi(length(x),1));      %i = random index in A corresponding to a +/- 1
        [ix,iy] = ind2sub(size(A),i);


        %if this amino acid is not near an edge
        if(ix>=3 && ix<=L-2 && iy>=3 && iy<=L-2)
            %find the 8 indices in A corresponding to the 8 nhbrs of i
            nhbrsi = [(i-L-1:i-L+1)',(i-1:i+1)',(i+L-1:i+L+1)'];
            zeroEntries = find(A==0);

            %pick j - zero nhbr of i
            idx=1;
            zeroNhbrs = inf;
            for(k = 1:length(zeroEntries))
                if(any(any(zeroEntries(k) == nhbrsi)))
                    zeroNhbrs(idx)=nhbrsi(find(zeroEntries(k) == nhbrsi));
                    idx=idx+1;
                end
            end
            if(zeroNhbrs(1)==inf)
                clear zeroNhbrs;
                continue;
            end
            r = randi(length(zeroNhbrs),1); %random index in nhbrs
            j = zeroNhbrs(r);  % j = random index in A corresponding to 0 neighbor of i
            
            clear zeroNhbrs;


            %Update J
            J = J0;
            nhbrsj = [j-L-1:j-L+1,j-1,j+1,j+L-1:j+L+1];
            clear n2;
            idx=1;
            for(k=nhbrsj)
                if(A(k)==1 || A(k)==-1)
                    n2(idx) = k;
                    idx = idx+1;
                end
            end
            for(k=n2)
                %only if i,j were originally neighbors
                if(J0(i,k)==100)
                    J(j,k)=100;
                    J(k,j)=100;
                end
            end
            nhbrsi = [i-L-1:i-L+1,i-1,i+1,i+L-1:i+L+1];
            clear n2;
            idx=1;
            for(k=nhbrsi)
                J(i,k)=0;
                J(k,i)=0;
            end

            %Calc Energies
            dEs = -h(A,i,nhbrsi,intvec)+h(A,j,nhbrsj,intvec);     %change in spin energies
            dEa = h2(J0,i,nhbrsi)-h2(J,j,nhbrsj);         %energy penality for breaking AA bond
            dE = dEs+dEa;
            if(dE<=0)
                A(j) = A(i);
                A(i) = 0;
                J0=J;
                PlotMy(A);
            elseif(rand()<exp(-dE/T))
                A(j) = A(i);
                A(i) = 0;
                J0 = J;
                PlotMy(A);
            else
                %try recurstion to move all nhbrs at once
            end

        else
%             out = 'WARNING HIT EDGE OF MATRIX!'
%             PlotMy(A);
%             pause()
        end
        
    end
    PlotMy(A);
end

%Calculate local spin energies
function e = h(A,i,nhbrs,J)
    e = 0;
    for(j = nhbrs)
          if((A(i)==1&&A(j)==0) || (A(j)==1&&A(i)==0))
              e = e+2*J(1);
          elseif((A(i)==-1&&A(j)==0) || (A(j)==-1&&A(i)==0))
              e = e+2*J(2);
          elseif((A(i)==1&&A(j)==1) || (A(j)==-1&&A(i)==-1))
              e = e+2*J(3);
          elseif((A(i)==-1&&A(j)==1) || (A(j)==1&&A(i)==-1))
              e = e+2*J(4);
          end
    end
end

%Calculate local amino acid covalent bond energy
function e = h2(J,i,nhbrs)
    e = 0;
    for(j=nhbrs)
        e = e+J(i,j);
    end
end

function PlotMy(A)
    toPlot1 = A;
    toPlot2 = A;
    toPlot1(toPlot1<0) = 0;
    toPlot2(toPlot2>0) = 0;
    spy(toPlot1);
    hold on;
    spy(toPlot2,'r');
    legend('Hydrophobic','Hydrophilic');
    hold off;
    drawnow;
    pause(.1)
end