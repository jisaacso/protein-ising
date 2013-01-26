function ProteinIsing2(A0,T,tf)

%A0 should be square
A = A0;
L = length(A);
J0 = sparse(L^2,L^2);

%initialize J0 such that ALL adjacent Amino acids have a HUGE
%energy penality for switching (aka AA bonds are never broken)

%find all amino acids
x = [find(A==1);find(A==-1)];
%for each amino acid, initially set J0 to 1/2 for every neighbor 
for(i=1:length(x))
    J0(x(i),x(i)-1)=1/2;
    J0(x(i)-1,x(i))=1/2;
    
    J0(x(i),x(i)+1)=1/2;
    J0(x(i)+1,x(i))=1/2;
    
    J0(x(i),x(i)-L)=1/2;
    J0(x(i)-L,x(i))=1/2;
    
    J0(x(i),x(i)+L)=1/2;
    J0(x(i)+L,x(i))=1/2;
    
    J0(x(i),x(i)-L-1)=1/2;
    J0(x(i)-L-1,x(i))=1/2;
    
    J0(x(i),x(i)-L+1)=1/2;
    J0(x(i)-L+1,x(i))=1/2;
    
    J0(x(i),x(i)+L-1)=1/2;
    J0(x(i)+L-1,x(i))=1/2;
    
    J0(x(i),x(i)+L+1)=1/2;
    J0(x(i)+L+1,x(i))=1/2;
end
%for each amino acid, set J0 to be 5000 if your neighbor is another aa
for(i=1:length(x))
    %n1 = index in A of one AA
    n1 = x(i);
    %find all neighbors of n1
    nhbrs = [x(i)-L-1:x(i)-L+1,x(i)-1,x(i)+1,x(i)+L-1:x(i)+L+1];
    idx=1;
    %for each neighbor of n1
    for(j=1:length(nhbrs))
        %if the neighbor is an amino acid
        if(A(nhbrs(j))==1 || A(nhbrs(j))==-1)
            %mark it
            n2(idx) = nhbrs(j);
            idx = idx+1;
        end
    end
    %for each marked neighbor, change J0 from 1/2 -> 5000
    for(j=1:length(n2))
        J0(n1,n2)=5000;
    end
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
        nhbrs = [(i-L-1:i-L+1)',(i-1:i+1)',(i+L-1:i+L+1)'];
        zeroEntries = find(A==0);
        
        %randomly pick a neighboring zero, j, to exchange with
        %(only exchanges AA with 0s to preserve AA order: prevent
        %adjacent AAs from swapping)
        idx=1;
        zeroNhbrs = inf;
        for(k = 1:length(zeroEntries))
            if(any(any(zeroEntries(k) == nhbrs)))
                zeroNhbrs(idx)=nhbrs(find(zeroEntries(k) == nhbrs));
                idx=idx+1;
            end
        end
        if(zeroNhbrs(1)==inf)
            clear zeroNhbrs;
            continue;
        end
        r = randi(length(zeroNhbrs),1); %random index in nhbrs
        j = zeroNhbrs(r);  % j = random index in A corresponding to 0 neighbor of i
        
        %construct matrix of 25 neighbors about i
        nhbrs2 = [(i-2*L-2:i-2*L+2)',(i-L-2:i-L+2)',(i-2:i+2)',(i+L-2:i+L+2)',(i+2*L-2:i+2*L+2)'];
        %ai = submatrix of A: the 25 nhbrs about i
        ai = A(nhbrs2);
        %swap i,j to make af
        A2 = A;
        A2(i) = 0;
        A2(j) = A(i);
        af = A2(nhbrs2);
        %swap J(i,j)s
        J = J0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %i,j = indices in A that were switched
        
        %initially set J to be 1/2 for all neighbors of i and j
        for(idx=[i,j])
            J(idx,idx-1)=1/2;
            J(idx-1,idx)=1/2;

            J(idx,idx+1)=1/2;
            J(idx+1,idx)=1/2;

            J(idx,idx-L)=1/2;
            J(idx-L,idx)=1/2;

            J(idx,idx+L)=1/2;
            J(idx+L,idx)=1/2;

            J(idx,idx-L-1)=1/2;
            J(idx-L-1,idx)=1/2;

            J(idx,idx-L+1)=1/2;
            J(idx-L+1,idx)=1/2;

            J(idx,idx+L-1)=1/2;
            J(idx+L-1,idx)=1/2;

            J(idx,idx+L+1)=1/2;
            J(idx+L+1,idx)=1/2;
        end
        %now update J
        for(idx=[i,j])    %for each AA pair
            n1 = idx;
            nhbrs = [idx-L-1:idx-L+1,idx-1,idx+1,idx+L-1:idx+L+1];
            idx=1;
            for(k=1:length(nhbrs))
                if(A(nhbrs(k))==1 || A(nhbrs(k))==-1)
                    n2(idx) = nhbrs(k);
                    idx = idx+1;
                end
            end
            for(k=1:length(n2))
                %only if i,j were originally neighbors
                if(idx==i)
                    if(J0(i,n2)==5000)
                        J(j,n2)=5000;
                        J(n2,j)=5000;
                        J(i,n2)=1/2;
                        J(n2,i)=1/2;
                    else
                        J(j,n2)=10;
                        J(n2,j)=10;
                        J(i,n2)=1/2;
                        J(n2,i)=1/2;
                    end
                else
                    
                end
            end
        end
        
        
        
        
        
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dE = h(af,J,nhbrs2)-h(ai,J0,nhbrs2);
        if(dE<=0)
            A=A2;
            J0=J;
        elseif(rand()<exp(-dE/T))
            A=A2;
            J0=J;
        else
            dE
        end
        
        
        toPlot1 = A;
        toPlot2 = A;
        toPlot1(toPlot1<0) = 0;
        toPlot2(toPlot2>0) = 0;
        spy(toPlot1);
        legend('positive')
        hold on;
        spy(toPlot2,'r');
        hold off;
        drawnow;
        pause(.1)
        
    end
%     pause(.1)
end

end

function e = h(a,J,nhbrsInA)
    L = 5;
    nhbrs = [2:5,7:10,12:15,17:20];
    e = 0;
    for(i = nhbrs)
        e = e+J(nhbrsInA(i),nhbrsInA(i-1))*((a(i)==a(i-1))+(a(i)==-a(i-1)));
        e = e+J(nhbrsInA(i),nhbrsInA(i+L))*((a(i)==a(i+L))+(a(i)==-a(i+L)));
        e = e+J(nhbrsInA(i),nhbrsInA(i+L-1))*((a(i)==a(i+L-1))+(a(i)==-a(i+L-1)));
    end
    for(i=[1,6,11,16])
        e = e+J(nhbrsInA(i),nhbrsInA(i+L))*((a(i)==a(i+L))+(a(i)==-a(i+L)));
    end
    for(i=22:25)
        e = e+J(nhbrsInA(i),nhbrsInA(i-1))*((a(i)==a(i-1))+(a(i)==-a(i-1)));
    end
    %need energy at boundaries
    e = -e;
end