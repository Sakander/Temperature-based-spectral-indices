close all
clc
format short
n=32;            % Order of the graph

B=  [[0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]];

A= reshape(B,[n,n]);
if A==A'
      disp('Matrix is Symmetric')
else
      disp('Not Symmetric')
end

l=size(A,1);
e=[];
s=[];
t=[];
k=[];
for i=1:l
    for j=i+1:l
        if A(i,j)==1
            e=[e;[i,j]];
        end
    end
end

edges=size(e,1);
j=ones(edges,1);

k=sum(A);
d=k';
t=d./(n-d);

% The HT1 matrix
A1=[];
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=(t(i)+t(j))^2;
            A1=A;
        end
    end
end
q1=eig(A1);
rho1=max(q1);
E1=sum(abs(eig(A1)));
EE1=sum(exp(eig(A1)));

% The HT2 matrix
A2=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=(t(i)*t(j))^2;
            A2=A;
        end
    end
end
q2=eig(A2);
rho2=max(q2);
E2=sum(abs(eig(A2)));
EE2=sum(exp(eig(A2)));

% The ST matrix
A3=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=1./sqrt(t(i)+t(j));
            A3=A;
        end
    end
end
q3=eig(A3);
rho3=max(q3);
E3=sum(abs(eig(A3)));
EE3=sum(exp(eig(A3)));

% The PT matrix
A4=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=1./sqrt(t(i)*t(j));
            A4=A;
        end
    end
end
q4=eig(A4);
rho4=max(q4);
E4=sum(abs(eig(A4)));
EE4=sum(exp(eig(A4)));

% The RST matrix
A5=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=sqrt(t(i)+t(j));
            A5=A;
        end
    end
end
q5=eig(A5);
rho5=max(q5);
E5=sum(abs(eig(A5)));
EE5=sum(exp(eig(A5)));

% The RPT matrix
A6=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=sqrt(t(i)*t(j));
            A6=A;
        end
    end
end
q6=eig(A6);
rho6=max(q6);
E6=sum(abs(eig(A6)));
EE6=sum(exp(eig(A6)));


% The FT matrix
A7=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=((t(i))^2+(t(j))^2);
            A7=A;
        end
    end
end
q7=eig(A7);
rho7=max(q7);
E7=sum(abs(eig(A7)));
EE7=sum(exp(eig(A7)));

% The AGT matrix
A8=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=((t(i)+t(j))/2*sqrt(t(i)*t(j)));
            A8=A;
        end
    end
end
q8=eig(A8);
rho8=max(q8);
E8=sum(abs(eig(A8)));
EE8=sum(exp(eig(A8)));

% The TABC matrix
A9=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=sqrt(abs((t(i)+t(j)-2)/(t(i)*t(j))));
            A9=A;
        end
    end
end
q9=eig(A9);
rho9=max(q9);
E9=sum(abs(eig(A9)));
EE9=sum(exp(eig(A9)));

% The TH matrix
A10=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=2/(t(i)+t(j));
            A10=A;
        end
    end
end
q10=eig(A10);
rho10=max(q10);
E10=sum(abs(eig(A10)));
EE10=sum(exp(eig(A10)));

% The TSO matrix
A11=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=sqrt((t(i))^2+(t(j))^2);
            A11=A;
        end
    end
end
q11=eig(A11);
rho11=max(q11);
E11=sum(abs(eig(A11)));
EE11=sum(exp(eig(A11)));

% The mTSO matrix
A12=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=1./sqrt((t(i))^2+(t(j))^2);
            A12=A;
        end
    end
end
q12=eig(A12);
rho12=max(q12);
E12=sum(abs(eig(A12)));
EE12=sum(exp(eig(A12)));

% The RRPT matrix
A13=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=sqrt((t(i)-1)*(t(j)-1));
            A13=A;
        end
    end
end
q13=eig(A13);
rho13=max(q13);
E13=sum(abs(eig(A13)));
EE13=sum(exp(eig(A13)));

% The GAT matrix
A14=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=(2*sqrt(t(i)*t(j)))/((t(i)+t(j)));
            A14=A;
        end
    end
end
q14=eig(A14);
rho14=max(q14);
E14=sum(abs(eig(A14)));
EE14=sum(exp(eig(A14)));




% disp('Temp-based spectral descriptors:')
% fprintf('The HT1-SR is %4.4f\n',rho1');
% fprintf('The HT1-E is %4.4f\n',E1');
% fprintf('The HT1-EE is %4.4f\n',EE1');
% fprintf('The HT2-SR is %4.4f\n',rho2');
% fprintf('The HT2-E is %4.4f\n',E2');
% fprintf('The HT2-EE is %4.4f\n',EE2');
% fprintf('The ST-SR is %4.4f\n',rho3');
% fprintf('The ST-E is %4.4f\n',E3');
% fprintf('The ST-EE is %4.4f\n',EE3');
% fprintf('The PT-SR is %4.4f\n',rho4');
% fprintf('The PT-E is %4.4f\n',E4');
% fprintf('The PT-EE is %4.4f\n',EE4');
% fprintf('The RST-SR is %4.4f\n',rho5');
% fprintf('The RST-E is %4.4f\n',E5');
% fprintf('The RST-EE is %4.4f\n',EE5');
% fprintf('The RPT-SR is %4.4f\n',rho6');
% fprintf('The RPT-E is %4.4f\n',E6');
% fprintf('The RPT-EE is %4.4f\n',EE6');
% fprintf('The FT-SR is %4.4f\n',rho7');
% fprintf('The FT-E is %4.4f\n',E7');
% fprintf('The FT-EE is %4.4f\n',EE7');
% fprintf('The AGT-SR is %4.4f\n',rho8');
% fprintf('The AGT-E is %4.4f\n',E8');
% fprintf('The AGT-EE is %4.4f\n',EE8');
% fprintf('The TABC-SR is %4.4f\n',rho9');
% fprintf('The TABC-E is %4.4f\n',E9');
% fprintf('The TABC-EE is %4.4f\n',EE9');
% fprintf('The TH-SR is %4.4f\n',rho10');
% fprintf('The TH-E is %4.4f\n',E10');
% fprintf('The TH-EE is %4.4f\n',EE10');
% fprintf('The TSO-SR is %4.4f\n',rho11');
% fprintf('The TSO-E is %4.4f\n',E11');
% fprintf('The TSO-EE is %4.4f\n',EE11');
% fprintf('The mTSO-SR is %4.4f\n',rho12');
% fprintf('The mTSO-E is %4.4f\n',E12');
% fprintf('The mTSO-EE is %4.4f\n',EE12');
% fprintf('The RRPT-SR is %4.4f\n',rho13');
% fprintf('The RRPT-E is %4.4f\n',E13');
% fprintf('The RRPT-EE is %4.4f\n',EE13');
% fprintf('The GAT-SR is %4.4f\n',rho14');
% fprintf('The GAT-E is %4.4f\n',E14');
% fprintf('The GAT-EE is %4.4f\n',EE14');

disp('Temp-based spectral descriptors:')
fprintf('%4.4f\n',rho1');
fprintf('%4.4f\n',E1');
fprintf('%4.4f\n',EE1');
fprintf('%4.4f\n',rho2');
fprintf('%4.4f\n',E2');
fprintf('%4.4f\n',EE2');
fprintf('%4.4f\n',rho3');
fprintf('%4.4f\n',E3');
fprintf('%4.4f\n',EE3');
fprintf('%4.4f\n',rho4');
fprintf('%4.4f\n',E4');
fprintf('%4.4f\n',EE4');
fprintf('%4.4f\n',rho5');
fprintf('%4.4f\n',E5');
fprintf('%4.4f\n',EE5');
fprintf('%4.4f\n',rho6');
fprintf('%4.4f\n',E6');
fprintf('%4.4f\n',EE6');
fprintf('%4.4f\n',rho7');
fprintf('%4.4f\n',E7');
fprintf('%4.4f\n',EE7');
fprintf('%4.4f\n',rho8');
fprintf('%4.4f\n',E8');
fprintf('%4.4f\n',EE8');
fprintf('%4.4f\n',rho9');
fprintf('%4.4f\n',E9');
fprintf('%4.4f\n',EE9');
fprintf('%4.4f\n',rho10');
fprintf('%4.4f\n',E10');
fprintf('%4.4f\n',EE10');
fprintf('%4.4f\n',rho11');
fprintf('%4.4f\n',E11');
fprintf('%4.4f\n',EE11');
fprintf('%4.4f\n',rho12');
fprintf('%4.4f\n',E12');
fprintf('%4.4f\n',EE12');
fprintf('%4.4f\n',rho13');
fprintf('%4.4f\n',E13');
fprintf('%4.4f\n',EE13');
fprintf('%4.4f\n',rho14');
fprintf('%4.4f\n',E14');
fprintf('%4.4f\n',EE14');