A=[0 0 1;  0 0 2; 0 0 3; 0 1 1; 0 1 -1; 1 1 1;...
    1 1 -1;1 -1 -1; 0 2 1; 0 2 -1; 0 -2 1;0 2 2; 0 2 -2; ...
    1 2 2;1 2 -2;1 -2 -2;-1 2 2;-1 2 -2];
% 0 1 3;0 2 3;1 2 3;...
%     1 1 3; 2 2 3;1 3 3;2 3 3; 3 3 3;0 -1 3;0 1 -3;0 -2 3;0 2 -3;...
%     1 2 3;-1 2 3;1 -2 3;1 2 -3;-1 -2 3;-1 2 -3;1 -2 -3;...
%     1 1 3;-1 1 3;-1 -1 3;-1 1 -3; 2 2 3;-2 2 3;-2 -2 3;-2 2 -3;...
%     1 3 3;-1 3 3;-1 -3 3;1 -3 -3;2 3 3;-2 3 3;-2 -3 3;2 -3 -3; 3 3 3];
B=[];
for ii=1:size(A,1)
B=[B;unique(perms(A(ii,:)),'rows')];
end
unique(B,'rows');C=[702 862 1326];C=C*B';
C(C<390)=0;C(C>3500)=0;[C,ia,ic]=unique(C','rows');
B=B(ia,:);B(C==0,:)=[];C(C==0)=[];

[P,L]=findpeaks(abs(real(Full_Spec(:,end))));
for ii=1:length(L)
    [~,I(ii)]=min(abs(w(L(ii))-C));
end
D=[B(I,:) C(I)];
