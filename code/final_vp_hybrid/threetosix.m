function J=threetosix(JJ)
J=zeros(3,6);
J(1:2,1:2)=JJ(1:2,1:2);
J(1:2,3)=JJ(1:2,6);
J(3,1:2)=JJ(6,1:2);
J(3,3)=JJ(6,6);
%
J(1:2,4:5)=JJ(1:2,7:8);
J(1:2,6)=JJ(1:2,12);
J(3,4:5)=JJ(6,7:8);
J(3,6)=JJ(6,12);
end