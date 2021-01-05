% clear all
% close all
angles_of_members=[40 62.62 130 62.28 40 90 0 0];

%The Matrix has to be modified
    connection=[1 1 0 0 0 0 0 0 ;
    1 0 1 1 0 0 0 1 ;
    0 1 1 0 1 0 0 0 ;
    0 0 0 1 1 1 1 0 ;
    ];
    
%A is the coefficient matrix for the equation system. There are a total of
%8 unknowns, so the dimension of A is 8 by 8.
A=zeros(8);
%External Forces exerted on the joints, given in terms of x and y
%components: Fx on Joint 1, Fy on Joint 1, Fx on Joint 2, Fy on Joint 2,
%..., Fx on Joint 4, Fy on Joint 4
b = [0 -8829 0 0 0 0 0 0]';

%Code to fill the coefficient matrix
for i=1:4
    for j=1:2
        row_num=2*(i-1)+j;
        if j==1
            A(row_num,:)=connection(i,:).*cosd(angles_of_members);
        else
            A(row_num,:)=connection(i,:).*sind(angles_of_members);
        end
    end
    angles_of_members=angles_of_members+180*connection(i,:); % 180+degrees equals -1*radians
end

F=linsolve(A,-b);
%Solve for the unknown forces
%{
fig = uifigure;
fig.Position = [500 500 490 180];

h = uihtml(fig);
h.Position = [20 20 450 130];
h.HTMLSource = fullfile('C:\','Users','13369','Downloads','Vector_Calc_Proj','displayDataFromMATLAB.html');
h.Data = 'Hello World';
%}
