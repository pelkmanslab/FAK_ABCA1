function savepoints=MF_DrawLinesBetwnNodes_coord( X1, Y1, X2, Y2)
% Connect two pixels in a matrice with 1
%[MF] gives now normal coordinates, faster as doesn't use unique anymore
%and doesn't require a structure as input
% Command line
% ------------
% result=linept(matrice, X1, Y1, X2, Y2)
%   matrice : matrice where I'll write
%   (X1, Y1), (X2, Y2) : points to connect
%   result : matrix + the line
%
% Note
% ----
%   matrice can contents anything
%   (X1, Y1), (X2, Y2) can be out of the matrice
%
% Example
% -------
% a = linept(zeros(5, 10), 2, 2, 3, 9)
% a =
% 
%      0     0     0     0     0     0     0     0     0     0
%      0     1     1     1     1     0     0     0     0     0
%      0     0     0     0     0     1     1     1     1     0
%      0     0     0     0     0     0     0     0     0     0
%      0     0     0     0     0     0     0     0     0     0
%
% Georges Cubas 20/11/03
% georges.c@netcourrier.com
% Version 1.0
%matrice=zeros(16,16);

%X1=10;Y1=10;X2=2;Y2=3;
savepoints1=zeros(abs(1+(max(1,X2)-max(1,X1))),2);
savepoints2=zeros(abs(1+(max(1,Y2)-max(1,Y1))),2);

count1=0;count2=0;
for x=max(1, X1):sign(X2 - X1):max(1, X2)
    y = round(f_testfuntion(x, X1, Y1, X2, Y2));
    if y > 0
        count1=count1+1;%result(x, y) = 300;
        savepoints1(count1,:)=[x,y];
    end
end
savepoints1(savepoints1(:,1)==0,:)=[];
for y=max(1, Y1):sign(Y2 - Y1):max(1, Y2)
   x = round(f_testfuntion2(y, X1, Y1, X2, Y2));
   if x > 0
       %result(x, y) = 300;
       count2=count2+1;
       savepoints2(count2,:)=[x,y];
   end
end
savepoints2(savepoints2(:,1)==0,:)=[];



savepoints=[savepoints1;savepoints2];


function y=f_testfuntion(x, X1, Y1, X2, Y2)
a = (Y2 - Y1)/(X2 - X1);
b = Y1 - X1 * a;
y = a * x + b;

function x=f_testfuntion2(y, X1, Y1, X2, Y2)
if X1==X2
    x = X1;
else
	a = (Y2 - Y1)/(X2 - X1);
	b = Y1 - X1 * a;
	x = (y - b)/a;
end


