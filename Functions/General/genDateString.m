function dateString = genDateString(n)
%Written by David Rothlein
X=clock;
yr=num2str(X(1));

if X(2)<10
    month=['0',num2str(X(2))];
else
    month=num2str(X(2));
end;

if X(3)<10
    day=['0',num2str(X(3))];
else
    day=num2str(X(3));
end;

dateString=[yr,month,day];
