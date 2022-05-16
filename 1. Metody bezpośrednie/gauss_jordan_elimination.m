% function x = gauss_jordan_elimination(A,b)
% 
%  n=size(b,1);
%  x=zeros(n,1);
%  Ab=[A b];
%  for i=1:n
%      Ab(i,:)=Ab(i,:)./Ab(i,i);
%      for j=1:n
%          if i~= j
%              Ab(j,:)=Ab(j,:)-Ab(i,:)*Ab(i,j);
%          end
%      end
%  end
% end
function x = gauss_jordan_elimination(A,b)
s = length(A);
for j = 1:(s-1)
    for i = s:-1:j+1
        m = A(i,j)/A(j,j);
        A(i,:) = A(i,:) - m*A(j,:);
        b(i) = b(i) - m*b(j);
    end
end 
x = zeros(s,1);
x(s) = b(s)/A(s,s);               
for i = s-1:-1:1                    
    sum = 0;
    for j = s:-1:i+1                
        sum = sum + A(i,j)*x(j);    
    end 
    x(i) = (b(i)- sum)/A(i,i);
end 