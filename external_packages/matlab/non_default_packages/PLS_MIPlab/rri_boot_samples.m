function B=permutation(n, dif_subjs)
% This script is used for permutations.

B=[];

bootstrap_num=0;
k=n;

for i = 1:n,
    A(i) = 1;
end


while A(1)<n,

  for i = 1:n, 
      array(i) = 0;
  end

  for i = 1:n, 
      array(A(i))=1;
  end


  sum=0;

  for i = 1:n,
      sum=sum+array(i);
  end

  

  if sum == dif_subjs
     for i = 1:n,
%         disp(A(i));
B=[B,A(i)];
     end

%     disp(''\n'');
     bootstrap_num = bootstrap_num+1;
  end 



  if A(k) == n
    i=k;

    while A(i) == n,
        i = i-1;
    end

    A(i)=A(i)+1;

    for j = (i+1):n,
        A(j)=A(i);
    end

  else
    A(k)=A(k)+1;
  end

end


if dif_subjs == 1
    for i= 1:n, 
%      disp(n);
B=[B,n];
    end
    bootstrap_num=bootstrap_num+1;
end

%disp(bootstrap_num);
B=[reshape(B,[n,length(B)/n])]';

