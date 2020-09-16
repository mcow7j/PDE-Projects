function coarse=restrict3(fine);

[r,c] = size(fine);
   m = c-1;
   m2 = m/2;
   
   

coarse = zeros(2,m2+1);

coarse(1,2:m2)=0.5*fine(1,3:2:m-1)+0.25*fine(1,2:2:m-2)+0.25*fine(1,4:2:m);
coarse(1,1)=fine(1,1)*0.75+fine(1,2)*0.25;
coarse(1,m2+1)=fine(1,m+1)*0.75+fine(1,m)*0.25;

coarse(2,2:m2)=0.5*fine(2,3:2:m-1)+0.25*fine(2,2:2:m-2)+0.25*fine(2,4:2:m);
coarse(2,1)=fine(2,1)*0.75+fine(2,2)*0.25;
coarse(2,m2+1)=fine(2,m+1)*0.75+fine(2,m)*0.25;

end