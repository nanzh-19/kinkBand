clc
%先把材料旋正，然后将应变载荷转化为应力载荷画图
count=758;
e1_unrotated=zeros(count,1);
e2_unrotated=zeros(count,1);
e12_unrotated=zeros(count,1);
defgra11_nonzero=zeros(count,1);
defgra22_nonzero=zeros(count,1);
defgra12_nonzero=zeros(count,1);
for i=1:count
    jiaodu_=jiaodu(i)*pi/180;
    e1_unrotated(i)=(e1(i)+e2(i))/2+(e2(i)-e1(i))*cos(pi-2*jiaodu_)/2;
    e2_unrotated(i)=(e1(i)+e2(i))/2+(e2(i)-e1(i))*cos(2*jiaodu_)/2;
    e12_unrotated(i)=(e2(i)-e1(i))*sin(2*jiaodu_)/2;
    defgra11_nonzero(i)=e1_unrotated(i)+1;
    defgra22_nonzero(i)=e2_unrotated(i)+1;
    defgra12_nonzero(i)=e12_unrotated(i);
end
xaxis=zeros(count,1);
yaxis=zeros(count,1);
temp=min(energyvalue);
for i=1:count
   xaxis(i)=e12_unrotated(i)/e1_unrotated(i);
   yaxis(i)=e2_unrotated(i)/e1_unrotated(i);
   energyvalue(i)=energyvalue(i)/temp;
end
plot(xaxis,yaxis,'.');