%先把材料旋正，然后将应变载荷转化为应力载荷画图

ef=324;  %纤维的杨氏模量
em=3;    %基体的杨氏模量
vf=0.3;  %纤维的泊松比
vm=0.3;  %基体的泊松比
wf=0.6;  %纤维的体积分数
wm=1-wf; %基体的体积分数
gf=0.5*ef/(1+vf);  %纤维的剪切模量 
gm=0.5*em/(1+vm);  %基体的剪切模量
s=zeros(3,3);      %初始化柔度矩阵
s(1,1)=(wf*ef+wm*em+(wf*wm*ef*em*(vf-vm)^2)/(wf*ef*(1-vm^2)+wm*em*(1-vf^2)))^-1;     %柔度矩阵分量S11
s(2,2)=wf/ef+wm/em-(2*wf*wm*(vf*em-vm*ef)^2)/((1-vf)*wm*ef*em^2+(1-vm)*wf*em*ef^2);  %柔度矩阵分量S22
s(3,3)=wf/gf+wm/gm;%柔度矩阵分量S33
s(1,2)=(wf*vf+wm*vm-vf*vm)/(wf*vm*ef+wm*vf*em-wf*ef-wm*em);  %柔度矩阵分量S12
s(2,1)=s(1,2);     %柔度矩阵分量S21
c=s^-1;            %求解刚度矩阵
count=758;
cigm11=zeros(count,1);
cigm22=zeros(count,1);
cigm12=zeros(count,1);
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
%small deformation
for i=1:count
        cigm11(i)=c(1,1)*e1_unrotated(i)+c(1,2)*e2_unrotated(i)+2*c(1,3)*e12_unrotated(i);
        cigm22(i)=c(1,2)*e1_unrotated(i)+c(2,2)*e2_unrotated(i)+2*c(2,3)*e12_unrotated(i);
        cigm12(i)=c(1,3)*e1_unrotated(i)+c(2,3)*e2_unrotated(i)+2*c(3,3)*e12_unrotated(i); 
end
%larde deformation
% for i=1:count
%    Stress=zeros(2,2);
%    Stress(1,1)=c(1,1)*0.5*(defgra11_nonzero(i)^2+defgra12_nonzero(i)^2-1)+c(1,2)*0.5*(defgra12_nonzero(i)^2+defgra22_nonzero(i)^2-1)+c(1,3)*(defgra11_nonzero(i)*defgra12_nonzero(i)+defgra22_nonzero(i)*defgra12_nonzero(i));
%    Stress(2,2)=c(1,2)*0.5*(defgra11_nonzero(i)^2+defgra12_nonzero(i)^2-1)+c(2,2)*0.5*(defgra12_nonzero(i)^2+defgra22_nonzero(i)^2-1)+c(2,3)*(defgra11_nonzero(i)*defgra12_nonzero(i)+defgra22_nonzero(i)*defgra12_nonzero(i));
%    Stress(1,2)=c(1,3)*0.5*(defgra11_nonzero(i)^2+defgra12_nonzero(i)^2-1)+c(2,3)*0.5*(defgra12_nonzero(i)^2+defgra22_nonzero(i)^2-1)+c(3,3)*(defgra11_nonzero(i)*defgra12_nonzero(i)+defgra22_nonzero(i)*defgra12_nonzero(i));
%    Stress(2,1)=Stress(1,2);
%    F=[defgra11_nonzero(i) defgra12_nonzero(i);defgra12_nonzero(i) defgra22_nonzero(i)];
%    cigmtem=F*Stress*F'/det(F);
%    cigm11(i)=cigmtem(1,1);
%    cigm22(i)=cigmtem(2,2);
%    cigm12(i)=cigmtem(1,2);
% end
plot3(cigm11,cigm22,cigm12,'.');
xlabel('\sigma_1');
ylabel('\sigma_2');
zlabel('\sigma_1_2');
%% 提取剪应力为0时的应力加载
sigma11=zeros(count,1);
sigma22=zeros(count,1);
jd=zeros(count,1);
energy=zeros(count,1);
for i=1:count
   if abs(cigm12(i))<10e-6
       if cigm11(i)~=0||cigm22(i)~=0
      sigma11(i)=cigm11(i);
      sigma22(i)=cigm22(i);
      jd(i)=jiaodu(i);
      energy(i)=energyvalue(i);
       end
   end
end
plot(sigma11,sigma22,'.')
title('没有切应力时的失稳加载');
xlabel('\sigma_1');
ylabel('\sigma_2');
%% 提取2方向应力为0时的应力加载
sigma11=zeros(1520/2,1);
sigma12=zeros(1520/2,1);
jd=zeros(1520/2,1);
energy=zeros(1520/2,1);
k=1;
for i=1:count
   if abs(cigm22(i))<10e-6
       if cigm11(i)~=0||cigm22(i)~=0
      sigma11(k)=cigm11(i);
      sigma12(k)=cigm12(i);
      jd(k)=jiaodu(i);
      energy(k)=energyvalue(i);
      k=k+1;
       end
   end
end
plot(sigma11,sigma12,'.')
title('\sigma_2<10^{-3}时的失稳加载');
xlabel('\sigma_1');
ylabel('\sigma_{12}');
%% 将三维数据通过比值压缩为二维
%plot x:sigma12/sigma11 y:sigma22/sigma11
yaxis=zeros(count,1);
xaxis=zeros(count,1);
for i=1:count
   xaxis(i)=cigm12(i)/cigm11(i);
   yaxis(i)=cigm22(i)/cigm11(i);
end
figure(2)
plot(xaxis,yaxis,'.');
xlabel('\sigma_{12}/\sigma_{11}');
ylabel('\sigma_{22}/\sigma_{11}');
% 进一步提取22和12较小时的图像
k=0;
for i=1:count 
    if xaxis(i)<=0.1 && abs(yaxis(i))<=0.1
        k=k+1;
    end 
end
xaxis_mini=zeros(k,1);
yaxis_mini=zeros(k,1);
energy_mini=zeros(k,1);
j=1;
for i=1:count 
    if xaxis(i)<=0.1 && abs(yaxis(i))<=0.1
        xaxis_mini(j)=xaxis(i);
        yaxis_mini(j)=yaxis(i);
        energy_mini(j)=energyvalue(i);
        j=j+1;
    end 
end
tmp=min(energy_mini);
energy_mini=energy_mini/tmp;
figure(3)
plot(xaxis_mini,yaxis_mini,'.');
xlabel('\sigma_{12}/\sigma_{11}');
ylabel('\sigma_{22}/\sigma_{11}');