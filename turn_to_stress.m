%% 
ef=324;  %��ά������ģ��
em=3;    %���������ģ��
vf=0.3;  %��ά�Ĳ��ɱ�
vm=0.3;  %����Ĳ��ɱ�
wf=0.6;  %��ά���������
wm=1-wf; %������������
gf=0.5*ef/(1+vf);  %��ά�ļ���ģ�� 
gm=0.5*em/(1+vm);  %����ļ���ģ��
s=zeros(3,3);      %��ʼ����Ⱦ���
s(1,1)=(wf*ef+wm*em+(wf*wm*ef*em*(vf-vm)^2)/(wf*ef*(1-vm^2)+wm*em*(1-vf^2)))^-1;     %��Ⱦ������S11
s(2,2)=wf/ef+wm/em-(2*wf*wm*(vf*em-vm*ef)^2)/((1-vf)*wm*ef*em^2+(1-vm)*wf*em*ef^2);  %��Ⱦ������S22
s(3,3)=wf/gf+wm/gm;%��Ⱦ������S33
s(1,2)=(wf*vf+wm*vm-vf*vm)/(wf*vm*ef+wm*vf*em-wf*ef-wm*em);  %��Ⱦ������S12
s(2,1)=s(1,2);     %��Ⱦ������S21
c=s^-1;            %���նȾ���
count=758;
cigm11=zeros(count,1);
cigm22=zeros(count,1);
cigm12=zeros(count,1);
cigm1=zeros(count,1);
cigm2=zeros(count,1);
arfa=zeros(count,1);
defgra11_nonzero=e1+1;
defgra22_nonzero=e2+1;
for i=1:count
        jiaodu_=jiaodu(i)*pi/180;
        c_rotated(1,1)=c(1,1)*cos(jiaodu_)^4+2*(c(1,2)+2*c(3,3))*sin(jiaodu_)^2*cos(jiaodu_)^2+c(2,2)*sin(jiaodu_)^4;
        c_rotated(1,2)=(c(1,1)+c(2,2)-4*c(3,3))*sin(jiaodu_)^2*cos(jiaodu_)^2+c(1,2)*(sin(jiaodu_)^4+cos(jiaodu_)^4);
        c_rotated(2,2)=c(1,1)*sin(jiaodu_)^4+2*(c(1,2)+2*c(3,3))*sin(jiaodu_)^2*cos(jiaodu_)^2+c(2,2)*cos(jiaodu_)^4;
        c_rotated(1,3)=(c(1,1)-c(1,2)-2*c(3,3))*sin(jiaodu_)*cos(jiaodu_)^3+(c(1,2)-c(2,2)+2*c(3,3))*sin(jiaodu_)^3*cos(jiaodu_);
        c_rotated(2,3)=(c(1,1)-c(1,2)-2*c(3,3))*sin(jiaodu_)^3*cos(jiaodu_)+(c(1,2)-c(2,2)+2*c(3,3))*sin(jiaodu_)*cos(jiaodu_)^3;
        c_rotated(3,3)=(c(1,1)+c(2,2)-2*c(1,2)-2*c(3,3))*sin(jiaodu_)^2*cos(jiaodu_)^2+c(3,3)*(sin(jiaodu_)^4+cos(jiaodu_)^4);
        c_rotated(2,1)=c_rotated(1,2);
        c_rotated(3,1)=c_rotated(1,3); 
        c_rotated(3,2)=c_rotated(2,3);
%         cigm11(i)=defgra11_nonzero(i)*(c_rotated(1,1)*(defgra11_nonzero(i)-1)+c_rotated(1,2)*(defgra22_nonzero(i)-1))/defgra22_nonzero(i);
%         cigm22(i)=defgra22_nonzero(i)*(c_rotated(1,2)*(defgra11_nonzero(i)-1)+c_rotated(2,2)*(defgra22_nonzero(i)-1))/defgra11_nonzero(i);
%         cigm12(i)=c_rotated(1,3)*(defgra11_nonzero(i)-1)+c_rotated(2,3)*(defgra22_nonzero(i)-1);
        cigm11(i)=defgra11_nonzero(i)*(0.5*(c_rotated(1,1)*(defgra11_nonzero(i)^2-1)+0.5*c_rotated(1,2)*(defgra22_nonzero(i)^2-1)))/defgra22_nonzero(i);
        cigm22(i)=defgra22_nonzero(i)*(0.5*c_rotated(1,2)*(defgra11_nonzero(i)^2-1)+0.5*c_rotated(2,2)*(defgra22_nonzero(i)^2-1))/defgra11_nonzero(i);
        cigm12(i)=0.5*c_rotated(1,3)*(defgra11_nonzero(i)^2-1)+0.5*c_rotated(2,3)*(defgra22_nonzero(i)^2-1);       
%         if cigm11(i)>cigm22(i)
        cigm1(i)=(cigm11(i)+cigm22(i))/2+sqrt(((cigm11(i)-cigm22(i))/2)^2+cigm12(i)^2);
        cigm2(i)=(cigm11(i)+cigm22(i))/2-sqrt(((cigm11(i)-cigm22(i))/2)^2+cigm12(i)^2);
        arfa(i)=abs(0.5*atan(-2*cigm12(i)/(cigm22(i)-cigm11(i))))*180/pi;
%         else 
%         cigm2(i)=(cigm11(i)+cigm22(i))/2+sqrt(((cigm11(i)-cigm22(i))/2)^2+cigm12(i)^2);
%         cigm1(i)=(cigm11(i)+cigm22(i))/2-sqrt(((cigm11(i)-cigm22(i))/2)^2+cigm12(i)^2);
%         arfa(i)=abs(0.5*atan(-2*cigm12(i)/(cigm11(i)-cigm22(i))))*180/pi;    
%         end
% 
end
%%
figure(1)
plot3(cigm11,cigm22,cigm12,'.');
xlabel('\sigma_1');ylabel('\sigma_2');
zlabel('\sigma_{12}');

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
% ��һ����ȡ22��12��Сʱ��ͼ��
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
%% ��ȡ��Ӧ��Ϊ0ʱ��Ӧ������
sigma11=zeros(1520/2,1);
sigma22=zeros(1520/2,1);
jd=zeros(1520/2,1);
energy=zeros(1520/2,1);
for i=1:count
   if abs(cigm12(i))<10e-6
       if cigm11(i)~=0||cigm22(i)~=0
      sigma11(i)=cigm11(i);
      sigma22(i)=cigm22(i);
      jd(i)=arfa(i);
      energy(i)=energyvalue(i);
       end
   end
end
plot(sigma11,sigma22,'.')
title('û����Ӧ��ʱ��ʧ�ȼ���');
xlabel('\sigma_1');
ylabel('\sigma_2');
%% ��ȡ2����Ӧ��Ϊ0ʱ��Ӧ������
sigma11=zeros(1520/2,1);
sigma12=zeros(1520/2,1);
jd=zeros(1520/2,1);
energy=zeros(1520/2,1);
k=1;
for i=1:count
   if abs(cigm22(i))<10e-3
       if cigm11(i)~=0||cigm22(i)~=0
      sigma11(k)=cigm11(i);
      sigma12(k)=cigm12(i);
      jd(k)=arfa(i);
      energy(k)=energyvalue(i);
      k=k+1;
       end
   end
end
plot(sigma11,sigma12,'.')
title('\sigma_2<10^{-3}ʱ��ʧ�ȼ���');
xlabel('\sigma_1');
ylabel('\sigma_{12}');
