clc;     %清除命令行窗口中的数据
clear;   %清楚工作区数据
ef=230;  %纤维的杨氏模量
em=3;    %基体的杨氏模量
vf=0.3;  %纤维的泊松比
vm=0.3;  %基体的泊松比
wf=0.6;  %纤维的体积分数
wm=1-wf; %基体的体积分数
gf=0.5*ef/(1+vf);  %纤维的剪切模量 
gm=0.5*em/(1+vm);  %基体的剪切模量
s=zeros(3,3);      %初始化柔度矩阵
s(1,1)=(wf*ef+wm*em+(wf*wm*ef*em*(vf-vm)^2)/(wf*ef*(1-vm^2)+wm*em*(1-vf^2)))^-1;     %柔度矩阵分量S11
s(3,3)=wf/ef+wm/em-(2*wf*wm*(vf*em-vm*ef)^2)/((1-vf)*wm*ef*em^2+(1-vm)*wf*em*ef^2);  %柔度矩阵分量S22
s(5,5)=wf/gf+wm/gm;%柔度矩阵分量S33
s(1,3)=(wf*vf+wm*vm-vf*vm)/(wf*vm*ef+wm*vf*em-wf*ef-wm*em);  %柔度矩阵分量S12
s(1,2)=(wf*vf*ef+wm*em*vm-wf*vf*ef*vm^2-wm*em*vm*vf^2)/((wf*vm*ef+wm*vf*em)^2-(wf*ef+wm*em)^2);
s(3,1)=s(1,3);     %柔度矩阵分量S21
s(2,1)=s(1,2);
s(2,2)=s(1,1);
s(2,3)=s(1,3);
s(3,2)=s(2,3);
s(4,4)=s(5,5);
s(6,6)=2*(s(1,1)-s(1,2));
c=s^-1;            %求解刚度矩阵


