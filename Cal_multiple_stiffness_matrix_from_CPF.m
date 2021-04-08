% ��CPF�������õ��ն���
clear all;
clc;
global c
calC();
f11=0.9852;
f22=sqrt(-c(2,2)*(c(1,2)*f11^2-c(1,2)-c(2,2)))/c(2,2);
f12=0;
f21=0;
F=[f11 f12;f21 f22];
global duf
calCPF(F);
%% ��ɵ���
node1=[1,1];
node2=[0 1];
node3=[1 0];
node4=[0 0];
elem_num=2;
degree=(elem_num+2)*2;
elem1=[1 2 4];
elem2=[4 3 1];
elem1_=[node1 node2 node4];
elem2_=[node4 node3 node1];
elem1_stiff=calstiff(elem1_);
elem2_stiff=calstiff(elem2_);
%% ����ܸ�
global stiffness_matrix
stiffness_matrix=zeros(degree,degree);
assemble(elem1_stiff,elem1);
assemble(elem2_stiff,elem2);
elem_eig=eig(elem1_stiff([2,4],[2,4]))
global_eig=eig(stiffness_matrix([2,4,6],[2,4,6]))

%%
%������ϸն���
function calC()
global c
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
end

%����ƫPƫF
function calCPF(F)
global duf
global c
f11=F(1,1);
f12=F(1,2);
f21=F(2,1);
f22=F(2,2);
duf=zeros(4,4);    %��ʼ��ƫP/ƫF����
duf(1,1)=0.5*c(1,1)*(3*f11^2+f21^2-1)+f12^2*c(3,3)+0.5*c(1,2)*(f12^2+f22^2-1)+c(1,3)*(3*f11*f12+f21*f22);  %����ΪƫP/ƫF����ķ�������
duf(1,2)=f21*f12*c(3,3)+f11*f22*c(1,2)+f11*f21*c(1,3)+f12*f22*c(2,3);
duf(1,3)=(2*f11*f12+f21*f22)*c(3,3)+f11*f12*c(1,2)+0.5*(3*f11^2+f21^2-1)*c(1,3)+0.5*(3*f12^2+f22^2-1)*c(2,3);
duf(1,4)=f11*f21*c(1,1)+f12*f22*c(3,3)+f11*f22*c(1,3)+f21*f12*c(1,3);
duf(2,2)=0.5*(f12^2+3*f22^2-1)*c(2,2)+f21^2*c(3,3)+0.5*(f11^2+f21^2-1)*c(1,2)+(f11*f12+3*f21*f22)*c(2,3);
duf(2,3)=f12*f22*c(2,2)+f11*f21*c(3,3)+(f21*f12+f11*f22)*c(2,3);
duf(2,4)=(f11*f12+2*f21*f22)*c(3,3)+f21*f22*c(1,2)+0.5*(f11^2+3*f21^2-1)*c(1,3)+0.5*(f12^2+3*f22^2-1)*c(2,3);
duf(3,3)=0.5*(3*f12^2+f22^2-1)*c(2,2)+f11^2*c(3,3)+0.5*(f11^2+f21^2-1)*c(1,2)+(3*f11*f12+f21*f22)*c(2,3);
duf(3,4)=f11*f22*c(3,3)+f21*f12*c(1,2)+f21*f11*c(1,3)+f22*f12*c(2,3);
duf(4,4)=0.5*(f11^2+3*f21^2-1)*c(1,1)+f22^2*c(3,3)+0.5*(f12^2+f22^2-1)*c(1,2)+(f11*f12+3*f21*f22)*c(1,3);
duf(2,1)=duf(1,2); 
duf(3,1)=duf(1,3);
duf(4,1)=duf(1,4);
duf(3,2)=duf(2,3);
duf(4,2)=duf(2,4);
duf(4,3)=duf(3,4);
end

%���ϵ���
function assemble(K,jiedian)
global stiffness_matrix
poi1=jiedian(1);
poi2=jiedian(2);
poi3=jiedian(3);
% �Խ�
stiffness_matrix(poi1*2-1:poi1*2,poi1*2-1:poi1*2)=stiffness_matrix(poi1*2-1:poi1*2,poi1*2-1:poi1*2)+K(1:2,1:2);
stiffness_matrix(poi2*2-1:poi2*2,poi2*2-1:poi2*2)=stiffness_matrix(poi2*2-1:poi2*2,poi2*2-1:poi2*2)+K(3:4,3:4);
stiffness_matrix(poi3*2-1:poi3*2,poi3*2-1:poi3*2)=stiffness_matrix(poi3*2-1:poi3*2,poi3*2-1:poi3*2)+K(5:6,5:6);
%�ǶԽ�
stiffness_matrix(poi1*2-1:poi1*2,poi2*2-1:poi2*2)=stiffness_matrix(poi1*2-1:poi1*2,poi2*2-1:poi2*2)+K(1:2,3:4);
stiffness_matrix(poi1*2-1:poi1*2,poi3*2-1:poi3*2)=stiffness_matrix(poi1*2-1:poi1*2,poi3*2-1:poi3*2)+K(1:2,5:6);
stiffness_matrix(poi2*2-1:poi2*2,poi3*2-1:poi3*2)=stiffness_matrix(poi2*2-1:poi2*2,poi3*2-1:poi3*2)+K(3:4,5:6);

stiffness_matrix(poi2*2-1:poi2*2,poi1*2-1:poi1*2)=stiffness_matrix(poi2*2-1:poi2*2,poi1*2-1:poi1*2)+K(3:4,1:2);
stiffness_matrix(poi3*2-1:poi3*2,poi1*2-1:poi1*2)=stiffness_matrix(poi3*2-1:poi3*2,poi1*2-1:poi1*2)+K(5:6,1:2);
stiffness_matrix(poi3*2-1:poi3*2,poi2*2-1:poi2*2)=stiffness_matrix(poi3*2-1:poi3*2,poi2*2-1:poi2*2)+K(5:6,3:4);
end

%���㵥��
function elemstiff=calstiff(x) % x1 y1 x2 y2 x3 y3
global duf
y23=x(4)-x(6);
x32=x(5)-x(3);
y31=x(6)-x(2);
x13=x(1)-x(5);
y12=x(2)-x(4);
x21=x(3)-x(1);
B=[y23 0 y31 0 y12 0;0 x32 0 x13 0 x21;x32 0 x13 0 x21 0;0 y23 0 y31 0 y12]';
elemstiff=1/2*B*duf*B';
end
