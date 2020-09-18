%% ����ʧ�ȵ�
clc;     %��������д����е�����
clear;   %�������������
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
n=45;    %�ض����ؽǶ���Ӧ�������ת��cita�зִ���
step=30; %�ض����ؽǶ���Ӧ����ط�ֵA������
times=46;%���ؽǶ��зִ���
total=n*times;     %�洢��������Ŀռ��С
e1=zeros(total,1); %����ʧ��ʱ1�����Ӧ�����
e2=zeros(total,1); %����ʧ��ʱ2�����Ӧ�����
energy=zeros(total,1);     %����ʧ��ʱ��Ӧ����
def_gra11=zeros(total,1);  %����ʧ��ʱ�ı����ݶȷ���F11
def_gra22=zeros(total,1);  %����ʧ��ʱ�ı����ݶȷ���F22
jiaodu=zeros(total,1);     %����ʧ��ʱ�ļ��ؽǶ�
c_rotated=zeros(3,3);      %��ʼ��������ת֮��ĸն���
f12=0;                     %�����ݶ�F12=0
f21=0;                     %�����ݶ�F21=0
for k=1:times      %���ؽǶ�ѭ������ѭ��91�Σ���0�ȵ�90��
    angle=(k-1)*0.25*pi/(times-1); %���ؽǶȻ��ȱ�ʾ
    c_rotated(1,1)=c(1,1)*cos(angle)^4+2*(c(1,2)+2*c(3,3))*sin(angle)^2*cos(angle)^2+c(2,2)*sin(angle)^4;  %��ת��ն���ķ���C11
    c_rotated(1,2)=(c(1,1)+c(2,2)-4*c(3,3))*sin(angle)^2*cos(angle)^2+c(1,2)*(sin(angle)^4+cos(angle)^4);  %��ת��ն���ķ���C12
    c_rotated(2,2)=c(1,1)*sin(angle)^4+2*(c(1,2)+2*c(3,3))*sin(angle)^2*cos(angle)^2+c(2,2)*cos(angle)^4;  %��ת��ն���ķ���C22
    c_rotated(1,3)=(c(1,1)-c(1,2)-2*c(3,3))*sin(angle)*cos(angle)^3+(c(1,2)-c(2,2)+2*c(3,3))*sin(angle)^3*cos(angle);  %��ת��ն���ķ���C13
    c_rotated(2,3)=(c(1,1)-c(1,2)-2*c(3,3))*sin(angle)^3*cos(angle)+(c(1,2)-c(2,2)+2*c(3,3))*sin(angle)*cos(angle)^3;  %��ת��ն���ķ���C23
    c_rotated(3,3)=(c(1,1)+c(2,2)-2*c(1,2)-2*c(3,3))*sin(angle)^2*cos(angle)^2+c(3,3)*(sin(angle)^4+cos(angle)^4);     %��ת��ն���ķ���C33
    c_rotated(2,1)=c_rotated(1,2);  %��ת��ն���ķ���C21
    c_rotated(3,1)=c_rotated(1,3);  %��ת��ն���ķ���C31 
    c_rotated(3,2)=c_rotated(2,3);  %��ת��ն���ķ���C32
%     if k==46       %��if�ж���Ϊ���ò�����ת45��ʱӦ�������ת��cita�зִ������
%         n=45;
%     else
%         n=30;
%     end
for i=1:n     %Ӧ�������ת��alphaѭ�����ܹ�ѭ��n��    
%     if k>46   %��if�ж��ǻ���֮ǰ�ó��Ĵ��Խ������Ӧ�������ת��cita���Ż���Ϊ�˼���������ʧ�ȵ�
%        cita=i*pi/n+pi;     %�����ؽǶȴ���45��ʱ��ֻ����3��4���޵�Ӧ�����
%     elseif k<46
%        cita=i*pi/n+pi/2;   %�����ؽǶȴ���45��ʱ��ֻ����2��3���޵�Ӧ�����
%     else
%       cita=i*1.5*pi/n+pi/2;%�����ؽǶȵ���45��ʱ������2��3��4���޵�Ӧ�����
%     end
    cita=pi/2+i*1.5*pi/n;
    for j=1:step   %ѭ��Ӧ����ط�ֵ��ѭ��step��
        A=0.001*j; %Ӧ����ֵ��0.001��ʼ����
        f11=A*cos(cita)+1; %��1�����Ӧ�����ת���ɱ����ݶ�
        f22=A*sin(cita)+1; %��2�����Ӧ�����ת���ɱ����ݶ�      
        duf=zeros(4,4);    %��ʼ��ƫP/ƫF����
        duf(1,1)=0.5*c_rotated(1,1)*(3*f11^2+f21^2-1)+f12^2*c_rotated(3,3)+0.5*c_rotated(1,2)*(f12^2+f22^2-1)+c_rotated(1,3)*(3*f11*f12+f21*f22);  %����ΪƫP/ƫF����ķ�������
        duf(1,2)=f21*f12*c_rotated(3,3)+f11*f22*c_rotated(1,2)+f11*f21*c_rotated(1,3)+f12*f22*c_rotated(2,3);
        duf(1,3)=(2*f11*f12+f21*f22)*c_rotated(3,3)+f11*f12*c_rotated(1,2)+0.5*(3*f11^2+f21^2-1)*c_rotated(1,3)+0.5*(3*f12^2+f22^2-1)*c_rotated(2,3);
        duf(1,4)=f11*f21*c_rotated(1,1)+f12*f22*c_rotated(3,3)+f11*f22*c_rotated(1,3)+f21*f12*c_rotated(1,3);
        duf(2,2)=0.5*(f12^2+3*f22^2-1)*c_rotated(2,2)+f21^2*c_rotated(3,3)+0.5*(f11^2+f21^2-1)*c_rotated(1,2)+(f11*f12+3*f21*f22)*c_rotated(2,3);
        duf(2,3)=f12*f22*c_rotated(2,2)+f11*f21*c_rotated(3,3)+(f21*f12+f11*f22)*c_rotated(2,3);
        duf(2,4)=(f11*f12+2*f21*f22)*c_rotated(3,3)+f21*f22*c_rotated(1,2)+0.5*(f11^2+3*f21^2-1)*c_rotated(1,3)+0.5*(f12^2+3*f22^2-1)*c_rotated(2,3);
        duf(3,3)=0.5*(3*f12^2+f22^2-1)*c_rotated(2,2)+f11^2*c_rotated(3,3)+0.5*(f11^2+f21^2-1)*c_rotated(1,2)+(3*f11*f12+f21*f22)*c_rotated(2,3);
        duf(3,4)=f11*f22*c_rotated(3,3)+f21*f12*c_rotated(1,2)+f21*f11*c_rotated(1,3)+f22*f12*c_rotated(2,3);
        duf(4,4)=0.5*(f11^2+3*f21^2-1)*c_rotated(1,1)+f22^2*c_rotated(3,3)+0.5*(f12^2+f22^2-1)*c_rotated(1,2)+(f11*f12+3*f21*f22)*c_rotated(1,3);
        duf(2,1)=duf(1,2);
        duf(3,1)=duf(1,3);
        duf(4,1)=duf(1,4);
        duf(3,2)=duf(2,3);
        duf(4,2)=duf(2,4);
        duf(4,3)=duf(3,4);
        Cri11=duf(1,1);  %��ƫP/ƫF���������ֵ��ʱ�洢���������Ŵ��㷨��������
        Cri12=duf(1,2);
        Cri13=duf(1,3);
        Cri14=duf(1,4);
        Cri22=duf(2,2);
        Cri23=duf(2,3);
        Cri24=duf(2,4);
        Cri33=duf(3,3);
        Cri34=duf(3,4);
        Cri44=duf(4,4);
        Lb1=[0;0;Cri11;Cri12;Cri13;Cri14;Cri22;Cri23;Cri24;Cri33;Cri34;Cri44]; %�Ŵ��㷨�������½�  
        Ub1=[2*pi;2*pi;Cri11;Cri12;Cri13;Cri14;Cri22;Cri23;Cri24;Cri33;Cri34;Cri44]; %�Ŵ��㷨�������Ͻ�
        options=gaoptimset('Generations',1500,'PopulationSize',50,'TolFun',1e-6);  %�Ŵ��㷨�Ĳ�������
        [xx,AE] = ga(@test,12,[],[],[],[],Lb1,Ub1,[],options);  %�����Ŵ��㷨�����䣬�õ��Ż��⣬����xx���Ż���ĺͱ���ֵ��AEΪ�Ŷ�����    
        if AE<0  %����Ŷ������״γ��ָ����������ǡ��ʧ�ȣ���¼������������ʧ����ͼ
           position=(k-1)*n+i;  %�����������ݴ洢��λ��
           e1(position)=f11-1;  %�洢1�����Ӧ���غ�
           e2(position)=f22-1;  %�洢2�����Ӧ���غ�
           energy(position)=0.125*(f11^2+f21^2-1)^2*c_rotated(1,1)+0.125*(f12^2+f22^2-1)^2*c_rotated(2,2)+0.5*(f11*f12+f21*f22)^2*c_rotated(3,3)+0.25*(f11^2+f21^2-1)*(f12^2+f22^2-1)*c_rotated(1,2)+0.5*(f11^2+f21^2-1)*(f11*f12+f21*f22)*c_rotated(1,3)+0.5*(f11*f12+f21*f22)*(f12^2+f22^2-1)*c_rotated(2,3);%�洢ʧ��ʱ������
           def_gra11(position)=f11;  %�洢ʧ��ʱ��Ӧ���ݶȷ���F11
           def_gra22(position)=f22;  %�洢ʧ��ʱ��Ӧ���ݶȷ���F22
           jiaodu(position)=angle*180/pi; %�洢���ؽǶȣ��Ƕ��ƣ�
           break;
        end               
    end
end
jiaodu(position)
end
%% ȥ����Ч���ݵ�
count=0;  %���³���Ϊ��ȥԭ�г����0ֵ
for i=1:total  %���������ݱ���
    if e1(i)~=0||e2(i)~=0  %�ҳ����еķ���Ӧ�����
        count=count+1;     %����
    end
end
strain1=zeros(count,1);    %1����Ӧ�����
strain2=zeros(count,1);    %2����Ӧ�����
angle=zeros(count,1);      %���ؽǶ�
energyvalue=zeros(count,1);%ʧ��ʱ������ֵ
bianxingtidu11=zeros(count,1);%�����ݶȷ���F11
bianxingtidu22=zeros(count,1);%�����ݶȷ���F22
k=0; %����
for i=1:total  %��ѭ������Ч�����ݴ洢����
    if e1(i)~=0||e2(i)~=0 
        k=k+1;
        strain1(k)=e1(i);
        strain2(k)=e2(i);
        angle(k)=jiaodu(i);
        energyvalue(k)=energy(i);
        bianxingtidu11(k)=def_gra11(i);
        bianxingtidu22(k)=def_gra22(i);
    end
end

