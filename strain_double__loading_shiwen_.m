%% 计算失稳点
clc;     %清除命令行窗口中的数据
clear;   %清楚工作区数据
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
n=45;    %特定加载角度下应变加载旋转角cita切分次数
step=30; %特定加载角度下应变加载幅值A步长数
times=46;%加载角度切分次数
total=n*times;     %存储数据所需的空间大小
e1=zeros(total,1); %材料失稳时1方向的应变加载
e2=zeros(total,1); %材料失稳时2方向的应变加载
energy=zeros(total,1);     %材料失稳时的应变能
def_gra11=zeros(total,1);  %材料失稳时的变形梯度分量F11
def_gra22=zeros(total,1);  %材料失稳时的变形梯度分量F22
jiaodu=zeros(total,1);     %材料失稳时的加载角度
c_rotated=zeros(3,3);      %初始化材料旋转之后的刚度阵
f12=0;                     %变形梯度F12=0
f21=0;                     %变形梯度F21=0
for k=1:times      %加载角度循环，从循环91次，从0度到90度
    angle=(k-1)*0.25*pi/(times-1); %加载角度弧度表示
    c_rotated(1,1)=c(1,1)*cos(angle)^4+2*(c(1,2)+2*c(3,3))*sin(angle)^2*cos(angle)^2+c(2,2)*sin(angle)^4;  %旋转后刚度阵的分量C11
    c_rotated(1,2)=(c(1,1)+c(2,2)-4*c(3,3))*sin(angle)^2*cos(angle)^2+c(1,2)*(sin(angle)^4+cos(angle)^4);  %旋转后刚度阵的分量C12
    c_rotated(2,2)=c(1,1)*sin(angle)^4+2*(c(1,2)+2*c(3,3))*sin(angle)^2*cos(angle)^2+c(2,2)*cos(angle)^4;  %旋转后刚度阵的分量C22
    c_rotated(1,3)=(c(1,1)-c(1,2)-2*c(3,3))*sin(angle)*cos(angle)^3+(c(1,2)-c(2,2)+2*c(3,3))*sin(angle)^3*cos(angle);  %旋转后刚度阵的分量C13
    c_rotated(2,3)=(c(1,1)-c(1,2)-2*c(3,3))*sin(angle)^3*cos(angle)+(c(1,2)-c(2,2)+2*c(3,3))*sin(angle)*cos(angle)^3;  %旋转后刚度阵的分量C23
    c_rotated(3,3)=(c(1,1)+c(2,2)-2*c(1,2)-2*c(3,3))*sin(angle)^2*cos(angle)^2+c(3,3)*(sin(angle)^4+cos(angle)^4);     %旋转后刚度阵的分量C33
    c_rotated(2,1)=c_rotated(1,2);  %旋转后刚度阵的分量C21
    c_rotated(3,1)=c_rotated(1,3);  %旋转后刚度阵的分量C31 
    c_rotated(3,2)=c_rotated(2,3);  %旋转后刚度阵的分量C32
%     if k==46       %该if判断是为了让材料旋转45度时应变加载旋转角cita切分次数变多
%         n=45;
%     else
%         n=30;
%     end
for i=1:n     %应变加载旋转角alpha循环，总共循环n次    
%     if k>46   %该if判断是基于之前得出的粗略结果，对应变加载旋转角cita的优化，为了计算出更多的失稳点
%        cita=i*pi/n+pi;     %当加载角度大于45度时，只计算3、4象限的应变加载
%     elseif k<46
%        cita=i*pi/n+pi/2;   %当加载角度大于45度时，只计算2、3象限的应变加载
%     else
%       cita=i*1.5*pi/n+pi/2;%当加载角度等于45度时，计算2、3、4象限的应变加载
%     end
    cita=pi/2+i*1.5*pi/n;
    for j=1:step   %循环应变加载幅值，循环step次
        A=0.001*j; %应力幅值从0.001开始增长
        f11=A*cos(cita)+1; %将1方向的应变加载转换成变形梯度
        f22=A*sin(cita)+1; %将2方向的应变加载转换成变形梯度      
        duf=zeros(4,4);    %初始化偏P/偏F矩阵
        duf(1,1)=0.5*c_rotated(1,1)*(3*f11^2+f21^2-1)+f12^2*c_rotated(3,3)+0.5*c_rotated(1,2)*(f12^2+f22^2-1)+c_rotated(1,3)*(3*f11*f12+f21*f22);  %以下为偏P/偏F矩阵的分量计算
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
        Cri11=duf(1,1);  %将偏P/偏F矩阵分量的值临时存储，用来做遗传算法的输入量
        Cri12=duf(1,2);
        Cri13=duf(1,3);
        Cri14=duf(1,4);
        Cri22=duf(2,2);
        Cri23=duf(2,3);
        Cri24=duf(2,4);
        Cri33=duf(3,3);
        Cri34=duf(3,4);
        Cri44=duf(4,4);
        Lb1=[0;0;Cri11;Cri12;Cri13;Cri14;Cri22;Cri23;Cri24;Cri33;Cri34;Cri44]; %遗传算法变量的下界  
        Ub1=[2*pi;2*pi;Cri11;Cri12;Cri13;Cri14;Cri22;Cri23;Cri24;Cri33;Cri34;Cri44]; %遗传算法变量的上界
        options=gaoptimset('Generations',1500,'PopulationSize',50,'TolFun',1e-6);  %遗传算法的参数设置
        [xx,AE] = ga(@test,12,[],[],[],[],Lb1,Ub1,[],options);  %调用遗传算法工具箱，得到优化解，其中xx是优化后的和变量值，AE为扰动能量    
        if AE<0  %如果扰动能量首次出现负数，则材料恰好失稳，记录以下数据来画失稳相图
           position=(k-1)*n+i;  %计算以下数据存储的位置
           e1(position)=f11-1;  %存储1方向的应变载荷
           e2(position)=f22-1;  %存储2方向的应变载荷
           energy(position)=0.125*(f11^2+f21^2-1)^2*c_rotated(1,1)+0.125*(f12^2+f22^2-1)^2*c_rotated(2,2)+0.5*(f11*f12+f21*f22)^2*c_rotated(3,3)+0.25*(f11^2+f21^2-1)*(f12^2+f22^2-1)*c_rotated(1,2)+0.5*(f11^2+f21^2-1)*(f11*f12+f21*f22)*c_rotated(1,3)+0.5*(f11*f12+f21*f22)*(f12^2+f22^2-1)*c_rotated(2,3);%存储失稳时的能量
           def_gra11(position)=f11;  %存储失稳时的应变梯度分量F11
           def_gra22(position)=f22;  %存储失稳时的应变梯度分量F22
           jiaodu(position)=angle*180/pi; %存储加载角度（角度制）
           break;
        end               
    end
end
jiaodu(position)
end
%% 去除无效数据点
count=0;  %以下程序为除去原有程序的0值
for i=1:total  %将所有数据遍历
    if e1(i)~=0||e2(i)~=0  %找出其中的非零应变加载
        count=count+1;     %计数
    end
end
strain1=zeros(count,1);    %1方向应变加载
strain2=zeros(count,1);    %2方向应变加载
angle=zeros(count,1);      %加载角度
energyvalue=zeros(count,1);%失稳时的能量值
bianxingtidu11=zeros(count,1);%变形梯度分量F11
bianxingtidu22=zeros(count,1);%变形梯度分量F22
k=0; %计数
for i=1:total  %该循环将有效的数据存储起来
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

