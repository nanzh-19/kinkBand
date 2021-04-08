%% 四边形单元的刚度阵
clear all
clc
f='D:\SIMULIA\Temp\abc.txt';
M=csvread(f);
dof=8;
ele1=zeros(dof,dof);
k=dof/2+1;
for i=1:dof
   if i<=dof/2
       ele1(i,1:dof/2)=M(i,:);
   else 
       ele1(i,1:dof/2)=M(k,:);
       k=k+1;
       ele1(i,dof/2+1:dof)=M(k,:);
       k=k+1;
   end
end
a=diag(ele1);
ele1=ele1+ele1'-diag(a);
elematrix=ele1([4,6,8],[4,6,8]);
DET=det(elematrix)
EIG=eig(elematrix)
%% 三角形单元的刚度阵
clear all
clc
f='D:\SIMULIA\Temp\abc.txt';
M=csvread(f);
dof=6;
ele1=zeros(dof,dof);
k=4+1;
for i=1:dof
   if i<=4
       ele1(i,1:4)=M(i,:);
   else 
       ele1(i,1:4)=M(k,:);
       k=k+1;
       ele1(i,5:dof)=M(k,1:dof-4);
       k=k+1;
   end
end
a=diag(ele1);
ele1=ele1+ele1'-diag(a);
%elematrix=matrix([1,2,4],[1,2,4]);  %节点顺序为213 沿x方向加载
elematrix=ele1([1,2,4],[1,2,4]);
DET=det(elematrix)
EIG=eig(elematrix)
%% 两个三角形单元
clear all
clc
f='D:\SIMULIA\Temp\abc.txt';
M=csvread(f);
dof=6;
ele1=zeros(dof,dof);
ele2=zeros(dof,dof);
k=4+1;
for i=1:dof
   if i<=4
       ele1(i,1:4)=M(i,:);
       ele2(i,1:4)=M(i+8,:);
   else 
       ele1(i,1:4)=M(k,:);
       ele2(i,1:4)=M(k+8,:);
       k=k+1;
       ele1(i,5:dof)=M(k,1:dof-4);
       ele2(i,5:dof)=M(k+8,1:dof-4);
       k=k+1;
   end
end
a=diag(ele1);
ele1=ele1+ele1'-diag(a);
b=diag(ele2);
ele2=ele2+ele2'-diag(b);
matrix=zeros(8,8);
matrix(1:2,1:2)=ele1(1:2,1:2)+ele2(5:6,5:6); %组装总体刚度矩阵
matrix(3:4,3:4)=ele1(3:4,3:4);
matrix(5:6,5:6)=ele2(3:4,3:4);
matrix(7:8,7:8)=ele1(5:6,5:6)+ele2(1:2,1:2);
matrix(1:2,3:4)=ele1(1:2,3:4);
matrix(1:2,5:6)=ele2(5:6,3:4);
matrix(1:2,7:8)=ele1(1:2,5:6)+ele2(5:6,1:2);
matrix(3:4,7:8)=ele1(3:4,5:6);
matrix(5:6,7:8)=ele2(3:4,1:2);
matrix(3:4,1:2)=matrix(1:2,3:4)';
matrix(5:6,1:2)=matrix(1:2,5:6)';
matrix(7:8,1:2)=matrix(1:2,7:8)';
matrix(7:8,3:4)=matrix(3:4,7:8)';
matrix(7:8,5:6)=matrix(5:6,7:8)';
minimatrix=matrix([4,6,8],[4,6,8]);%去除固定约束的行和列
DET=det(minimatrix)
EIG=eig(minimatrix)

