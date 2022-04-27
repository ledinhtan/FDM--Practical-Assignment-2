%% Solve equation -u''(x)=f(x) with the Dirichlet boundary condition 
clear all
close all
clc
%% Initial informations
ax=0.0;
bx=1.0;
% ay=0.0;
% by=0.0;
cases=1;
N=10;% number of mesh points of first mesh
number_mesh=4;
number_mesh_point=zeros(number_mesh,1);
norm_max=zeros(number_mesh,1);
norm_l2=zeros(number_mesh,1);
norm_maxh1=zeros(number_mesh,1);
norm_h1=zeros(number_mesh,1);
%% Solve discrite solution and refine mesh

for inumber_mesh=1:number_mesh
    
    delta=(bx-ax)/N;
    
%% Create mesh point    
    x=zeros(N+1,1);
    for i=1:N+1
        x(i)=(i-1)*delta;
    end
    
    y=zeros(N+1,1);
    for j=1:N+1
        y(j)=(j-1)*delta;
    end
 %% Creat matrix B
    B=sparse(N-1,N-1);
    for i=1:N-1
        if (i==1)
            B(i,i)=4;
            B(i,i+1)=-1;
        elseif (i==N-1)
            B(i,i)=4;
            B(i,i-1)=-1;
        else
            B(i,i-1)=-1;
            B(i,i)=4;
            B(i,i+1)=-1;
        end
    end
 %% Creat matrix C
    C=sparse(N-1,N-1);
    for i=1:N-1
        if (i==1)
            C(i,i)=3;
            C(i,i+1)=-1;
        elseif (i==N-1)
            C(i,i)=3;
            C(i,i-1)=-1;
        else
            C(i,i-1)=-1;
            C(i,i)=3;
            C(i,i+1)=-1;
        end
    end

 %% Create matrix A   
    A=sparse((N-1)*(N-1));
    A(1:(N-1),1:(N-1))=B;
    A(1:(N-1),N:2*(N-1))=-eye(N-1);
    
    A((N-2)*(N-1)+1:(N-1)*(N-1),(N-2)*(N-1)+1:(N-1)*(N-1))=C;
    A((N-2)*(N-1)+1:(N-1)*(N-1),(N-3)*(N-1)+1:(N-2)*(N-1))=-eye(N-1);
    
    for i=2:N-2
        A((i-1)*(N-1)+1:i*(N-1),(i-2)*(N-1)+1:(i-1)*(N-1))=-eye(N-1);
        A((i-1)*(N-1)+1:i*(N-1),(i-1)*(N-1)+1:i*(N-1))=B;
        A((i-1)*(N-1)+1:i*(N-1),i*(N-1)+1:(i+1)*(N-1))=-eye(N-1);
    end
   A = 1/(delta^2)*A;
%% Create vector F    
    F=zeros((N-1)*(N-1),1);
    for j=1:N-1
        if(j==1)
            for i =1:N-1
                if(i==1)
                    F((j-1)*(N-1)+i)=functionf(j*delta,i*delta,cases)+(1/delta^2)*bC(j*delta,(i-1)*delta,cases,1)...
                                                                         +(1/delta^2)*bC((j-1)*delta,i*delta,cases,2);
                elseif(i==N-1)
                    F((j-1)*(N-1)+i)=functionf(j*delta,i*delta,cases)+(1/delta^2)*bC(j*delta,(i+1)*delta,cases,3)...
                                                                         +(1/delta^2)*bC((j-1)*delta,i*delta,cases,2);
                else
                    F((j-1)*(N-1)+i)=functionf(j*delta,i*delta,cases)+(1/delta^2)*bC((j-1)*delta,i*delta,cases,2);
                end
            end
        elseif(j==N-1)
            for i =1:N-1
                if(i==1)
                    F((j-1)*(N-1)+i)=functionf(j*delta,i*delta,cases)+(1/delta^2)*bC(j*delta,(i-1)*delta,cases,1)...
                                                                         +(1/delta)*bC((j+1)*delta,i*delta,cases,4);
                elseif(i==N-1)
                    F((j-1)*(N-1)+i)=functionf(j*delta,i*delta,cases)+(1/delta^2)*bC(j*delta,(i+1)*delta,cases,3)...
                                                                         +(1/delta)*bC((j+1)*delta,i*delta,cases,4);
                else
                    F((j-1)*(N-1)+i)=functionf(j*delta,i*delta,cases)+(1/delta)*bC((j+1)*delta,i*delta,cases,4);
                end
            end
        else
            for i=1:N-1
                if(i==1)
                    F((j-1)*(N-1)+i)=functionf(j*delta,i*delta,cases)+(1/delta^2)*bC(j*delta,(i-1)*delta,cases,1);
                elseif(i==N-1)
                    F((j-1)*(N-1)+i)=functionf(j*delta,i*delta,cases)+(1/delta^2)*bC(j*delta,(i+1)*delta,cases,3);
                else
                    F((j-1)*(N-1)+i)=functionf(j*delta,i*delta,cases);
                end
            end
        end
    end

%% Solve discrete solution
    u=A\F;
%% Get exact solution    
    u_ex=zeros(N+1,N+1);
    for j=1:N+1
        for i=1:N+1
            u_ex(i,j)=u_exact(x(i),y(j),cases);
        end
    end
%% Create discrete solution with boundary 
    u_dis=u_ex;
    for i=1:N-1
        for j=1:N-1
            u_dis(i+1,j+1)=u((i-1)*(N-1)+j);
        end
    end
  
    
%% Calculate the error on L^infinity
    norm_max(inumber_mesh)=0.0;
    for i=1:N+1
        for j=1:N+1
            if (abs(u_dis(i,j)-u_ex(i,j)) > norm_max(inumber_mesh))
                norm_max(inumber_mesh)=abs(u_dis(i,j)-u_ex(i,j));
            end
        end
    end
    
    norm_max(inumber_mesh)
%%  Calculate the error on L^2 

    norm_l2(inumber_mesh)=0;
    for i=1:N+1
        for j=1:N+1
            norm_l2(inumber_mesh)=norm_l2(inumber_mesh)+(u_dis(i,j)-u_ex(i,j))^2*delta^2;
        end
    end
    
    norm_l2(inumber_mesh)=(norm_l2(inumber_mesh))^(1/2);
    norm_l2(inumber_mesh)
%% Calculate the error on maxH1    

    norm_maxh1(inumber_mesh)=0;
    for i=1:N
        for j=1:N
            if (abs(((u_dis(i+1,j)-u_ex(i+1,j))-(u_dis(i,j)-u_ex(i,j)))/delta)+...
                abs(((u_dis(i,j+1)-u_ex(i,j+1))-(u_dis(i,j)-u_ex(i,j)))/delta)> norm_maxh1(inumber_mesh))
                norm_maxh1(inumber_mesh)=abs(((u_dis(i+1,j)-u_ex(i+1,j))-(u_dis(i,j)-u_ex(i,j)))/delta)+...
                                         abs(((u_dis(i,j+1)-u_ex(i,j+1))-(u_dis(i,j)-u_ex(i,j)))/delta);
            end
        end
    end
    norm_maxh1(inumber_mesh)

%% Calculate the error on H1

    norm_h1(inumber_mesh)=0;
    for i=1:N
        for j=1:N
            norm_h1(inumber_mesh)=norm_h1(inumber_mesh)+(((u_dis(i+1,j)-u_ex(i+1,j))-(u_dis(i,j)-u_ex(i,j)))/delta)^2*delta...
                                                       +(((u_dis(i,j+1)-u_ex(i,j+1))-(u_dis(i,j)-u_ex(i,j)))/delta)^2*delta;
        end
    end
    norm_h1(inumber_mesh)=(norm_h1(inumber_mesh))^(1/2);
%% Figure exact and dicrete solutions    
    figure
    subplot(1,2,1)
    surf(x,y,u_dis)
    title({['solution discrese khi Nx = ',num2str(N), ' Ny = ',num2str(N)]});
    subplot(1,2,2)
    surf(x,y,u_ex)
    title('solution exact')
    
%% Refine mesh (increse mesh point)    
    N=2*N;
    number_mesh_point(inumber_mesh)=N;
end

%% Figure for errors respect to number of mesh point
figure
plot(log(number_mesh_point), -log(norm_max),'blue', log(number_mesh_point), -log(norm_l2), 'red',...
    log(number_mesh_point), -log(norm_maxh1), 'cyan', log(number_mesh_point), -log(norm_h1), 'magenta',...
    log(number_mesh_point), 2*log(number_mesh_point),'black', log(number_mesh_point), 3/2*log(number_mesh_point),'green',...
    log(number_mesh_point), 2*log(number_mesh_point)+2, 'yellow'); % cases=1
% plot(log(number_mesh_point), -log(norm_max),'blue', log(number_mesh_point), -log(norm_l2), 'red',...
%     log(number_mesh_point), -log(norm_maxh1), 'cyan', log(number_mesh_point), -log(norm_h1), 'magenta',...
%     log(number_mesh_point), -2*log(number_mesh_point),'black', log(number_mesh_point), -3/2*log(number_mesh_point),'green',...
%     log(number_mesh_point), -2*log(number_mesh_point)-2, 'yellow'); % cases=2 and cases=3

xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
title('Errors');
legend('norm_max','norm_l2','norm_maxh1','norm_h1','2x', '3/2x', '2x+2','Location','NorthEastOutside');  
