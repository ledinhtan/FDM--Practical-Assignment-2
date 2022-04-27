%% Solve equation -u''(x)=f(x) with the Dirichlet boundary condition 
clear all
close all
clc
%% Initial informations
ax=0.0;
bx=1.0;
ay=0.0;
by=1.0;
cases=2;
Nx=14;
Ny=16;
number_mesh=4;
number_mesh_point=zeros(number_mesh,1);
norm_max=zeros(number_mesh,1);
norm_l2=zeros(number_mesh,1);
norm_maxh1=zeros(number_mesh,1);
norm_h1=zeros(number_mesh,1);
%% Solve discrite solution and refine mesh

for inumber_mesh=1:number_mesh
    
    delta_x=(bx-ax)/Nx;
    delta_y=(by-ay)/Ny;
%% Create mesh point    
    x=zeros(Nx+1,1);
    for i=1:Nx+1
        x(i)=(i-1)*delta_x;
    end
    
    y=zeros(Ny+1,1);
    for j=1:Ny+1
        y(j)=(j-1)*delta_y;
    end
 %% Creat matrix B
    B=sparse(Nx-1,Nx-1);
    for i=1:Nx-1
        if (i==1)
            B(i,i)=2*(delta_x^2+delta_y^2);
            B(i,i+1)=-delta_y^2;
        elseif (i==Nx-1)
            B(i,i)=2*(delta_x^2+delta_y^2);
            B(i,i-1)=-delta_y^2;
        else
            B(i,i-1)=-delta_y^2;
            B(i,i)=2*(delta_x^2+delta_y^2);
            B(i,i+1)=-delta_y^2;
        end
    end
 %% Create matrix A   
    A=sparse((Nx-1)*(Ny-1));
    A(1:(Nx-1),1:(Nx-1))=B;
    A(1:(Nx-1),Nx:2*(Nx-1))=-eye(Nx-1)*delta_x^2;
    
    A((Ny-2)*(Nx-1)+1:(Ny-1)*(Nx-1),(Ny-2)*(Nx-1)+1:(Ny-1)*(Nx-1))=B;
    A((Ny-2)*(Nx-1)+1:(Ny-1)*(Nx-1),(Ny-3)*(Nx-1)+1:(Ny-2)*(Nx-1))=-eye(Nx-1)*delta_x^2;
    
    for i=2:Ny-2
        A((i-1)*(Nx-1)+1:i*(Nx-1),(i-2)*(Nx-1)+1:(i-1)*(Nx-1))=-eye(Nx-1)*delta_x^2;
        A((i-1)*(Nx-1)+1:i*(Nx-1),(i-1)*(Nx-1)+1:i*(Nx-1))=B;
        A((i-1)*(Nx-1)+1:i*(Nx-1),i*(Nx-1)+1:(i+1)*(Nx-1))=-eye(Nx-1)*delta_x^2;
    end
   A = 1/(delta_x^2 * delta_y^2)*A;
%% Create vector F    
    F=zeros((Nx-1)*(Ny-1),1);
    for j=1:Ny-1
        if j==1
            for i=1:Nx-1
                if i==1
                    F((j-1)*(Nx-1)+i)=functionf(x(i+1),y(j+1),cases)+1/delta_y^2*bC(x(i+1),y(j),cases)+1/delta_x^2*bC(x(i),y(j+1),cases);
                elseif i==Nx-1
                    F((j-1)*(Nx-1)+i)=functionf(x(i+1),y(j+1),cases)+1/delta_y^2*bC(x(i+1),y(j),cases)+1/delta_x^2*bC(x(i+2),y(j+1),cases);
                else
                    F((j-1)*(Nx-1)+i)=functionf(x(i+1),y(j+1),cases)+1/delta_y^2*bC(x(i+1),y(j),cases);
                end
            end
        elseif j==Ny-1
            for i=1:Nx-1
                if i==1
                    F((j-1)*(Nx-1)+i)=functionf(x(i+1),y(j),cases)+1/delta_y^2*bC(x(i+1),y(j),cases)+1/delta_x^2*bC(x(i),y(j+1),cases);
                elseif i==Nx-1
                    F((j-1)*(Nx-1)+i)=functionf(x(i+1),y(j+1),cases)+1/delta_y^2*bC(x(i+1),y(j+2),cases)+1/delta_x^2*bC(x(i+2),y(j+1),cases);
                else
                    F((j-1)*(Nx-1)+i)=functionf(x(i+1),y(j+1),cases)+1/delta_y^2*bC(x(i+1),y(j+2),cases);
                end
            end
        else
            for i=1:Nx-1
                if i==1
                    F((j-1)*(Nx-1)+i)=functionf(x(i+1),y(j+1),cases)+1/delta_x^2*bC(x(i),y(j+1),cases);
                elseif i==Nx-1
                    F((j-1)*(Nx-1)+i)=functionf(x(i+1),y(j+1),cases)+1/delta_x^2*bC(x(i+2),y(j+1),cases);
                else
                    F((j-1)*(Nx-1)+i)=functionf(x(i+1),y(j+1),cases);
                end
            end
        end
    end

%% Solve discrete solution
    u=A\F;
%% Get exact solution    
    u_ex=zeros(Nx+1,Ny+1);
    for j=1:Ny+1
        for i=1:Nx+1
            u_ex(i,j)=u_exact(x(i),y(j),cases);
        end
    end
%% Create discrete solution with boundary 
    u_dis=zeros(Nx+1,Ny+1);
    k=1;
    for j=1:Ny+1
        for i=1:Nx+1
            if (i==1||i==Nx+1)
                u_dis(i,j)=u_ex(i,j);
            elseif(j==1||j==Ny+1)
                u_dis(i,j)=u_ex(i,j);
            else
                u_dis(i,j)=u(k);
                k=k+1;
            end
        end
    end
  
    
%% Calculate the error on L^infinity
    norm_max(inumber_mesh)=0.0;
    for i=1:Nx+1
        for j=1:Ny+1
            if (abs(u_dis(i,j)-u_ex(i,j)) > norm_max(inumber_mesh))
                norm_max(inumber_mesh)=abs(u_dis(i,j)-u_ex(i,j));
            end
        end
    end
    
    norm_max(inumber_mesh)
%%  Calculate the error on L^2 

    norm_l2(inumber_mesh)=0;
    for i=1:Nx+1
        for j=1:Ny+1
            norm_l2(inumber_mesh)=norm_l2(inumber_mesh)+(u_dis(i,j)-u_ex(i,j))^2*delta_x*delta_y;
        end
    end
    
    norm_l2(inumber_mesh)=(norm_l2(inumber_mesh))^(1/2);
    norm_l2(inumber_mesh)
%% Calculate the error on maxH1    

    norm_maxh1(inumber_mesh)=0;
    for i=1:Nx
        for j=1:Ny
            if (abs(((u_dis(i+1,j)-u_ex(i+1,j))-(u_dis(i,j)-u_ex(i,j)))/delta_y)+...
                abs(((u_dis(i,j+1)-u_ex(i,j+1))-(u_dis(i,j)-u_ex(i,j)))/delta_x)> norm_maxh1(inumber_mesh))
                norm_maxh1(inumber_mesh)=abs(((u_dis(i+1,j)-u_ex(i+1,j))-(u_dis(i,j)-u_ex(i,j)))/delta_y)+...
                                         abs(((u_dis(i,j+1)-u_ex(i,j+1))-(u_dis(i,j)-u_ex(i,j)))/delta_x);
            end
        end
    end
    norm_maxh1(inumber_mesh)

%% Calculate the error on H1

    norm_h1(inumber_mesh)=0;
    for i=1:Nx
        for j=1:Ny
            norm_h1(inumber_mesh)=norm_h1(inumber_mesh)+(((u_dis(i+1,j)-u_ex(i+1,j))-(u_dis(i,j)-u_ex(i,j)))/delta_y)^2*delta_x*delta_y...
                                                       +(((u_dis(i,j+1)-u_ex(i,j+1))-(u_dis(i,j)-u_ex(i,j)))/delta_x)^2*delta_x*delta_y;
        end
    end
    norm_h1(inumber_mesh)=(norm_h1(inumber_mesh))^(1/2);
%% Figure exact and dicrete solutions    
    figure
    subplot(1,2,1)
    surf(y,x,u_dis)
    title({['solution discrese khi Nx = ',num2str(Nx), ' Ny = ',num2str(Ny)]});
    subplot(1,2,2)
    surf(y,x,u_ex)
    title('solution exact')
    
%% Refine mesh (increse mesh point)    
    Nx=2*Nx;
    Ny=2*Ny;
    number_mesh_point(inumber_mesh)=Ny;
end

%% Figure for errors respect to number of mesh point
figure
plot(log(number_mesh_point), -log(norm_max),'blue', log(number_mesh_point), -log(norm_l2), 'red',...
    log(number_mesh_point), -log(norm_maxh1), 'cyan', log(number_mesh_point), -log(norm_h1), 'magenta',...
    log(number_mesh_point), 2*log(number_mesh_point),'black', log(number_mesh_point), 3/2*log(number_mesh_point),'green',...
    log(number_mesh_point), 2*log(number_mesh_point)+2, 'yellow');xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
title('Errors');
legend('norm_max','norm_l2','norm_maxh1','norm_h1','2x', '3/2x', '2x+2','Location','NorthEastOutside');