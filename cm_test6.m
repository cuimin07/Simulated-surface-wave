clc; clear all

snapshot_time = 150; % time--ms

tic 
%%% 参数设置 %%%
N = 200;
fm = 28;   % 中心频率
space_x = 2; space_z = 2;  % 网格步长
sample_t = 0.1/1000;  % 时间步长0.1ms
T = 300/1000;    % 总时间为250ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%不调%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = T/sample_t;   % 外层时间循环总次数:250/0.1 = 2500次
layer = 20*space_x;  % 吸收层的厚度==40
coefficientR = 0.0001;  % 理想的反射系数：基本不反射

C1 = 1.125;
C2 = -0.04166667;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 均匀半空间模型的初始化 %%%
% 预分配提速:
VP = zeros(200,200)+1300;   %纵波速度
VS = zeros(200,200)+800;    %横波速度
Density = zeros(200,200)+2700;  % 介质的密度
lame1 = zeros(200,200);
lame2 = zeros(200,200);
for i = 1:N
    j = 1:N;  % 矢量化加速
    lame1(i,j) = Density(i,j) .* (VP(i,j).^2 - 2*VS(i,j).^2);
    lame2(i,j) = Density(i,j) .* VS(i,j).^2;
end

%%% U、V、P、Q、R：速度、应力全空间初始化 %%%U:vx; V:vz; P:Pxx; Q:Pzz; R:Pxz
% 预分配:全x向参数
U_x = zeros(200,200);
V_x = zeros(200,200);
P_x = zeros(200,200);
Q_x = zeros(200,200);
R_x = zeros(200,200);
% 预分配:全y向的参数
U_y = zeros(200,200);
V_y = zeros(200,200);
P_y = zeros(200,200);
Q_y = zeros(200,200);
R_y = zeros(200,200);
% 预分配:总参数
U = zeros(200,200);
V = zeros(200,200);
P = zeros(200,200);
Q = zeros(200,200);
R = zeros(200,200);
% 预分配:间距
dx = zeros(1,200);
dy = zeros(1,200);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%震源剖面%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xxx = zeros(5000,200);
zzz = zeros(5000,200);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% AEA 自由表面边界条件 %%%
% 其实就是把原始空间的前4行拿出来，用来做自由空间层！
% 所以差分计算时，i和j都是从4开始
for i = 1:4
	j = 1:N  % 矢量化加速
	Density(i,j) = 0.5*Density(i,j);
	lame1(i,j) = 0;
	lame2(i,j) = lame2(i,j);
end

%%% 开始有限差分的计算 %%%
for k = 1:K  % 时间是最外层循环：相当于总的迭代次数/传播的时间
    fprintf('迭代计算开始,当前是第%d次"时间"大循环!\n',k)
	mibinbin = k;
	% AEA 自由边界设置：
	for i = 1:N
		Q(i,4) = 0;
	end

	% 差分计算：各点的速度更新
	for i = 4:(N-4)
		for j = 4:(N-4)
			% 区域边界的设置
			if i <= 4 + layer/space_x
				dx(i) = (2*VP(i,j)/layer)*log(1/coefficientR)*((4+layer/space_x-i)*space_x/layer)^4;
			end
			if i >= N-4-layer/space_x
				dx(i) = (2*VP(i,j)/layer)*log(1/coefficientR)*((i-(N-4-layer/space_x))*space_x/layer)^4;
			end
			if j >= N-4-layer/space_x  %  400 - 4 - 10/0.5 = 346  每一个i下，j要在346之后才会不等于0
                %fprintf('第%d次"行循环"已进入j>=346阶段,当前:j=%d\n',i,j);
                %pause(1);
				dy(j) = (2*VP(i,j)/layer)*log(1/coefficientR)*((j-(N-4-layer/space_x))*space_x/layer)^4;
			end
			% 速度迭代计算
			U_x(i,j) = ( (1-0.5*sample_t*dx(i))/(1+0.5*sample_t*dx(i)))*U_x(i,j) + ...
			           (1/(1+0.5*sample_t*dx(i)))*(1/Density(i,j))*(sample_t/space_x)*...
			           ( C1*(P(i+1,j)-P((i),j))+C2*( P((i+2),j)-P((i-1),j)));  % 已检查
			U_y(i,j) = ((1-0.5*sample_t*dy(j))/(1+0.5*sample_t*dy(j)))*U_y(i,j) + ...
			           (1/(1+0.5*sample_t*dy(j)))*(1/Density(i,j))*(sample_t/space_z)* ...
			           ( C1*(R(i,j)-R(i,j-1))+C2*(R(i,(j+1))-R(i,(j-2))));
			U(i,j) = U_x(i,j) + U_y(i,j); 

			V_x(i,j) = ((1-0.5*sample_t*dx(i))/(1+0.5*sample_t*dx(i)))*V_x(i,j) + ...
			           (1/(1+0.5*sample_t*dx(i)))*(1/Density(i,j))*(sample_t/space_x)*...
			           ( C1*(R(i,j)-R(i-1,j))+C2*(R((i+1),j)-R(i-2,j)));
			V_y(i,j) = ((1-0.5*sample_t*dy(j))/(1+0.5*sample_t*dy(j)))*V_y(i,j) + ...
			           (1/(1+0.5*sample_t*dy(j)))*(1/Density(i,j))*(sample_t/space_z)* ...
			           ( C1*(Q(i,(j+1))-Q(i,j))+C2*( Q(i,(j+2))-Q(i,(j-1))));
			V(i,j) = V_x(i,j) + V_y(i,j);
		end
	end

    % 差分计算：各点的应力更新
    for i = 4:(N-4)
    	for j = 4:(N-4)
    		if i <= 4 + layer/space_x
    			dx(i) = (2*VP(i,j)/layer)*log(1/coefficientR)*((4+layer/space_x-i)*space_x/layer)^4;
    		end
    		if i >= N-4-layer/space_x
    			dx(i) = (2*VP(i,j)/layer)*log(1/coefficientR)*((i-(N-4-layer/space_x))*space_x/layer)^4;
    		end
    		if j >= N-4-layer/space_x
    			dy(j) = (2*VP(i,j)/layer)*log(1/coefficientR)*((j-(N-4-layer/space_x))*space_x/layer)^4;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    		% 中心震源的位置:
    		V(100,5) = 3*(1-2*(pi*fm*(k-1000)*sample_t)^2)*...
    		           exp(-(pi*fm*(k-1000)*sample_t)^2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    		% 震源的设置:
    		P_x(i,j) = ((1-0.5*sample_t*dx(i))/(1+0.5*sample_t*dx(i)))*P_x(i,j) + ...
    		           (1/(1+0.5*sample_t*dx(i)))*(lame1(i,j)+2*lame2(i,j))*(sample_t/space_x)*...
    		           (C1*(U((i),j)-U(i-1,j))+C2*(U((i+1),j)-U((i-2),j)));
    		P_y(i,j) = ((1-0.5*sample_t*dy(j))/(1+0.5*sample_t*dy(j)))*P_y(i,j) + ...
    		           (1/(1+0.5*sample_t*dy(j)))*lame1(i,j)*(sample_t/space_z)*...
    		           (C1*(V(i,j)-V(i,(j-1))) + C2*( V(i,(j+1))-V(i,(j-2))));
    		P(i,j) = P_x(i,j) + P_y(i,j);

    		Q_x(i,j) = ((1-0.5*sample_t*dx(i))/(1+0.5*sample_t*dx(i)))*Q_x(i,j) + ...
    		           (1/(1+0.5*sample_t*dx(i)))*lame1(i,j)*(sample_t/space_x)*...
    		           (C1*(U((i),j)-U(i-1,j))+C2*(U((i+1),j)-U((i-2),j)));
    		Q_y(i,j) = ((1-0.5*sample_t*dy(j))/(1+0.5*sample_t*dy(j)))*Q_y(i,j) + ...
    		           (1/(1+0.5*sample_t*dy(j)))*(lame1(i,j)+2*lame2(i,j))*(sample_t/space_z)*...
                       (C1*(V(i,j)-V(i,(j-1)))+C2*(V(i,j+1)-V(i,(j-2))));
            Q(i,j) = Q_x(i,j) + Q_y(i,j);

            R_x(i,j) = ((1-0.5*sample_t*dx(i))/(1+0.5*sample_t*dx(i)))*R_x(i,j) + ...
                       (1/(1+0.5*sample_t*dx(i)))*lame2(i,j)*(sample_t/space_x)*...
                       (C1*(V(i+1,j)-V((i),j)) + C2*( V((i+2),j)-V((i-1),j)));
            R_y(i,j) = ((1-0.5*sample_t*dy(j))/(1+0.5*sample_t*dy(j)))*R_y(i,j) + ...
                       (1/(1+0.5*sample_t*dy(j)))*lame2(i,j)*(sample_t/space_z)*...
                       (C1*(U(i,(j+1))-U(i,j)) + C2*(U(i,(j+2))-U(i,(j-1))));
            R(i,j) = R_x(i,j) + R_y(i,j);
        end
    end
   
    xxx(k,:)=U(:,5)';
    zzz(k,:)=V(:,5)';
end
toc

figure( 'position',[200 150 370 400]),pcolor(zzz(500:2500,:)),shading interp,axis ij
title('Z方向模拟记录');
set(gca,'FontSize',10);
xlabel('Trace Number','fontsize',12),ylabel('Time(ms)','fontsize',12);
set(gca,'xaxislocation','top','box','off','TickDir','out');
colorbar
hcb=colorbar
title(hcb,' ','fontsize',12)
colormap('gray');

% figure( 'position',[200 150 370 350]),pcolor(xxx(500:2500,:)),shading interp,axis ij
% title('X方向模拟记录');
% set(gca,'FontSize',10);
% xlabel('Trace Number','fontsize',12),ylabel('Time(ms)','fontsize',12);
% set(gca,'xaxislocation','top','box','off','TickDir','out');
% colorbar
% hcb=colorbar
% title(hcb,' ','fontsize',12)
% colormap('gray');