# main_vision.m

```matlab
%% global define
global L1;
global Lr;
global L2;
global Lg;
global Lcp;
global gamma1;
global gamma2;
global gamma3;
global L_start;
global config;
global target_handle;
global segs_handle;


%% var define
L1 = 40;
Lr = 20;
L2 = 60; % 20-60
Lg = 60;
Lcp = 60;
gamma1 = 0;
gamma2 = 0;%-3*pi/4;
gamma3 = 0;
L_start = 20;

Lmin = 20;
Lmax = 60 + 40;
radius_lim = 80/pi;


%% init
phi  = 0;
L = L_start - L1 - Lr - L2;
theta1 = 0;
delta1 = 0;
theta2 = 0;%pi/4;
delta2 = 0;
ori_pos = [0;0;0];
ori_rot = eye(3);

config = [phi; L; theta1; delta1; theta2; delta2];
[R0b_g, P0b_g] = forward_kinematics(config);
config_updated = config;
length_updated = [L1; Lr; L2; Lg];
dP = [0;0;0];
dR = 0;

%% setup
figure(2)
view(45, 20);
grid on
axis equal

%% process
while 1
    P0b_g = P0b_g + R0b_g * dP;
    [Lt, theta, delta] = inverse_solver(P0b_g);
    
    if Lt <= L2
        L2_draw = Lt;
        Lr_draw = Lr + L2 - Lt;
        length_updated = [L1; Lr_draw; L2_draw; Lg];
        
        if theta > Lt / radius_lim
            theta = Lt / radius_lim;
        end
        
        theta2 = theta;
        phi = phi + dR;
        delta2 = delta + phi;
        
        config_updated = [phi, Lt - L1 - Lr - L2, theta1, delta1, theta2, delta2];
        disp('topology 1')
    else
        % TODO 当前构型算一次
        [Lr_out, theta_new, ~] = inverse_config_2(P0b_g, theta);
        length_updated = [L1; Lr; L2; Lg];
        
        if theta_new > 60 / radius_lim
            theta_new = 60 / radius_lim;
        end
        
        theta2 = theta_new;
        phi = phi + dR;
        delta2 = delta + phi;
        
        config_updated = [phi, (L2 + Lr_out) - L1 - Lr - L2, theta1, delta1, theta2, delta2];
        disp('topology 2')
    end
    
    [R0b_g, P0b_g] = forward_kinematics(config_updated);
    draw_coordinate_system(30, R0b_g, P0b_g', 'rgb');
    segs_handle = segments_visual_s3(config_updated, ori_pos, ori_rot, length_updated, 1);
    draw_target(P0b_g);
    pause(0.001);
    [dP, dR] = wait_input();
    %[R0b_g, P0b_g] = forward_kinematics(config_updated);
    
    delete_segs_fun(segs_handle);
    delete(target_handle);
    
end

```



# inverse_solver.m

```matlab
function [ L, theta, delta ] = inverse_solver( target_pos )

%% var define
global Lg;

pos_x = target_pos(1);
pos_y = target_pos(2);
pos_z = target_pos(3);

%% calc delta
if abs(pos_x) < 1e-7
    delta = pi / 2;
else
    delta = atan2(-pos_y, pos_x);
end


%% calc theta
ibegin = 0;
iend = pi;
imid = (ibegin + iend) / 2;

max_iter_times = 50;
iter_torl = 0.001;

for i = 0 : max_iter_times
    value = f_theta(imid, target_pos);
    
    if(abs(value) < iter_torl)
        break;
    end
    
    if value < 0
        ibegin = imid;
    else
        iend = imid;
    end
    
    imid = (ibegin + iend) / 2;
end
theta = imid;

%% calc L
if abs(theta) < 1e-7
    L = pos_z - Lg;
else
    L = (pos_z - Lg * cos(theta)) * theta / sin(theta);
end


end
```



# inverse_config_2.m

```matlab
function [ Lr, theta, delta ] = inverse_config_2( target_pos, theta_ori )

%% var define
global Lg;
global L2;

pos_x = target_pos(1);
pos_y = target_pos(2);
pos_z = target_pos(3);

%% calc delta
if abs(pos_x) < 1e-7
    delta = pi / 2;
else
    delta = atan2(-pos_y, pos_x);
end

% %% calc theta
% ibegin = 0;
% iend = pi;
% imid = (ibegin + iend) / 2;
% 
% max_iter_times = 50;
% iter_torl = 0.001;
% 
% for i = 0 : max_iter_times
%     value = f_theta_config2(imid, target_pos);
%     
%     if(abs(value) < iter_torl)
%         break;
%     end
%     
%     if value < 0
%         ibegin = imid;
%     else
%         iend = imid;
%     end
%     
%     imid = (ibegin + iend) / 2;
% end
% theta = imid;

%% calc theta with Newton method
theta_i0 = theta_ori;
theta_i1 = theta_i0 - f_theta_config2(theta_i0, target_pos) / f_dot_theta_config2(theta_i0);
count = 0;
while abs(theta_i1 - theta_i0) > 0.001
    theta_i0 = theta_i1;
    theta_i1 = theta_i0 - f_theta_config2(theta_i0, target_pos) / f_dot_theta_config2(theta_i0);    
    count = count + 1;
    if count >= 100
        break;
    end
end
theta = theta_i1;

if count >= 100
    theta = theta_ori;
end

%% calc Lr
if abs(theta) < 1e-7
    Lr = pos_z - Lg - L2;
else
    Lr = pos_z - Lg * cos(theta) - L2 / theta * sin(theta);
end

end
```



# f_theta.m

```matlab
function [ value ] = f_theta( theta, target_pos )

global Lg;

pos_x = target_pos(1);
pos_y = target_pos(2);
pos_z = target_pos(3);

%% method
left = sqrt(pos_x * pos_x + pos_y * pos_y);
right = (pos_z - Lg * cos(theta)) * tan(theta/2) + Lg * sin(theta);
value = right - left;

end
```



# f_theta_config2.m

```matlab
function [ value ] = f_theta_config2( theta, target_pos )

global Lg;
global L2;

pos_x = target_pos(1);
pos_y = target_pos(2);

left = sqrt(pos_x * pos_x + pos_y * pos_y);
if theta < 1e-7
    right = 0;
else
    right = L2 / theta * (1 - cos(theta)) + Lg * sin(theta);
end
value = right - left;

end
```

