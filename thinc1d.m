% ������
clearvars;

% ���W�n�̐���
cor.Lx = 6;% �v�Z�̈�
cor.dLx = 0.02;% �v�f�T�C�Y
cor.x = 0 : cor.dLx : cor.Lx; % ���W�̐���
cor.nx = cor.Lx / cor.dLx + 1;% �v�f��

% �z��̊m��
u = zeros(cor.nx, 1);

% �p�����[�^�[
pt.dt = 0.01;
pt.t_max = 3.5;
pt.iter_max = pt.t_max / pt.dt;
pt.dt_pos = 0.02;% �摜�o�͊Ԋu
u(:, 1) = 1;% �ڗ����x
nu = u(1, 1) * pt.dt / cor.dLx;% �N�[������
beta = 0.2;

% ���������̐ݒ�
xmid = 0.5 * (cor.x(end) + cor.x(1));
C = exp(-20 * (cor.x - xmid).^2);% Gaussian wave
Cex = C;
Cup = C;

for iter = 1 : pt.iter_max
    
    % ���ԍX�V
    pt.iter = iter;
    pt.time = pt.iter * pt.dt;
    
    % ������
    Cex = exp(-20 * (cor.x - xmid - u(1,1) * pt.time).^2);
    
    % 1�����x���㍷���@
    Cup_old = Cup;
    for i = 2 : cor.nx
        Cup(i) = Cup_old(i) - nu * (Cup_old(i) - Cup_old(i - 1));% 1�����x���㍷��
    end
    Cup(1) = Cup(2);
    Cup(end) = Cup(end - 1);
    
    % thinc�@
    C = thinc1d(C, u, beta, pt, cor);
    
    % ���A���^�C������
    filename = ['thinc', ', beta = ', num2str(beta), ', dt = ', num2str(pt.dt), ', dx = ', num2str(cor.dLx),', u = ', num2str(u(1, 1)),', nu = ', num2str(nu(1)),'.gif'];
    visplot(filename, C, Cex, Cup, beta, pt, cor,u, nu)
    
    
end

%% �ȉ��֐�

function [C] = thinc1d(C, u, beta, pt, cor)

% �z��m��

e1 = zeros(cor.nx, 1);
e3 = zeros(cor.nx, 1);
e13 = zeros(cor.nx, 1);
e4 = zeros(cor.nx, 1);
e5 = zeros(cor.nx, 1);
d = zeros(cor.nx - 3, 1);
F = zeros(cor.nx - 3, 1);

for i =  1 : cor.nx - 3 % F��ł̃��[�v
    
    % �����̌����ɉ����ĎQ�Ƃ�������֐��̈ʒu��؂�ւ���B
    if u(i) >= 0
        is = i + 1;
        gamma = 1;
    else
        is = i + 2 ;
        gamma = 0;
    end
    
    % �t���b�N�XF_{i-1/2}�����Ƃ߂�B�iF_{i+1/2}�ł͂Ȃ����Ƃɒ��ӁB�j
    if abs(C(is) - 1) < 10^-5 % C(is)��1�Ȃ��
        F(i) = u(i) * pt.dt;
    elseif  C(is) < 10^-10 % C(is)��0�Ȃ��
        F(i) = 0;
    else
        
        if C(is + 1) >= C(is - 1)
            alpha = 1;
        else
            alpha = -1;
        end
        
        e1(is) = exp(beta * (2 * C(is) - 1)/alpha);
        e3(is) = exp(beta);
        e13(is) = e3(is)*(e3(is)-e1(is))/(e3(is) * e1(is)-1);
        d(is) = 0.5 / beta * log(e13(is));
        e4(is) = cosh(beta * (gamma - u(i) * pt.dt / cor.dLx - d(is)));
        e5(is) = cosh(beta * (gamma - d(is)));
        F(i) = 0.5 * (u(i) * pt.dt + alpha * cor.dLx / beta * log(e5(is) / e4(is)));
        
    end
    
    
end

% �ڗ�������������
for i =  3 : cor.nx - 2% C��ł̃��[�v
    
    Cf = (F(i - 1) - F(i - 2)) / cor.dLx;% �Z��i���甲���o��VOF�l
    C(i) = C(i) - Cf;% �����̗N���o�����͖�������
    
end

% ���E����
C(1) = C(3);
C(2) = C(3);
C(end - 1) = C(end - 2);
C(end) = C(end - 2);

% ����
% �R�}���h�E�B���h�E�ւ̏o��
txt = ['iter = ', num2str(pt.iter), ' / ', num2str(pt.iter_max)];
disp(txt);

end

function [] = visplot(filename, C, Cex, Cup, beta, pt, cor,u, nu)

%plot(cor.x, Cex, cor.x, C)
plot(cor.x, Cup, cor.x, C)

% ���̐ݒ�
title(['time = ', num2str(pt.time, '%.3f')]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
axis equal;
%axis tight; axis on;
fig=gcf;
fig.Position = [100 100 1000 300];
fig.Color='white';
xlim([cor.x(1) cor.x(end)]);
ylim([0 1.2]);
xlabel('x')
ylabel('C')

% �}��
%legtxt = ['thinc', ', beta = ', num2str(beta), ', dt = ', num2str(pt.dt), ', dx = ', num2str(cor.dLx),', u = ', num2str(u(1, 1)),', nu = ', num2str(nu(1))];
%legend({'exact', legtxt},'Location','southwest','FontSize', 10)
legtxt = ['1st order upwind', ', dt = ', num2str(pt.dt), ', dx = ', num2str(cor.dLx),', u = ', num2str(u(1, 1)),', nu = ', num2str(nu(1))];
legtxt2 = ['thinc', ', beta = ', num2str(beta), ', dt = ', num2str(pt.dt), ', dx = ', num2str(cor.dLx),', u = ', num2str(u(1, 1)),', nu = ', num2str(nu(1))];
legend({legtxt, legtxt2},'Location','southwest','FontSize', 10)
legend('boxoff')

frame = getframe(1);
im = frame2im(frame);
[imind, cm] = rgb2ind(im, 256);
if pt.time == pt.dt_pos
    imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.1, 'Loopcount', inf);
elseif rem(pt.time, pt.dt_pos) == 0
    imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.1, 'WriteMode', 'append');
end

end
