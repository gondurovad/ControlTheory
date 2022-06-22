function task2603
% ������ ��������� �������
   
    h = 0.1;
    N = 100;

% ������ ������� �������

    A = [0.0 1.0; 0.0  -1];
    Bu = [0.0; 1.0];
    Bw = [0.0; 1.0];
     
    Q = [1.0,  0.0;
         0.0,  1.0];
    R = 1;
    
% ���������� �������
    
    A_discr = expm(A*h);
    B_discr = integral(@(s) expm(-1*A*s)*Bu, 0,h, 'ArrayValued', true);
    
% ������� X � u, ������� � �����
    
    X_k2 = zeros(2,2);    % ��������� �������, X_N
    k=[];                 % ������ ������������ k: u_(k-1)=k_(k-1)*PHI_(k-1)
   
    for i=N:-1:1
         X_k1=Q+A_discr'*X_k2*A_discr-A_discr'*X_k2*B_discr*(1/(R+B_discr'*X_k2*B_discr))*B_discr';
         k=[k; -inv(R+B_discr'*X_k1*B_discr)*(B_discr'*X_k1*A_discr)];
         X_k2=X_k1;
    end
    
% ���������� ��������� ��������
   
    sigma_w = 1;                                                                              % �� ����� ���������� ��� w
    sigma_v = sigma_w*integral(@(s) expm(A*s)*Bw*Bw'*expm(A'*s), 0,h, 'ArrayValued', true);   % �������� ������� ��� v
    
    E=[0 0];                   % ��� �������� ��� v
    v = mvnrnd(E, sigma_v, N); % ���������� v_k
    
    w=[];
    W = mvnrnd(0, sigma_w, N); % ���������� w_k
    
% ������� ��, � ������ 

    PHI_0=[0.2;0.2];
    PHI_k1=PHI_0;
    PHI=PHI_k1;
    
    for i=2:N
        PHI_k2 = (A_discr+B_discr*k(i, :))*PHI_k1+v(i, :)';
        PHI=[PHI, PHI_k2];
        PHI_k1=PHI_k2;
    end
   
% ������� ����������

    U=[];
    for i=1:N
        U=[U; k(N+1-i, :)*PHI(:, i)];
    end

% ����������� � �������� �������

    t_k2=0;
    t_all=[];
    phi_all=[];
    for i=1:N
        t_k1=t_k2;
        t_k2=t_k2+h;
        [t_tmp, phi_tmp]=ode45(@(t,PH) A*PH+Bu*U(i, :)+Bw*W(i, :), [t_k1 t_k2], PHI(:, i));
        t_all=[t_all; t_tmp];
        phi_all=[phi_all; phi_tmp];
    end
    
% ������ �������

    figure;
    subplot(2,2,1);
        plot(1:1:N, PHI(1,:),'-');
        xlabel('t', 'FontSize', 10);
        ylabel('phi_1(t)', 'FontSize', 10);
    subplot(2,2,2);
        plot(1:1:N, PHI(2,:),'-');
        xlabel('t', 'FontSize', 10);
        ylabel('phi_2(t)', 'FontSize', 10);   
end