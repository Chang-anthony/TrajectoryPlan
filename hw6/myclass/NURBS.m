classdef NURBS
    %NURBS Summary of this class goes here
    %Detailed explanation goes here
    properties
    %given
        % degree (order - 1 => in this program is (k-1)) 
        % 階數 (為階數-1 => 為(k-1))
        p
        % order (order of NURBS curves)
        % NURBS 曲線之級數 
        k
        % control points (Determine the position of the curve)
        % 控制點 決定曲線的位置，通常不會在曲線上 
        P 
        % weight (Affects the trend of the curve towards the control point)
        % 權重 影響曲線朝向控制點的趨勢
        W
        % knot vector (Variation of the curve that determines the design with respect to the parameter u)
        % 節點向量 決定設計曲線對於參數u的變化
        % number of konts = n + k + 1 =>(order k + control point n+1)
        knot

    %calculate
        % basis funtion (基底函數)
        N
        % single rational B-spline
        % 單一有理式B型曲線
        R
        % Same knot
        Same_knot
    end
    
    methods
        function NURBS = NURBS(knots,control_points,weight,order)
            %NURBS Construct an instance of this class
            %   Detailed explanation goes here
            NURBS.p = order-1;
            NURBS.k = order;
            NURBS.P = control_points;
            NURBS.W = weight;
            NURBS.knot = knots;
            if knots(end) > 1
                NURBS.knot = NURBS.normalized();
            end
%             NURBS.Same_knot = NURBS.Get_Same_knot();
              NURBS.Same_knot = NURBS.Create_Same_knot();
        end

        function [Px,Py] = Get_Control_points_value(NURBS,SampleSize)
            %SampleSize is parameter u size
            %C(u) =  i=0 -> n sum( Pi*Ri,k(u))
            %R(i,k)(u) = i=0 -> n sum(N(i,k)(u)*Wi)/j=0 -> n sum(N(j,k)(u)*wj)
            n = length(NURBS.P);
            Px = [];
            Py = [];
            u = 1/SampleSize;

           for j = 1:SampleSize+1
                x = 0;
                y = 0;
                for i = 1:n                   
                    R_ik = NURBS.Get_Rational_basis(NURBS.k,i,u);
                    x = x + R_ik*NURBS.P(i,1);
                    y = y + R_ik*NURBS.P(i,2);
                end
                Px = [Px,x];
                Py = [Py,y];
                u = 1/SampleSize*j;
           end
        end

        function [Px,Py] = Get_Control_points_value_param_list(NURBS,param_list)
            %SampleSize is parameter u size
            %C(u) =  i=0 -> n sum( Pi*Ri,k(u))
            %R(i,k)(u) = i=0 -> n sum(N(i,k)(u)*Wi)/j=0 -> n sum(N(j,k)(u)*wj)
           n = length(NURBS.P);
           Px = [];
           Py = [];

           for j = 1:length(param_list)
                u = param_list(j);
                x = 0;
                y = 0;
                for i = 1:n                   
                    R_ik = NURBS.Get_Rational_basis(NURBS.k,i,u);
                    x = x + R_ik*NURBS.P(i,1);
                    y = y + R_ik*NURBS.P(i,2);
                end
                Px = [Px,x];
                Py = [Py,y];  
           end
        end
  
        function [dPx,dPy] = Get_Derivatives_Control_points(NURBS,k,u,m)
           %C(u) =  i=0 -> n sum( Pi*Ri,k(u))
           %R(i,k)(u) = i=0 -> n sum(N(i,k)(u)*Wi)/j=0 -> n sum(N(j,k)(u)*wj)
           x = 0;
           y = 0;
           n = length(NURBS.P);
           for i = 1:n
               R_ik = NURBS.Get_Derivatives_Rbasis(k,i,u,m);
               x = x + R_ik*NURBS.P(i,1);
               y = y + R_ik*NURBS.P(i,2);
           end
           dPx = x;
           dPy = y;
        end

        function R = Get_Rational_basis(NURBS,k,i,u)
           R_d = NURBS.Get_Rationnal_basis_denominator(k,u);
           R = NURBS.Get_Basis(k,i,u) * NURBS.W(i)/R_d;
        end

        function R_d = Get_Rationnal_basis_denominator(NURBS,k,u)
           w_length = length(NURBS.W);% w length = control points length
           R_d = 0;

           for w = 1:w_length
               R_d = R_d + NURBS.W(w)*NURBS.Get_Basis(k,w,u);
           end

        end

        function dR = Get_Derivatives_Rbasis(NURBS,k,i,u,m)
            if NURBS.Get_Derivatives_basis(k,i,u,m) ~= false
                dR_d = NURBS.Get_Derivatives_Rational_basis_denominator(k,u,m);
                if dR_d == 0
                   dR = false;
                else
                   dR = NURBS.Get_Derivatives_basis(k,i,u,m)* NURBS.W(i)/dR_d;
                end
            else
                dR = false;
            end
        end

        function dR_d = Get_Derivatives_Rational_basis_denominator(NURBS,k,u,m)
           w_length = length(NURBS.W);% w length = control points length
           dR_d = 0;

           for w = 1:w_length
               dR_d = dR_d + NURBS.W(w)*NURBS.Get_Derivatives_basis(k,w,u,m);
           end

        end

        function R_basis = Create_Rational_order_functions(NURBS,k,SampleSize)
            %k is now order
            %SampleSize is parameter u size
            %R(i,k)(u) = N(i,k)(u)*Wi/(j = 0 -> n (n = control points length))N(j,k)*Wj
            w_length = length(NURBS.W);% w length = control points length
            n_length = length(NURBS.knot);% n+1 length = knot vector length
            R_basis = zeros(n_length-k+1,SampleSize);
            
            for i = 1:w_length
                for j = 1:SampleSize
                   u = 1/SampleSize * j;
                   temp = NURBS.Get_Rationnal_basis_denominator(k,u);
                   R_basis(i,j) = NURBS.Get_Basis(k,i,u) * NURBS.W(i)/temp;
                end
            end  
        end

        function R_basis = Create_Rational_order_param_list(NURBS,k,u_list)
            %k is now order
            %length(u_list) = Sample Size 
            %R(i,k)(u) = N(i,k)(u)*Wi/(j = 0 -> n (n = control points length))N(j,k)*Wj
            w_length = length(NURBS.W);% w length = control points length
            n_length = length(NURBS.knot);% n+1 length = knot vector length
            R_basis = zeros(n_length-k+1,length(u_list));

            for i = 1:w_length
                for j = 1:length(u_list)
                   u = u_list(j);
                   R_d = NURBS.Get_Rationnal_basis_denominator(k,u);
                   R_basis(i,j) = NURBS.Get_Basis(k,i,u) * NURBS.W(i)/R_d;
                end
            end 
        end

        function R_basis = Create_Rational_functions(NURBS,SampleSize)
            %SampleSize is parameter u size
            %if k = 1 will return R0,k ... R1,k...R2,k...R(n+1),k
            n_length = length(NURBS.knot); % knot points length n+1
            R_basis = zeros(NURBS.k,n_length,SampleSize);
             for order = 1:NURBS.k
                 for i = 1:n_length-order+1
                     for j = 1:SampleSize
                         u = 1/SampleSize *j;
                         R_basis(order,i,j) = NURBS.Get_Rational_basis(order,i,u); 
                     end
                 end
             end
        end

        function Create_RBasis_description(NURBS,SampleSize)
            NURBS.R = NURBS.Create_Rational_functions(SampleSize);
        end

        function Create_NBasis_description(NURBS,SampleSize)
            NURBS.N = NURBS.Create_Basis_funtions(SampleSize);
        end
      
        function N = Get_Basis(NURBS,k,i,u)
            %k is now order
            %u is paremeter we generator
            %i is current knot position
            %N(i,k)(u) = u-u(i)/u(i+k-1)*N(i,k-1)(u) + u(i+k)-u/u(i+k)-u(i+1)*N(i+1,k-1)
            knot_i = NURBS.knot(i);
            knot_i_1 = NURBS.knot(i+1);
            knot_i_k_1 = NURBS.knot(k+i-1);
            knot_i_k = NURBS.knot(k+i);
            if k == 1 %order = 1
                if knot_i <= u && u <= knot_i_1
                    N = 1;
                else
                    N = 0;
                end
            else
                if (knot_i_k_1 - knot_i == 0)
                    cofe1 = 0;
                else
                    cofe1 = (u-knot_i)/(knot_i_k_1-knot_i);
                end
                if (knot_i_k - knot_i_1 == 0)
                    cofe2 = 0;
                else
                    cofe2 = (knot_i_k - u)/(knot_i_k - knot_i_1);
                end
                N = cofe1 * NURBS.Get_Basis(k-1,i,u) + cofe2 * NURBS.Get_Basis(k-1,i+1,u);
            end
        end

        function dN = Get_Derivatives_basis(NURBS,k,i,u,m)
            %k is now order
            %u is paremeter we generator
            %i is current knot position
            %m is Derivatives times
            %dN = (k-1)*[N(i,k-1)'(m-1)(u)/u(i+k-1)-u(i) - N(i+1,k-1)'(m-1)(u)/u(i+k)-u(i)]
            if NURBS.Can_Derivatives(u) == true
                knot_i = NURBS.knot(i);
                knot_i_1 = NURBS.knot(i+1);
                knot_i_k_1 = NURBS.knot(k+i-1);
                knot_i_k = NURBS.knot(k+i);
                if m >= 1
                    temp_N = NURBS.Get_Derivatives_basis(k-1,i,u,m-1);
                    temp_N_1 = NURBS.Get_Derivatives_basis(k-1,i+1,u,m-1);
                    if (knot_i_k_1 - knot_i) == 0
                        cofe1 =0;
                    else
                        cofe1 = (temp_N * (k-1))/(knot_i_k_1 - knot_i);
                    end
                    if (knot_i_k - knot_i_1) == 0
                        cofe2 = 0;
                    else
                        cofe2 = (temp_N_1*(k-1))/(knot_i_k - knot_i_1);
                    end
                    dN = cofe1 - cofe2;
                elseif m == 0
                     dN = NURBS.Get_Basis(k,i,u);
                else
                    disp("some thing went wrong in Get_Derivatives_basis");
                    dN = false;
                end
            else
                disp("this parameter u can't Derivative dN return false");
                dN = false;
            end
        end

        function N_basis = Create_Basis_funtions(NURBS,SampleSize)
         %SampleSize is parameter u size
         %if k = 1 will return N0,k ... N1,k...N2,k...N(n+1),k
         %return All Basis functions 
         n_length = length(NURBS.knot); % knot points length n+1
         N_basis = zeros(NURBS.k,n_length,SampleSize);
            for order = 1:NURBS.k
                for j = 1:(n_length-order)
                    for i = 1:SampleSize
                        u = 1/SampleSize * i;
                        N_basis(order,j,i) = NURBS.Get_Basis(order,j,u);
                    end
                end
            end
        end

        function N_basis = Create_Basis_order_funtions(NURBS,k,SampleSize)
         %SampleSize is parameter u size
         %k is now we want to create basis functions 
         %if k = 1 will return N0,k ... N1,k...N2,k...N(n+1),k
         n_length = length(NURBS.knot); % knot points length n+1
         %basis functions size is knot length + k-1 
         %if k = 1  Nbasis = size(knot_length,SampleSize)
         %if k = 2  Nbasis = size(knot_length-1,SampleSize)
         N_basis = zeros(n_length-k+1,SampleSize);
             for j = 1:(n_length-k)
                 for i = 1:SampleSize
                     u = 1/SampleSize * i;
                     N_basis(j,i) = NURBS.Get_Basis(k,j,u);
                 end
             end
        end

        function N_basis = Create_Basis_order_functions_paramlist(NURBS,k,u_list)
         %length(u_list) is sample size
         %k is now we want to create basis functions 
         %if k = 1 will return N0,k ... N1,k...N2,k...N(n+1),k
         n_length = length(NURBS.knot); % knot points length n+1
         %basis functions size is knot length + k-1 
         %if k = 1  Nbasis = size(knot_length,SampleSize)
         %if k = 2  Nbasis = size(knot_length-1,SampleSize)
         N_basis = zeros(n_length-k+1,length(u_list));
             for j = 1:(n_length-k)
                 for i = 1:length(u_list)
                     u = u_list(i);
                     N_basis(j,i) = NURBS.Get_Basis(k,j,u);
                 end
             end
        end

        function param = Get_Interpolator_param_list(NURBS,SLine,T,u)
          param = [];
          param = [param,u];
          while param(end) < 1 
              V = SLine.Get_normlized_time_vel_value(u);
              if NURBS.Can_Derivatives(u) ~= false
                 [dx,dy] = NURBS.Get_Derivatives_Control_points(NURBS.k,u,1);
                 du = V/sqrt(dx^2+dy^2);
                 if du <= 0.001
                     u = u + 0.01;
                 end
                 u = u + du*T;
                 disp(u);
                 disp(V);
                 if u >= 1
                    u = 1;
                 end
                 param = [param,u];
               else
                  u = u + 0.01;
               end
            end
        end

        function u = Get_Interpolator_param(NURBS,SLine,T,u)
          V = SLine.Get_normlized_time_vel_value(u);
          if NURBS.Can_Derivatives(u) ~= false
              [dx,dy] = NURBS.Get_Derivatives_Control_points(NURBS.k,u,1);
              du = V/sqrt(dx^2+dy^2);
              if du <= 0.001
                     u = u + 0.01;
              end
              u = u + du*T;
          
              if u >= 1
                 u = 1;
              end
           else
               u = u + 0.01;
           end
        end

        function param = Create_Interpolator_Vlist(NURBS,T,V,u)
            %u(k+1) = u(k) + T*du(k)/dt
            %u(k+1) = u(k) + T*V(t(k))/sqrt((x')^2+(y')^2)
            param = [];
            param = [param,u];
            for i = 1:length(V)
                %Avoid the non-differentiable point
                if NURBS.Get_Derivatives_Control_points(NURBS.k,u,1) ~= false
                    [dx,dy] = NURBS.Get_Derivatives_Control_points(NURBS.k,u,1);
                     du = V(i)/sqrt(dx^2+dy^2);
                     if du <= 0.001
                        u = u + 0.01;
                     end
                     u = u + du*T;
                     if u >= 1
                         u = 1;
                     end
                     param = [param,u];
                else
                    u = u + 0.01;
                end
            end
        end

        function T = Create_Sample_param(NURBS,SampleSize)
            %need to use normalized knot to create 
            %will get SampleSize+1 data for 0 to 1
            T = 0:1/SampleSize:NURBS.knot(end);
        end

        function norm = normalized(NURBS)
            maximum = max(NURBS.knot);
            norm = NURBS.knot./ maximum;
        end

        function Same_knot = Get_Same_knot(NURBS)
            %find same parameter on knots
            %ex knot = [0,0,0,0.25,0.5,0.5,1,1,1]
            %return = [0,0.5,1]
            Same_knot = [];
            for i = 1:length(NURBS.knot)
                temp = NURBS.knot(i);
                count = 1;
                for j = i+1:length(NURBS.knot)
                    if NURBS.knot(j) == temp
                        count = count+1;
                    end
                end
                if NURBS.k - count <= 1
                    Same_knot = [Same_knot,temp];
                end
            end
            Same_knot = unique(Same_knot);
        end

        function Same_knot = Create_Same_knot(NURBS)
            %find same parameter on knots
            %ex knot = [0,0,0,0.25,0.5,0.5,1,1,1]
            %return = [0,0.5,1]
            Same_knot = [];
            for i = 1:length(NURBS.knot)
                count = 1;
                for j = 1:length(NURBS.knot)
                    if i == j
                        continue;
                    end
                    if NURBS.knot(i) == NURBS.knot(j)
                        count = count + 1;
                    end
                end
                if NURBS.k - count <= 1
                    Same_knot = [Same_knot,NURBS.knot(i)];
                end
            end
            Same_knot = unique(Same_knot);
        end

        function logical = Can_Derivatives(NURBS,u)
            % if now u parameter in Same_knot it can't be Derivative so will return false
            % else now u parameter not in Same_knot will return true
              if(find(u == NURBS.Same_knot))
                  logical = false;
              else
                  logical = true;
              end
        end
    end

    methods(Static)

    end
end

