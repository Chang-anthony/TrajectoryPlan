classdef Line <handle
     %% MEMBERS
      properties
        bound1  %Boundary conditions1 
        bound2  %Boundary conditions2
        
        t0 %start time 
        tf %final time
        
        du %normal du = tf-t0
        
        R % out mat
        B % cubic function mat
        A % coefficient mat
        
        R_n % out mat normal
        B_n % cubic function mat normal
        A_n % coefficient mat normal
      end
   %% METHODS
   methods
       %%constructor
       function line = Line(bound1,bound2)
           line.bound1 =  bound1;
           line.bound2 =  bound2;
           
           line.t0 = bound1(1);
           line.tf = bound2(1);
           line.du = line.tf - line.t0;
           
           line.R = line.Create_Rmat(bound1,bound2);
           line.B = line.Create_Bmat(line.t0,line.tf);
           line.A = line.Clc_Amat(line.R,line.B);
           
           line.R_n = line.Create_Rmat_normalization(line,bound1,bound2);
           line.B_n = line.Create_Bmat_constant();
           line.A_n = line.Clc_Amat_normalization(line);
       end

      %計算送入該點的t值出來的值會多少
      function [x,y] = clc_value(line,t)
          now_cubic = Line.clc_cubic(t);
          out = line.A * now_cubic;
          x = out(1);
          y = out(2);
      end
      
      %計算送入該點的t值出來的值會多少,送入過後會計算正規畫區間
      function [x,y] = clc_value_use_time_normalization(line,t)
          u = (t-line.t0)/line.du;
          now_cubic = Line.clc_cubic(u);
          out = line.A_n * now_cubic;
          x = out(1);
          y = out(2);
      end
      
      function t =  clc_time_normalization(line,u)
          t = line.t0 + (line.tf-line.t0)*u;
      end
      
      function [x,y,t] = clc_value_normalization(line,u)
            now_cubic = Line.clc_cubic(u);
            out = line.A_n * now_cubic;
            t = line.clc_time_normalization(u);
            x = out(1);
            y = out(2);
      end
      
      %step 總分割點數 now t = t0 + (tf - t0 )/step * i
      function [pointsx,pointsy,pointst] = create_line_trajectory(line,step)
            pointsx = [];
            pointsy = [];
            pointst = [];
            
            %因為分割出來後會為 step+1 才會到終點
            for i = 1: step 
                t = line.t0 +(line.tf-line.t0)/step*i;
                now_cubic = Line.clc_cubic(t);
                out = line.A * now_cubic;
                pointsx = [pointsx,out(1)];
                pointsy = [pointsy,out(2)];
                pointst = [pointst,t];
            end  
      end
      
      %step 總分割點數  因normal完 範圍介於 0~1之間
      function [pointsx,pointsy,pointst] = create_line_trajectory_normalization(line,step)
            pointsx = [];
            pointsy = [];
            pointst = [];
            dt = 1/step;
            
            %因為分割出來後會為 step+1 才會到終點
            for i = 1: step
                u = 0 + dt*i;
                t = line.clc_time_normalization(u);
                now_cubic = Line.clc_cubic(u);
                out = line.A_n * now_cubic;
                pointsx = [pointsx,out(1)];
                pointsy = [pointsy,out(2)];
                pointst = [pointst,t];
            end  
      end
   end
   
   
   methods(Static)
       
     function line = Line_bound(t0,tf,xi_1,dxi_1,yi_1,dyi_1,xi,dxi,yi,dyi)
           bound1 = [t0,xi_1,dxi_1,yi_1,dyi_1];
           bound2 = [tf,xi,dxi,yi,dyi];
           
           line = Line(bound1,bound2);
     end
      
     function B = Create_Bmat_constant()
            B = [[0,0,1,3];
                 [0,0,1,2];
                 [0,1,1,1];
                 [1,0,1,0]];
      end
      
     function R = Create_Rmat(bound1,bound2)
           R = [[bound1(2),bound1(3),bound2(2),bound2(3)];
                [bound1(4),bound1(5),bound2(4),bound2(5)]];
       end
       
     function R = Create_Rmat_normalization(obj,bound1,bound2)
           R = [[bound1(2),bound1(3)*obj.du,bound2(2),bound2(3)*obj.du];
                [bound1(4),bound1(5)*obj.du,bound2(4),bound2(5)*obj.du]];
       end
       
     function B = Create_Bmat(t0,tf)
           
           B = [[t0*t0*t0,t0*t0*3,tf*tf*tf,tf*tf*3];
                [t0*t0,t0*2,tf*tf,2*tf];
                [t0,1,tf,1];
                [1,0,1,0]];
       end
       
     function A = Clc_Amat(R_mat,B_mat)
           inv_B = inv(B_mat);
           A = R_mat * inv_B ;
     end
       
     function A = Clc_Amat_normalization(line)
           inv_B = inv(line.B_n);
           A = line.R_n * inv_B ;
      end
      
     function n_cubic = clc_cubic(t)
                n_cubic = [t*t*t;
                            t*t;
                            t;
                            1];
     end
   end  
end