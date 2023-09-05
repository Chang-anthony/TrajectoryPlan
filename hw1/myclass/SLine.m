classdef SLine <handle
     %% MEMBERS
      properties
        Smax %Max position  
        Vmax %Max velocity 
        Amax %Max acceleration 
        Aavg %average acceleration 
        Jerk_max %Mat Jerk
        
        Ta
        Tb
        Tc
        Ts
        
        t1
        t2
        t3
        t4
        t5
        t6
        t7
      end
   %% METHODS
   methods
       %%constructor
       function SLine = SLine(Smax,Vmax,Amax,Aavg)
           if Aavg >= Amax/2 && Aavg < Amax
               SLine.Smax = Smax;
               SLine.Vmax = Vmax;
               SLine.Amax = Amax;
               SLine.Aavg = Aavg;
               SLine.Ta = Vmax/Aavg;
               SLine.Tb = 2*Vmax/Amax-SLine.Ta;
               SLine.Tc = (SLine.Ta-SLine.Tb)/2;
               SLine.Ts = (Smax-Vmax*SLine.Ta)/Vmax;

               SLine.t1 = SLine.Tc;
               SLine.t2 = SLine.Tc+SLine.Tb;
               SLine.t3 = SLine.Ta;
               SLine.t4 = SLine.Ta+SLine.Ts;
               SLine.t5 = SLine.Ta+SLine.Ts+SLine.Tc;
               SLine.t6 = SLine.Ta+SLine.Ts+SLine.Tc+SLine.Tb;
               SLine.t7 = 2*SLine.Ta+SLine.Ts;
           else
               disp("Aavg is not match SLine,can't create");
           end
       end
       
       function [points_acc,points_t] = clc_Acc_value(step)
           points_t = [];
           points_acc = [];
           dt = 0;
          if  dt < SLine.t1
              for i = 1:step
                  dt = dt + SLine.t1/step*i;
                  points_acc = [points_acc,(SLine.Amax-0)/SLine.Tc*dt];
                  points_t = [points_t,dt];
              end
          elseif dt >= SLine.t1 && dt < SLine.t2
             for i = 1:step
                  dt = dt + (SLine.t2-SLine.t1)/step*i;
                  points_acc = [points_acc,(SLine.Amax+0*dt)];
                  points_t = [points_t,dt];
             end
         elseif dt >= SLine.t2 && dt < SLine.t3
             for i = 1:step
                  dt = dt + (SLine.t3-SLine.t2)/step*i;
                  points_acc = [points_acc,(SLine.Amax+(0-SLine.Amax)/SLine.Tc*dt)];
                  points_t = [points_t,dt];
             end
         elseif dt >= SLine.t3 && dt < SLine.t4
             for i = 1:step
                  dt = dt + (SLine.t4-SLine.t3)/step*i;
                  points_acc = [points_acc,0];
                  points_t = [points_t,dt];
             end
         elseif dt >= SLine.t4 && dt < SLine.t5
             for i = 1:step
                  dt = dt + (SLine.t5-SLine.t4)/step*i;
                  points_acc = [points_acc,(0-SLine.Amax)/SLine.Tc*dt];
                  points_t = [points_t,dt];
             end 
         elseif dt >= SLine.t5 && dt < SLine.t6
             for i = 1:step
                  dt = dt + (SLine.t6-SLine.t5)/step*i;
                  points_acc = [points_acc,(-SLine.Amax+0*dt)];
                  points_t = [points_t,dt];
             end
         elseif dt >= SLine.t6 && dt < SLine.t7
             for i = 1:step
                  dt = dt + (SLine.t7-SLine.t6)/step*i;
                  points_acc = [points_acc,(-SLine.Amax+(0-SLine.Amax)/SLine.Tc*dt)];
                  points_t = [points_t,dt];
             end   
          end    
       end
       
      %�p��e�J���I��t�ȥX�Ӫ��ȷ|�h��
      function [x,y] = clc_value(line,t)
          now_cubic = Line.clc_cubic(t);
          out = line.A * now_cubic;
          x = out(1);
          y = out(2);
      end
      
      %�p��e�J���I��t�ȥX�Ӫ��ȷ|�h��,�e�J�L��|�p�⥿�W�e�϶�
      function [x,y] = clc_value_use_time_normalization(line,t)
          u = (t-line.t0)/line.du;
          now_cubic = Line.clc_cubic(u);
          out = line.A_n * now_cubic;
          x = out(1);
          y = out(2);
      end
         
      %step �`�����I�� now t = t0 + (tf - t0 )/step * i
      function [pointsx,pointsy,pointst] = create_line_trajectory(line,step)
            pointsx = [];
            pointsy = [];
            pointst = [];
            
            %�]�����ΥX�ӫ�|�� step+1 �~�|����I
            for i = 1: step 
                t = line.t0 +(line.tf-line.t0)/step*i;
                now_cubic = Line.clc_cubic(t);
                out = line.A * now_cubic;
                pointsx = [pointsx,out(1)];
                pointsy = [pointsy,out(2)];
                pointst = [pointst,t];
            end  
      end
      
      %step �`�����I��  �]normal�� �d�򤶩� 0~1����
      function [pointsx,pointsy,pointst] = create_line_trajectory_normalization(line,step)
            pointsx = [];
            pointsy = [];
            pointst = [];
            dt = 1/step;
            
            %�]�����ΥX�ӫ�|�� step+1 �~�|����I
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
       
 
   end  
end