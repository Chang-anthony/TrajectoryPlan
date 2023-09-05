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

        %normlized 
        du 
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

               SLine.Jerk_max = Amax/SLine.Tc;
               SLine.t1 = SLine.Tc;
               SLine.t2 = SLine.Tc+SLine.Tb;
               SLine.t3 = SLine.Ta;
               SLine.t4 = SLine.Ta+SLine.Ts;
               SLine.t5 = SLine.Ta+SLine.Ts+SLine.Tc;
               SLine.t6 = SLine.Ta+SLine.Ts+SLine.Tc+SLine.Tb;
               SLine.t7 = 2*SLine.Ta+SLine.Ts;

               SLine.du = SLine.t7;
           else
               disp("Aavg is not match SLine,can't create");
           end
       end
       
       function [points_acc,points_t] = Create_Acc_points(SLine,step)
           points_t = [];
           points_acc = [];  
           for i = 1:step
              dt = 0+(SLine.t1)/step*i;
              points_acc = [points_acc,(SLine.Amax-0)/SLine.Tc*dt];
              points_t = [points_t,dt];
           end

           for i = 1:step
               dt = (SLine.t2-SLine.t1)/step*i;
               points_acc = [points_acc,(SLine.Amax+0*dt)];
               t = dt + SLine.t1;
               points_t = [points_t,t];
           end
            
           for i = 1:step
               dt = (SLine.t3-SLine.t2)/step*i;
               points_acc = [points_acc,(SLine.Amax+(0-SLine.Amax)/(SLine.Tc)*dt)];
               t = dt + SLine.t2;
               points_t = [points_t,t];
           end
           
           for i = 1:step
               dt = SLine.t3+(SLine.t4-SLine.t3)/step*i;
               points_acc = [points_acc,0];
               points_t = [points_t,dt];
           end
            
           for i = 1:step
               dt = (SLine.t5-SLine.t4)/step*i;
               points_acc = [points_acc,(0-SLine.Amax)/SLine.Tc*dt];
               t = dt + SLine.t4;
               points_t = [points_t,t];
           end 
          
           for i = 1:step
              dt = (SLine.t6-SLine.t5)/step*i;
              points_acc = [points_acc,(-SLine.Amax+0*dt)];
              t = dt + SLine.t5;
              points_t = [points_t,t];
           end
           
           for i = 1:step
               dt = (SLine.t7-SLine.t6)/step*i;
               points_acc = [points_acc,(-SLine.Amax+(0+SLine.Amax)/SLine.Tc*dt)];
               t = dt + SLine.t6;
               points_t = [points_t,t];
           end     
       end
       
       function points_vel = Create_Vel_points(SLine,step)
           points_vel = [];
           for i = 1:step
              dt = 0+(SLine.t1)/step*i;
              points_vel = [points_vel,(0.5*SLine.Amax-0)/SLine.Tc*dt^2];
           end
            
           % v = v0 + at;
           temp = points_vel(end);
           for i = 1:step
               dt = (SLine.t2-SLine.t1)/step*i;
               points_vel = [points_vel,temp+(SLine.Amax*dt+0*dt)];
           end
           
           temp = points_vel(end);
           for i = 1:step
               dt = (SLine.t3-SLine.t2)/step*i;
               points_vel = [points_vel,temp+SLine.Amax*dt+(0-0.5*SLine.Amax)/(SLine.Tc)*dt^2];
           end
           
           temp = points_vel(end);
           for i = 1:step
               dt = (SLine.t4-SLine.t3)/step*i;
               points_vel = [points_vel,temp+ 0*dt];
           end
           
           temp = points_vel(end);
           for i = 1:step
               dt = (SLine.t5-SLine.t4)/step*i;
               points_vel = [points_vel,temp+(0-0.5*SLine.Amax)/SLine.Tc*dt^2];
           end 
          
           temp = points_vel(end);
           for i = 1:step
              dt = (SLine.t6-SLine.t5)/step*i;
              points_vel = [points_vel,temp+(-SLine.Amax*dt+0*dt)];
           end
           
           temp = points_vel(end);
           for i = 1:step
               dt = (SLine.t7-SLine.t6)/step*i;
               points_vel = [points_vel,(temp-SLine.Amax*dt+(0+0.5*SLine.Amax)/SLine.Tc*dt^2)];
           end 
       end
       
       function points_pos = Create_Pos_points(SLine,step)
           points_pos = [];
           v = [];
           for i = 1:step
              dt = 0+(SLine.t1)/step*i;
              v = [v,(0.5*SLine.Amax-0)/SLine.Tc*dt^2];
              %pos = ((SLine.Amax-0)/(6*SLine.Tc))*dt^3;
              %pos = SLine.Amax/(6*SLine.Tc) * dt;
              points_pos = [points_pos,v(end)*dt];
           end
            
           %S = S0+0.5at^2
           temp = points_pos(end);
           tempv = v(end);
           for i = 1:step
               dt = (SLine.t2-SLine.t1)/step*i;
               v = [v,tempv+(SLine.Amax*dt+0*dt)];
               %pos = (0.5*SLine.Amax*dt^2+0*dt);
               %pos = SLine.Amax*SLine.Tc*SLine.Tc/(6) -0.5*SLine.Amax*SLine.Tc*dt +0.5*SLine.Amax*dt*dt;
               points_pos = [points_pos,temp+v(end)*dt];
           end
           
           temp = points_pos(end);
           tempv = v(end);
           for i = 1:step
               dt = (SLine.t3-SLine.t2)/step*i;
               v = [v,tempv+SLine.Amax*dt+(0-0.5*SLine.Amax)/(SLine.Tc)*dt^2];
               %pos = 0.5*SLine.Amax*dt^2+((0-SLine.Amax)*dt^3)/(SLine.Tc*6);
               %pos = ((0.5*SLine.Amax*SLine.Tc+SLine.Amax*SLine.Tb)*(dt-SLine.Tc-SLine.Tb)+0.5*SLine.Amax*(dt-SLine.Tc-SLine.Tb)^2 ...
               %-SLine.Amax*(dt-SLine.Tc-SLine.Tb)^3/(6*SLine.Tc)+SLine.Amax*SLine.Tc*SLine.Tc/6+0.5*SLine.Amax*SLine.Tb*(SLine.Tb+SLine.Tc));
               points_pos = [points_pos,temp+v(end)*dt];
           end
           
           temp = points_pos(end);
           tempv = v(end);
           for i = 1:step
               dt = (SLine.t4-SLine.t3)/step*i;
               v = [v,tempv+ 0*dt];
               points_pos = [points_pos,temp+ tempv*dt];
           end
           
           temp = points_pos(end);
           tempv = v(end);
           for i = 1:step
               dt =(SLine.t5-SLine.t4)/step*i;
               v = [v,tempv+(-0.5*SLine.Amax)/SLine.Tc*dt^2];
               %pos = ((SLine.Amax-0)*dt^3/(SLine.Tc*6));
               %pos = -SLine.Amax/(6*SLine.Tc)*dt;
               pos = temp+v(end)*dt;
               points_pos = [points_pos,temp+v(end)*dt];
           end 
          
           temp = points_pos(end);
           tempv = v(end);
           for i = 1:step
              dt = (SLine.t6-SLine.t5)/step*i;
              v = [v,tempv+(-SLine.Amax*dt+0*dt)];
              pos = temp+v(end)*dt;
              %pos = (0.5*SLine.Amax*dt^2+0*dt);
              %pos = (SLine.Amax*SLine.Tc*SLine.Tc/(6) -0.5*SLine.Amax*SLine.Tc*dt +0.5*SLine.Amax*dt*dt);
              points_pos = [points_pos,temp+v(end)*dt];
           end
           
           temp = points_pos(end);
           tempv = v(end);
           for i = 1:step
               dt = (SLine.t7-SLine.t6)/step*i;
               v = [v,(tempv-SLine.Amax*dt+(0+0.5*SLine.Amax)/SLine.Tc*dt^2)];
               pos = temp+v(end)*dt;
               %pos = ((SLine.Amax-0)/(SLine.Tc*6))*dt^3;
               %pos = ((0.5*SLine.Amax*SLine.Tc+SLine.Amax*SLine.Tb)*(dt-SLine.Tc-SLine.Tb)+0.5*SLine.Amax*(dt-SLine.Tc-SLine.Tb)^2 ...
               %-SLine.Amax*(dt-SLine.Tc-SLine.Tb)^3/(6*SLine.Tc)+SLine.Amax*SLine.Tc*SLine.Tc/6+0.5*SLine.Amax*SLine.Tb*(SLine.Tb+SLine.Tc));
               points_pos = [points_pos,temp+v(end)*dt];
           end 
       end
       
       
       function [points_vel,points_t] = Create_Vel_points_Time(SLine,T)
           points_vel = [];
           points_t = [];
           dt = 0;
           while dt < SLine.t7
               if dt <= SLine.t1
                  points_vel = [points_vel,(0.5*SLine.Amax-0)/SLine.Tc*dt^2];
               elseif dt > SLine.t1 && dt <= SLine.t2
                   tempt = dt;
                   points_vel = [points_vel,(-0.5*SLine.Amax*SLine.Tc+SLine.Amax*tempt)];
               
               elseif dt > SLine.t2 && dt <= SLine.t3
                   tempt = dt;
                   v = SLine.Amax*(SLine.Tb+SLine.Tc)+(0-0.5*SLine.Amax)/(SLine.Tc)*(SLine.Tb+2*SLine.Tc)^2 +(SLine.Amax)*(SLine.Tb+2*SLine.Tc)*tempt/SLine.Tc-0.5*SLine.Amax/SLine.Tc*tempt^2;
                   points_vel = [points_vel,v];
               
               elseif dt > SLine.t3 && dt <= SLine.t4
                   points_vel = [points_vel,SLine.Vmax];
               
               elseif dt > SLine.t4 && dt <= SLine.t5
                  v = SLine.Vmax-0.5*SLine.Amax/SLine.Tc*(dt-SLine.Ta-SLine.Ts)^2;
                  points_vel = [points_vel,v];
              

               elseif dt > SLine.t5 && dt <= SLine.t6
                 tempt = SLine.t6-dt+SLine.t5;
                  v = abs(0.5*SLine.Amax*SLine.Tc-SLine.Amax*SLine.Tc-SLine.Amax*(tempt-SLine.Ta-SLine.Ts-SLine.Tc));
                  points_vel = [points_vel,v];
               else
                   v = 0.5*SLine.Amax*SLine.Tc-SLine.Amax*(dt-SLine.Ta-SLine.Ts-SLine.Tc-SLine.Tb)+0.5*SLine.Amax/SLine.Tc*(dt-SLine.Ta-SLine.Ts-SLine.Tc-SLine.Tb)^2;
                   points_vel = [points_vel,v];
               end
               dt = dt + T;
               points_t = [points_t,dt];
           end
       end
        
       function V = Get_normlized_time_vel_value(SLine,u)
           %取得 0~1 之間，給定u值的時間並回傳當下的V值
           %according to chapter two  u = (t - t0)/du
           %so t is t = u*du + t0
           dt = u*SLine.du + 0;
           if dt <= SLine.t1
               V = (0.5*SLine.Amax-0)/SLine.Tc*dt^2;
           elseif dt > SLine.t1 && dt <= SLine.t2
               tempt = dt;
               V = (-0.5*SLine.Amax*SLine.Tc+SLine.Amax*tempt); 
           elseif dt > SLine.t2 && dt <= SLine.t3
               tempt = dt;
               V = SLine.Amax*(SLine.Tb+SLine.Tc)+(0-0.5*SLine.Amax)/(SLine.Tc)*(SLine.Tb+2*SLine.Tc)^2 +(SLine.Amax)*(SLine.Tb+2*SLine.Tc)*tempt/SLine.Tc-0.5*SLine.Amax/SLine.Tc*tempt^2;
           elseif dt > SLine.t3 && dt <= SLine.t4
               V = SLine.Vmax;
           elseif dt > SLine.t4 && dt <= SLine.t5
               V = SLine.Vmax-0.5*SLine.Amax/SLine.Tc*(dt-SLine.Ta-SLine.Ts)^2;
           elseif dt > SLine.t5 && dt <= SLine.t6
               tempt = SLine.t6-dt+SLine.t5;
               V = abs(0.5*SLine.Amax*SLine.Tc-SLine.Amax*SLine.Tc-SLine.Amax*(tempt-SLine.Ta-SLine.Ts-SLine.Tc));
           else
               V = 0.5*SLine.Amax*SLine.Tc-SLine.Amax*(dt-SLine.Ta-SLine.Ts-SLine.Tc-SLine.Tb)+0.5*SLine.Amax/SLine.Tc*(dt-SLine.Ta-SLine.Ts-SLine.Tc-SLine.Tb)^2;
           end
       end
   end
       
      
   methods(Static)
       
 
   end  
end