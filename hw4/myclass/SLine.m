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

               SLine.Jerk_max = Amax/SLine.Tc;
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
           for i = 1:step
              dt = 0+(SLine.t1)/step*i;
              pos = ((SLine.Amax-0)/(6*SLine.Tc))*dt^3;
              points_pos = [points_pos,pos];
           end
            
           %S = S0+0.5at^2
           temp = points_pos(end);
           for i = 1:step
               dt = (SLine.t2-0)/step*i;
               pos = (SLine.Amax*SLine.Tc^2)/(6) -0.5*SLine.Amax*SLine.Tc*dt +0.5*SLine.Amax*dt^2;
               points_pos = [points_pos,pos];
           end
           
           temp = points_pos(end);
           for i = 1:step
               dt = SLine.t2+SLine.Tc/step*i;
               pos = ((0.5*(SLine.Amax)*SLine.Tc+SLine.Amax*SLine.Tb)*(dt-SLine.Tc-SLine.Tb)+0.5*SLine.Amax*(dt-SLine.Tc-SLine.Tb)^2 ...
                -SLine.Amax*(dt-SLine.Tc-SLine.Tb)^3/(6*SLine.Tc)+SLine.Amax*SLine.Tc^2/6+0.5*SLine.Amax*SLine.Tb*(SLine.Tb+SLine.Tc));
               points_pos = [points_pos,pos];
           end
           
           temp = points_pos(end);
           for i = 1:step
               dt = (SLine.t4-SLine.t3)/step*i;
               points_pos = [points_pos,temp+SLine.Vmax*dt];
           end
           
           temp = points_pos(end);
           for i = 1:step
               dt = SLine.t4+SLine.Tc/step*i;
               pos = SLine.Vmax*(dt-SLine.t4) -SLine.Amax*(dt-SLine.t4)^3/(6*SLine.Tc);
               points_pos = [points_pos,temp+pos];
           end 
          
           temp = points_pos(end);
           for i = 1:step
              dt = (SLine.t5)+SLine.Tb/step*i;
              pos = (0.5*SLine.Vmax*(2*SLine.Ts+SLine.Ta)+(5/6)*SLine.Amax*SLine.Tc^2+SLine.Amax*SLine.Tc*SLine.Tb)...
                  +(0.5*SLine.Amax*SLine.Tc+SLine.Amax*SLine.Tb)*(dt-SLine.Ta-SLine.Ts-SLine.Tc)-(SLine.Amax*(dt-SLine.Ta-SLine.Ts-SLine.Tc)^2/2);
              points_pos = [points_pos,pos];
           end

           temp = points_pos(end);
           for i = 1:step
               dt = (SLine.t6)+SLine.Tc/step*i; 
               pos = (SLine.Smax-(SLine.Amax*SLine.Tc^2)/6)+0.5*SLine.Amax*SLine.Tc*(dt-SLine.Ta-SLine.Ts-SLine.Tc-SLine.Tb)...
                   -(SLine.Amax*((dt-SLine.Ta-SLine.Ts-SLine.Tc-SLine.Tb)^2)/2)+SLine.Amax*(dt-SLine.Ta-SLine.Ts-SLine.Tc-SLine.Tb)^3/(6*SLine.Tc);
               points_pos = [points_pos,pos];
           end
       end 
   end
       
      
   methods(Static)
       
 
   end  
end