classdef Path <handle
    %% MEMBERS
      properties
        id
        lines
      end
   %% METHODS
   methods
       function path = Path()
           path.id = 0;
           path.lines = [];
       end
       
       function addline(path,line)
           path.id = length(path.lines)+1;
           path.lines = [path.lines,line];
       end
       
       function insert(path,line,id)
           path.id = length(path.lines)+1;
           path.lines = [path.lines(1:id-1) line path.lines(id:end)];
       end
       
       function remove(path,id)
           if path.id >= id
               path.lines(id) = [];
               path.id = length(path.lines);
           else
               disp("can not find this id line");
           end
       end
       
       function create_lines_table(path,table)
           line_amount = length(table);
           for i = 1:line_amount-1
               bound1 = table(i,:);%取出第i個行
               bound2 = table(i+1,:);%取出第i+1個行
               line = Line(bound1,bound2);
               path.addline(line);
           end
       end
       
       function create_lines_table_normalization(path,table)
           line_amount = length(table);
           for i = 1:line_amount-1
               bound1 = table(i,:);
               bound2 = table(i+1,:);
               line = Line(bound1,bound2);
               path.addline(line);
           end
       end
       
       function [xt,yt] = get_time_point_value(path,t,id)
           if t >= path.lines(id).t0  && t <= path.lines(id).tf
               [xt,yt] = path.lines(id).clc_value_use_time_normalization(t);
           else
               disp("can not match this line id");
               xt = 0;
               yt = 0;
           end
       end
       
        function [xt,yt,t] = get_time_point_value_normalization(path,u,id)
           %因為都已經確定被正規化了所以只要判斷這條線存不存在就好
           if path.id >= id
               [xt,yt,t] = path.lines(id).clc_value_normalization(u);
           else
               disp("can not match this line id");
               xt = 0;
               yt = 0;
               t = 0;
           end
       end
       
       function [all_pointx,all_pointy,all_pointt] = create_path_lines_point(path,step)
           all_pointx = [] ;
           all_pointy = [] ;
           all_pointt = [] ;
           for i = 1:length(path.lines)
               [tmpx,tmpy,tmpt] = path.lines(i).create_line_trajectory(step);
               all_pointx = [all_pointx,tmpx];
               all_pointy = [all_pointy,tmpy];
               all_pointt = [all_pointt,tmpt];
           end
       end
       
       function [all_pointx,all_pointy,all_pointt] = create_path_lines_point_normalization(path,step)
           all_pointx = [] ;
           all_pointy = [] ;
           all_pointt = [] ;
           for i = 1:length(path.lines)
               [tmpx,tmpy,tmpt] = path.lines(i).create_line_trajectory_normalization(step);
               all_pointx = [all_pointx,tmpx];
               all_pointy = [all_pointy,tmpy];
               all_pointt = [all_pointt,tmpt];
           end
       end
   end
   
    methods(Static)
        
    end
end
