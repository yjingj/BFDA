function [x, t_x, y, newx, new_tx, newy, invalid] = pre1(x, t_x, y, isNewSub)

   invalid=0;
   newx = []; new_tx= []; newy = 0; 
   if length(x)~=length(t_x)        % check whether x and t_x have the same length
       fprintf(1,'Error: invalid input, x and t_x should have equal length!\n');
       invalid=1;
       return;
   end       
   
   if isempty(isNewSub) || all(isNewSub == 0)
      if length(y)~=length(x)       % check whether y and x have the same length
          fprintf(1,'Error: invalid input, x and y should have equal length!\n');
          invalid=1;
          return;
      end  
      newx = []; new_tx = []; newy = []; 
          
   elseif length(isNewSub) == 1 && isNewSub > 0            %when isNewSub is scalar, the last number 
      if isNewSub == length(x)
         fprintf(1,'Error: invalid isNewSub, at least one or more subjects are not for predicton!\n');
         invalid=1;
         return;
      end
      if length(y)>length(x)-isNewSub             % check length(y)
          fprintf(1,'Warning: length(y) should be length(x)-isNewSub. Reset it to be so now!\n');
          y = y(1:(length(x)-isNewSub));
      elseif length(y)<length(x)-isNewSub
          fprintf(1,'Error: invalid input of y. Please check length(y)!\n');
          invalid=1;
          return;
      end
      newx = x((end-isNewSub+1):end);                       %of subjects ("isNewSub") are used for prediction
	  new_tx = t_x((end-isNewSub+1):end);
      x = x(1:(end-isNewSub));
      t_x = t_x(1:(end-isNewSub));              
      
 
   elseif length(isNewSub) == length(x)             %when isNewSub is an indicator,1 : new subject , 0: observed subject
     if all(isNewSub == 1)
        fprintf(1,'Error: invalid isNewSub, at least one or more subjects are not for predicton!\n');
        invalid=1;
        return;
     end
     if length(y)==length(x)
         fprintf(1,'Warning: length(y) should be length(x)-sum(isNewSub). Reset it to be so now!\n');
         y = y(isNewSub==0);
     elseif length(y)~=length(x)-sum(isNewSub)
         fprintf(1,'Error: invalid input of y. Please check length(y)!\n');
         invalid=1;
         return;
     end
     newx = x(isNewSub == 1);
     new_tx = t_x(isNewSub == 1);
      
      x = x(isNewSub == 0);
      t_x = t_x(isNewSub == 0);
      
 
   else
      fprintf(1,'Warning: isNewSub must be a positive integer or a 0-1 indicator vector!Reset isNewSub = [] now!\n');
      newx = []; new_tx = []; newy = []; 
   end