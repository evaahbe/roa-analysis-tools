
%_____________________________________________________________________
%
%   BUILD FILE 
%
%   This script uses trajectory data of full periodic orbit and generates 
%   a ppval-type spline function.
%
%   Note 1:Even though spline(...) should do a similar job this was hand
%   coded for a higher precision (derivative is taken from dynamics, not
%   from data)
%
%   Note 2: The data has to have the same start and end location point
%   
%   Note 3: For later usage the spline function is generated over several
%   periods ('rounds') of the orbit (there might be an option in the mkpp
%   to make the spline periodic replacing this part of the script)
%
%   Written by
%   Eva Ahbe, 15-05-2020
%_____________________________________________________________________


function [xspline, uinter, udotinter, tvecout] = construct_spline_control(sys,fns,tvec0,xvec,uvec,udotvec) 
    
    rounds = sys.rounds;

    
    tvec = tvec0;
    
    % do this for several 'rounds' around the limit cycle
    for i = 2:1:rounds 
        tvec = [tvec;tvec0(2:end)+(i-1)*tvec0(end)]; 
    end
    
    tvecout = tvec; % the output time vector
    tvec = [tvec;tvec0(2)+(i)*tvec0(end)]; % for the derivative at the end point..
    xvec = [xvec;repmat(xvec(2:end,:),rounds,1)];
    uvec = [uvec;repmat(uvec(2:end,:),rounds-1,1)];
    uvec = [uvec;uvec(2,:)];
    udotvec = [udotvec;repmat(udotvec(2:end,:),rounds-1,1)];
    udotvec = [udotvec;udotvec(2,:)];

    
    % get derivatives for each element in xvec
    for i=1:length(tvec)
        xdotvec(i,:) = fns.system_dyn(xvec(i,:)',uvec(i,:),sys);
    end
    
    k = 0;
    for i = 1:length(tvec)-1
        dt = tvec(i+1)-tvec(i);
        A = [dt^3, dt^2; ...
            3*dt^2, 2*dt];
        for j = 1:sys.xdim
            B = [xvec(i+1,j) - xdotvec(i,j)*dt - xvec(i,j);...
                xdotvec(i+1,j) - xdotvec(i,j)];
            k = k+1;
            coefsol(k,:) = [linsolve(A,B)',xdotvec(i,j),xvec(i,j)];
        end
    end
    
    xspline = mkpp(tvec,coefsol,sys.xdim);

%        if strcmp(sys.uoption,'spline') 
    uinter = spline(tvec,uvec);
    udotinter = spline(tvec,udotvec); 
%         elseif 
%             strcmp(obj.uoption,'zoh') 
%             uinter = mkpp(tvec,uvec(1:end-1));
%         elseif strcmp(obj.uoption,'pchip')
%             uinter = pchip(tvec,uvec);
%        else
%            error('Interpolation method for input not known!')
%        end

    
end
