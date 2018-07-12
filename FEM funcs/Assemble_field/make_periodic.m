function jjp= make_periodic(jj,a1,a2,a3,h)
% making the matrix suitable to use with periodic boundary condition s
%
% a2 and a3 are on periodic boundry and h is periodic coefficient
 jjp= [jj(a1,a1) jj(a1,a2)+h*jj(a1,a3);jj(a2,a1)+h*jj(a3,a1) jj(a2,a2)+jj(a3,a3)+h*jj(a2,a3)+h*jj(a3,a2)];
end 