function plot_mesh(msh,data,fignum)
figure(fignum); clf; hold on 
triplot(msh.t(1:3,:)',msh.p(1,:),msh.p(2,:)); % plotting the mesh 
axis equal
% filling nice colors for different materials

% windings 
mat_index_phaseA = [41 44];
mat_index_phaseB = [42 45];
mat_index_phaseC = [43 46];

for ele= 1: size(msh.t,2)
    if ismember(data.imat(ele),mat_index_phaseA)
        fill(msh.p(1,msh.t(1:3,ele)),msh.p(2,msh.t(1:3,ele)),'r')
    elseif ismember(data.imat(ele),mat_index_phaseB)
        fill(msh.p(1,msh.t(1:3,ele)),msh.p(2,msh.t(1:3,ele)),'y')
        elseif ismember(data.imat(ele),mat_index_phaseC)
        fill(msh.p(1,msh.t(1:3,ele)),msh.p(2,msh.t(1:3,ele)),'b')
    elseif data.imat(ele)>300 % permanent magnet 
        fill(msh.p(1,msh.t(1:3,ele)),msh.p(2,msh.t(1:3,ele)),'m')
    end 
      
end 

end 
