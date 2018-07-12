function detF = mappingDeterminant(F)
detF = F(1,:).*F(4,:) - F(2,:).*F(3,:);
end