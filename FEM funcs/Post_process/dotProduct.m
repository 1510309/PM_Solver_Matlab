function E = dotProduct(V1, V2)
E = sum(V1.*V2, 1);
%E = sum(bsxfun(@times, V1, V2), 1);
end