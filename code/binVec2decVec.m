function decVec = binVec2decVec(binVec, bits)
% convert binary vector BINVEC with bitlength BITS to normalized (-1, 1) decimal vector DECVEC

[popSize, nBases] = size(binVec);
nGenes = nBases/bits;
decVec = zeros(popSize, nGenes);
base = repmat(2.^(0:bits-1),popSize,nGenes);

tmp = binVec.*base;
for ind = 1:popSize
   decVec(ind,:) = sum(reshape(tmp(ind,:),[],nGenes));
end
decVec = (decVec - 2^(bits-1))/2^(bits-1);
