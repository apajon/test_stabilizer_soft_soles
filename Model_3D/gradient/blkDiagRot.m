function RBlock = blkDiagRot(angleTraj)
RBlock = [];
for i=1:length(angleTraj)
    R = [cos(angleTraj(i)),-sin(angleTraj(i));
        sin(angleTraj(i)),cos(angleTraj(i))];
    RBlock = blkdiag(RBlock,R);
end
end