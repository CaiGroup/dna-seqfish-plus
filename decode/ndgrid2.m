function Xall = ndgrid2(input)
%#codegen
nout = length(input);

    j = 1:nout;
    siz = zeros(1,nout);
    for i = 1:nout
        siz(i) = numel(input{i});
    end

Xall = input;
if nout == 2 % Optimized Case for 2 dimensions
    x = full(input{j(1)}(:));
    y = full(input{j(2)}(:)).';
    Xall{1} = repmat(x,size(y));
    Xall{2} = repmat(y,size(x));
else
    for i=1:nout
        x = full(input{j(i)});
        s = ones(1,20); 
        s(i) = numel(x);
        x1 = zeros(s);
        for k = 1:numel(x)
            x1(k) = x(k);
        end
        %x = reshape(x,s);
        s = ones(1,20); 
        for k = 1:numel(siz)
            s(k) = siz(k);
        end 
        s(i) = 1;
        Xall{i} = repmat(x1,s);
    end
end