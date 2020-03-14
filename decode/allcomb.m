function B = allcomb(args)
%#codegen

NC = length(args);

input = args;
for i = 1:length(input)
    input{i} = args{length(input)+1-i};
end

[XA] = ndgrid2(input) ;

A = input;
for i = NC:-1:1
    A{NC+1-i} = XA{i};
end

% concatenate
B = zeros(numel(A{1}),NC);
for i = 1:NC
    B(:,i) =reshape(A{i},[],1);
end

