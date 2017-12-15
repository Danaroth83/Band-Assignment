function I_out=changem1(I_in,value_new,value_old)

I_out=zeros(size(I_in));
for ii=1:length(value_old)
    I_out(I_in==value_old(ii))=value_new(ii);
end

end