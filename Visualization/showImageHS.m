function IMN=showImageHS(I_F,print,id,flag_cut_bounds,dim_cut,thvalues,L,bands)

if flag_cut_bounds
    I_F = I_F(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
end

if thvalues
    I_F(I_F > 2^L) = 2^L;
    I_F(I_F < 0) = 0;
end

IMN = viewimage(I_F(:,:,bands));

IMN = IMN(:,:,3:-1:1);

if print
    printImage(IMN,sprintf('../Outputs/%d.eps',id));
end

end