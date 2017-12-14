module myFun
    
    export smooth_tract_selective_PCA

    function smooth_tract_selective_PCA(tmparray::Array{Float32},tmparray2Dr::Array{Float32},magicnumber::Int64)
        tmparraydiff = abs.(tmparray2Dr.-tmparray)
        tmparraydiff = sum(tmparraydiff,1)                  
        tmpcorrI = sortperm(tmparraydiff[:])
        tmparray_chosen = tmparray2Dr[:,tmpcorrI[1:magicnumber]]
        mean_tmparray_chosen = mean(tmparray_chosen,1);
        tmparray_chosen = tmparray_chosen .- mean_tmparray_chosen;
        u,s,v = svd(tmparray_chosen);
        s[2]=s[2]/2f0;
        s[3:end]=0;
        tmparray_chosen2 = u*Diagonal(s)*v';
        tmparray_chosen2 = tmparray_chosen2 .+ mean_tmparray_chosen;
        seed_smoothed_pca = tmparray_chosen2[:,1];
        return seed_smoothed_pca
    end

end
