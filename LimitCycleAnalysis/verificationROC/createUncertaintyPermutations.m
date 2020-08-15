function [permutvector,uncercount] = createUncertaintyPermutations(uncer_dim)

    if uncer_dim == 0
        permutvector = 1;
        uncercount = 1;
    else
    
        %outputs a n_permutations x uncer_dim sized matrix of -1 and 1 entries.
        %Each row is a unique permutation of -1 and 1.
        permarray = repmat([-1,1],1,uncer_dim);
        allperms = nchoosek(permarray,uncer_dim);
        permutvector = unique(allperms, 'rows');

        uncercount = size(permutvector,1);
    end
end