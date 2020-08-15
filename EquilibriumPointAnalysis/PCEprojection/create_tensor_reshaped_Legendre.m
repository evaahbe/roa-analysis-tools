
function tensor = create_tensor_reshaped_Legendre(proj,order,dim)
%proj = basis onto which projected, order = number until which the sum
%goes, dim = amount of basis polynomials multiplied with each other before
%projected.


    if dim <1 || dim > 13 || mod(dim,1) ~=0
        error('dimension of tensor not supported')
    elseif dim == 1
        for j = 0:order
            tensor(j+1) = legendreProduct(proj,j)/legendreProduct(proj,proj);
        end
    elseif dim == 2
        for j = 0:order
            for k = 0:order
               tensor(j+1,k+1) = legendreProduct(proj,j,k)/legendreProduct(proj,proj);
            end
        end
    elseif dim == 3
        for j = 0:order
            for k = 0:order
                for l = 0:order
                    tensor(j+1,k+1,l+1) = legendreProduct(proj,j,k,l)/legendreProduct(proj,proj);
                end
            end
        end    
    elseif dim == 4
        for j = 0:order
            for k = 0:order
                for l = 0:order
                    for p = 0:order
                        tensor(j+1,k+1,l+1,p+1) = legendreProduct(proj,j,k,l,p)/legendreProduct(proj,proj);
                    end
                end
            end
        end   
    elseif dim == 5
        for j = 0:order
            for k = 0:order
                for l = 0:order
                    for p = 0:order
                        for r = 0:order
                            tensor(j+1,k+1,l+1,p+1,r+1) = legendreProduct(proj,j,k,l,p,r)/legendreProduct(proj,proj);
                        end
                    end
                end
            end
        end
    elseif dim == 6
        for j = 0:order
            for k = 0:order
                for l = 0:order
                    for p = 0:order
                        for r = 0:order
                            for q = 0:order
                                tensor(j+1,k+1,l+1,p+1,r+1,q+1) = legendreProduct(proj,j,k,l,p,r)/legendreProduct(proj,proj);
                            end
                        end
                    end
                end
            end
        end
    elseif dim == 7
        for j = 0:order
            for k = 0:order
                for l = 0:order
                    for p = 0:order
                        for r = 0:order
                            for q = 0:order
                                for s = 0:order
                                    tensor(j+1,k+1,l+1,p+1,r+1,q+1,s+1) = legendreProduct(proj,j,k,l,p,r,s)/legendreProduct(proj,proj);
                                end
                            end
                        end
                    end
                end
            end
        end
    elseif dim == 8
        for j = 0:order
            for k = 0:order
                for l = 0:order
                    for p = 0:order
                        for r = 0:order
                            for q = 0:order
                                for s = 0:order
                                    for t = 0:order
                                        tensor(j+1,k+1,l+1,p+1,r+1,q+1,s+1,t+1) = legendreProduct(proj,j,k,l,p,r,s,t)/legendreProduct(proj,proj);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end 
    elseif dim == 9
        for j = 0:order
            for k = 0:order
                for l = 0:order
                    for p = 0:order
                        for r = 0:order
                            for q = 0:order
                                for s = 0:order
                                    for t = 0:order
                                        for u = 0:order
                                            tensor(j+1,k+1,l+1,p+1,r+1,q+1,s+1,t+1,u+1) = legendreProduct(proj,j,k,l,p,r,s,t,u)/legendreProduct(proj,proj);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end 
    elseif dim == 10
        for j = 0:order
            for k = 0:order
                for l = 0:order
                    for p = 0:order
                        for r = 0:order
                            for q = 0:order
                                for s = 0:order
                                    for t = 0:order
                                        for u = 0:order
                                            for v = 0:order
                                                tensor(j+1,k+1,l+1,p+1,r+1,q+1,s+1,t+1,u+1,v+1) = legendreProduct(proj,j,k,l,p,r,s,t,u,v)/legendreProduct(proj,proj);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end 
    elseif dim == 11
        for j = 0:order
            for k = 0:order
                for l = 0:order
                    for p = 0:order
                        for r = 0:order
                            for q = 0:order
                                for s = 0:order
                                    for t = 0:order
                                        for u = 0:order
                                            for v = 0:order
                                                for w = 0:order
                                                    tensor(j+1,k+1,l+1,p+1,r+1,q+1,s+1,t+1,u+1,v+1,w+1) = legendreProduct(proj,j,k,l,p,r,s,t,u,v,w)/legendreProduct(proj,proj);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end 
    elseif dim == 12
        for j = 0:order
            for k = 0:order
                for l = 0:order
                    for p = 0:order
                        for r = 0:order
                            for q = 0:order
                                for s = 0:order
                                    for t = 0:order
                                        for u = 0:order
                                            for v = 0:order
                                                for w = 0:order
                                                    for y = 0:order
                                                        tensor(j+1,k+1,l+1,p+1,r+1,q+1,s+1,t+1,u+1,v+1,w+1,y+1) = legendreProduct(proj,j,k,l,p,r,s,t,u,v,w,y)/legendreProduct(proj,proj);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end 
    elseif dim == 13
        for j = 0:order
            for k = 0:order
                for l = 0:order
                    for p = 0:order
                        for r = 0:order
                            for q = 0:order
                                for s = 0:order
                                    for t = 0:order
                                        for u = 0:order
                                            for v = 0:order
                                                for w = 0:order
                                                    for y = 0:order
                                                        for z = 0:order
                                                            tensor(j+1,k+1,l+1,p+1,r+1,q+1,s+1,t+1,u+1,v+1,w+1,y+1,z+1) = legendreProduct(proj,j,k,l,p,r,s,t,u,v,w,y,z)/legendreProduct(proj,proj);
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end 
    end    
    
    if dim ==1
        return
    else
        tensor = clean(reshape(tensor,(order+1)^(dim-1),(order+1)),1e-6);
    end
end


