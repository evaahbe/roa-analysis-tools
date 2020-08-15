function [h_vec_sdp] = make_hvec_sdpvar(x,sys)

    hx = sym('hx',[1 sys.xdim]); 
    h_vec_sym = reshape(sys.h_vec,size(sys.h_vec,1)*size(sys.h_vec,2),1);
    for i = 1:length(h_vec_sym)
        h_vec_char = char(h_vec_sym(i,1)); %doing this in a loop is not possible with this syntax
        if sys.xdim ==1
            h_vec_char = strrep(h_vec_char,char(hx(1)),'x(1)');
        elseif sys.xdim == 2
            h_vec_char = strrep(h_vec_char,char(hx(1)),'x(1)');
            h_vec_char = strrep(h_vec_char,char(hx(2)),'x(2)');
        elseif sys.xdim == 3
            h_vec_char = strrep(h_vec_char,char(hx(1)),'x(1)');
            h_vec_char = strrep(h_vec_char,char(hx(2)),'x(2)');
            h_vec_char = strrep(h_vec_char,char(hx(3)),'x(3)');
        elseif sys.xdim == 4
            h_vec_char = strrep(h_vec_char,char(hx(1)),'x(1)');
            h_vec_char = strrep(h_vec_char,char(hx(2)),'x(2)');
            h_vec_char = strrep(h_vec_char,char(hx(3)),'x(3)');
            h_vec_char = strrep(h_vec_char,char(hx(4)),'x(4)');
        elseif sys.xdim == 5
            h_vec_char = strrep(h_vec_char,char(hx(1)),'x(1)');
            h_vec_char = strrep(h_vec_char,char(hx(2)),'x(2)');
            h_vec_char = strrep(h_vec_char,char(hx(3)),'x(3)');
            h_vec_char = strrep(h_vec_char,char(hx(4)),'x(4)');
            h_vec_char = strrep(h_vec_char,char(hx(5)),'x(5)');
        elseif sys.xdim == 6
            h_vec_char = strrep(h_vec_char,char(hx(1)),'x(1)');
            h_vec_char = strrep(h_vec_char,char(hx(2)),'x(2)');
            h_vec_char = strrep(h_vec_char,char(hx(3)),'x(3)');
            h_vec_char = strrep(h_vec_char,char(hx(4)),'x(4)');
            h_vec_char = strrep(h_vec_char,char(hx(5)),'x(5)');
            h_vec_char = strrep(h_vec_char,char(hx(6)),'x(6)');
        else 
            error('System dimension not supported for FB design, add the according script lines')
        end
        h_vec_char = eval(h_vec_char);
        h_vec_sdp(i,1) = h_vec_char;
    end
end