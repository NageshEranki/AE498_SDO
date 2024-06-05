function eqn_num = list_dofs( e , Nx )
    
    eqn_num = zeros( length(e) , 8);

    nodes = list_nodes( e , Nx );

    for i = 1:size(nodes,2)

        eqn_num( : ,  2*i-1 ) = 2 * nodes ( : ,  i ) - 1;

        eqn_num( : ,  2*i )  = 2 * nodes( : ,  i );

    end

end