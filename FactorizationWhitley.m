
matrix10 = [1 1 1 0 0 0 0 0 1 1;
            1 1 1 1 0 0 0 0 0 1;
            1 1 1 1 1 0 0 0 0 0; 
            0 1 1 1 1 1 0 0 0 0;
            0 0 1 1 1 1 1 0 0 0;
            0 0 0 1 1 1 1 1 0 0; 
            0 0 0 0 1 1 1 1 1 0; 
            0 0 0 0 0 1 1 1 1 1;
            1 0 0 0 0 0 1 1 1 1; 
            1 1 0 0 0 0 0 1 1 1];
        
        order = [1:10];
        [G, cliques, fill_ins] = triangulate(matrix10, order)
        [jtree, root, B, w] = cliques_to_jtree(cliques, ns)
        
        
        
