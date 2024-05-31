'''
    duplicate_begin_nodes_mask = np.isin(stacked_connections[0], unique_indices, invert=True)
    duplicate_end_nodes_mask = np.isin(stacked_connections[1], unique_indices, invert=True)

    # Find the coordinates of the duplicate begin nodes
    duplicate_begin_nodes_coordinates = rounded_stacked_coordinates[:, stacked_connections[0][duplicate_begin_nodes_mask]]
    duplicate_end_nodes_coordinates = rounded_stacked_coordinates[:, stacked_connections[1][duplicate_end_nodes_mask]]


    print(np.round(unique_nodes.T[0,:], decimals=2))
    print(duplicate_end_nodes_coordinates[0,:])
    indices = np.where(unique_nodes.T[0, :] == duplicate_begin_nodes_coordinates[0, :])[0]

    x_overlap_begin = np.where(np.isin(np.round(unique_nodes.T[0,:], decimals=2), duplicate_begin_nodes_coordinates[0,:]))[0]
    y_overlap_begin = np.where(np.isin(np.round(unique_nodes.T[1, :], decimals=2), duplicate_begin_nodes_coordinates[1, :]))[0]
    z_overlap_begin = np.where(np.isin(np.round(unique_nodes.T[2, :], decimals=2), duplicate_begin_nodes_coordinates[2, :]))[0]

    x_overlap_end = np.where(np.isin(np.round(unique_nodes.T[0, :], decimals=2), duplicate_end_nodes_coordinates[0, :]))[0]
    y_overlap_end = np.where(np.isin(np.round(unique_nodes.T[1, :], decimals=2), duplicate_end_nodes_coordinates[1, :]))[0]
    z_overlap_end = np.where(np.isin(np.round(unique_nodes.T[2, :], decimals=2), duplicate_end_nodes_coordinates[2, :]))[0]



    #a = np.round(unique_nodes.T[0,:], decimals=2)
    #b = duplicate_begin_nodes_coordinates[0,:]
    #print(a)
    #print(b)
    #print()
    #print(np.where(np.isin(a, b))[0])
    #print(np.where(np.all(a == b)))
    #print(x_overlap)
    #print(y_overlap)
    #print(z_overlap)
    print(x_overlap_begin)
    print(y_overlap_begin)
    print(z_overlap_begin)
    print(np.intersect1d(np.intersect1d(x_overlap_begin, y_overlap_begin), z_overlap_begin))

    unique_indices_to_replace_dupl_begins = np.intersect1d(np.intersect1d(x_overlap_begin, y_overlap_begin), z_overlap_begin)
    unique_indices_to_replace_dupl_ends = np.intersect1d(np.intersect1d(x_overlap_end, y_overlap_end), z_overlap_end)

    print(unique_indices_to_replace_dupl_begins.shape, stacked_connections[0][duplicate_begin_nodes_mask].shape)
    stacked_connections[0,:][duplicate_begin_nodes_mask] = unique_indices_to_replace_dupl_begins


    #print('a', duplicate_begin_nodes_mask)
    #print('b', stacked_connections[0][duplicate_begin_nodes_mask])
    #print('c', stacked_coordinates[:, stacked_connections[0][duplicate_begin_nodes_mask]])
    #print('d', np.where(stacked_coordinates == stacked_coordinates[:, stacked_connections[0][duplicate_begin_nodes_mask]]))
    #print(node_indices)
    print('----')

    #print(stacked_coordinates[:, duplicate_begin_nodes_mask])
    remapped_begin_nodes = np.where(stacked_coordinates[:, duplicate_begin_nodes_mask])


    remapped_begin_nodes = np.where(duplicate_begin_nodes_mask, stacked_connections[0], unique_indices[reverse_indices[stacked_connections[0]]])
    remapped_end_nodes = np.where(duplicate_end_nodes_mask, stacked_connections[1], unique_indices[reverse_indices[stacked_connections[1]]])

    # Stack remapped nodes to form remapped connections
    remapped_connections = np.vstack((remapped_begin_nodes, remapped_end_nodes))

    print(remapped_connections)
    print()
    print(self.n_per_hex)
    print(unique_edges.shape)
    '''


'''
        for idx, uc in enumerate(rounded_unique_nodes):
            x_equals = np.where(rounded_stacked_coordinates[0] == uc[0])
            y_equals = np.where(rounded_stacked_coordinates[1] == uc[1])
            z_equals = np.where(rounded_stacked_coordinates[2] == uc[2])
            common_indices = np.intersect1d(np.intersect1d(x_equals, y_equals), z_equals)
            #print(common_indices, common_indices[0])
            #print(np.where(stacked_connections[1] == common_indices[0]))
            print('====')
            stacked_connections[0][np.where(stacked_connections[0] == common_indices)] = idx
            stacked_connections[1][np.where(stacked_connections[1] == common_indices)] = idx


        #print(np.isin(stacked_connections, unique_indices))
        '''

global_coords = self.all_hex_XYZ[0]
global_indices = np.arange(0, self.n_per_hex)

for i in range(1, self.all_hex_connect.shape[0]):
    full_hex_coords = self.all_hex_XYZ[i]
    # print(full_hex_coords)
    # print(global_coords)
    full_hex_connections = self.all_hex_connect[i]
    full_hex_indices = np.arange(i * self.n_per_hex, (i + 1) * self.n_per_hex)
    boolean = np.isin(full_hex_coords, global_coords).all(0)

    b1 = np.isin(global_coords[0], full_hex_coords[0])
    b2 = np.isin(global_coords[0], full_hex_coords[0])
    b3 = np.isin(global_coords[0], full_hex_coords[0])
    b = np.logical_and(b3, np.logical_and(b1, b2))

    # print(global_coords)
    # print(full_hex_coords)
    # print(np.isin(global_coords[0], full_hex_coords[0]))
    # print(global_indices[np.isin(global_coords, full_hex_coords).all(0)])
    full_hex_indices[boolean] = global_indices[boolean]

    # print()
    # for j in range(full_hex_coords.shape[1]):

