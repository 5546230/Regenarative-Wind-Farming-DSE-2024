import numpy as np
import matplotlib.pyplot as plt


XYZ_coords = np.array([[-6, 12, 6, -12, 0],
                       [0, 0, 0, 0, 24],
                       [8, 8, -8, -8, 0]])

member_indices = np.array([[0, 1, 2, 3],
                           [4, 4, 4, 4]])
sections = np.array([0,0,0,0])
materials = np.array([0,0,0,0])


class Mesh:
    def __init__(self, XYZ_coords: np.array, member_indices: np.array, sections: np.array, materials: np.array):
        self.XYZ_coords = XYZ_coords
        self.X_coords = XYZ_coords[0,:]
        self.Y_coords = XYZ_coords[1, :]
        self.Z_coords = XYZ_coords[2, :]

        self.element_indices = member_indices
        self.section_indices = sections
        self.material_indices = materials
        self.element_lengths = self.calc_element_lengths()
        self.element_Es = np.ones_like(materials)
        self.element_As = np.ones_like(materials)

        self.transfer_matrix()

    def calc_element_lengths(self):
        start_points = self.XYZ_coords[:, member_indices[0, :]]
        end_points = self.XYZ_coords[:, member_indices[1, :]]

        differences = end_points - start_points
        lengths = np.linalg.norm(differences, axis=0)
        return lengths

    def transfer_matrix(self):
        I2 = np.eye(2)
        Ts = []
        Ls = self.element_lengths
        start_points = self.XYZ_coords[:, member_indices[0, :]]
        end_points = self.XYZ_coords[:, member_indices[1, :]]
        differences = end_points - start_points

        cosines = differences / Ls
        for i in range(cosines.shape[1]):
            semi_row = cosines[:,i]
            Ti = np.kron(I2, semi_row)
            Ts.append(Ti)
        Ts = np.array(Ts)
        return Ts

    def element_stiffness(self):
        k = np.array([[1, -1], [-1, 1]]) * (E * A) / L
        return k

    def plot_structure(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        for i in range(self.element_indices.shape[1]):
            member_ends =  self.element_indices[:,i]
            Xs = self.X_coords[member_ends]
            Ys = self.Y_coords[member_ends]
            Zs = self.Z_coords[member_ends]
            plt.plot(Xs, Ys, Zs, color='k')
        plt.show()


TEST = Mesh(XYZ_coords=XYZ_coords, member_indices=member_indices, sections=sections, materials=materials)
TEST.plot_structure()
print(TEST.element_lengths)
