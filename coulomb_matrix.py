import numpy as np 
from atomic_number import atomic_number

class Molecule:
    """Represents a molecule and generates its Coulomb matrix based on atomic numbers and interatomic distances."""
    def __init__(self, file_name):
        """Initializes a molecule from an .xyz file"""
         
        with open(file_name, 'r') as file:
            self.natoms = int(file.readline().strip()) 
            file.readline()

            atom_list = []
            for row in range(self.natoms):
                # Split the file to index elements and coordinates seperately
                line = file.readline()
                parts = line.split()
                
                element_symbol = parts[0]
                element_symbol = element_symbol.strip().capitalize()
                coords = np.array([float(parts[1]), float(parts[2]), float(parts[3])], dtype=float)
                atom_dict = {
                    'element': element_symbol,
                    'coordinates': coords
                }
                atom_list.append(atom_dict)
            self.xyz = tuple(atom_list)
    


    def rawCM(self):
        """Returns an array of the Coulomb Matrix of the molecule"""
        c_matrix = np.zeros((self.natoms, self.natoms))

        for i in range(self.natoms):
            elem_i = self.xyz[i]['element']
            Zi = atomic_number.get(elem_i)
            Ri = self.xyz[i]['coordinates']
        
            for j in range(self.natoms):
                elem_j = self.xyz[j]['element']
                Zj = atomic_number.get(elem_j)
                Rj = self.xyz[j]['coordinates']
            
                if i == j: 
                    c_matrix[i][j] = 0.5 * (Zi ** 2.4)      # Formula for when atom i = atom j
                else:
                    distance = np.linalg.norm(Ri - Rj) 
                    c_matrix[i][j] = (Zi * Zj) / distance       # Formula for when atom i != atom j
  
        return c_matrix

    def eigenCM(self):
        """Returns an array with the eigenvalues of the rawCM, in descending order."""
        eigenval = np.sort(np.linalg.eigh(self.rawCM())[0])[::-1]
        return eigenval

    def eigenCM_distance(self, other):
        """Computes the Euclidean distance between the eigenvalues of this molecule and another molecule"""
        # Create if loops for the different sized matrices
        if self.natoms != other.natoms:
            if self.natoms < other.natoms:
                smaller = self
                larger = other
            else:
                smaller = other
                larger = self

            pad_width = larger.natoms - smaller.natoms
            self_eigenval = np.pad(smaller.eigenCM(), (0, pad_width), 'constant')
            other_eigenval = larger.eigenCM()

        else:
            self_eigenval = self.eigenCM()
            other_eigenval = other.eigenCM()

        distance = np.linalg.norm(self_eigenval - other_eigenval)
        return distance

    def sortedCM(self):
        """Returns an array with the sorted Coulomb matrix of the molecule"""
        c_matrix = self.rawCM()
       
        row_norms = []

        for row in c_matrix:
            norm = np.linalg.norm(row)
            row_norms.append(norm)

        row_norms = np.array(row_norms)

        sorted_indices = np.argsort(row_norms)[::-1]

        sorted_cm = c_matrix[sorted_indices][:, sorted_indices]

        return sorted_cm

    def sortedCM_distance(self, other):
        """Computes the Frobenius Norm between the sorted coulomb matrices of both molecules"""
        # Similar to .eigenCM_distance method 
        if self.natoms != other.natoms:
            if self.natoms < other.natoms:
                smaller = self
                larger = other
            else:
                smaller = other
                larger = self
           
            pad_width = larger.natoms - smaller.natoms
            self_sort_eigenval = np.pad(smaller.sortedCM(), (0, pad_width), 'constant')
            other_sort_eigenval = larger.sortedCM()
        else:
            self_sort_eigenval = self.sortedCM()
            other_sort_eigenval = other.sortedCM()

        distance = np.linalg.norm(self_sort_eigenval - other_sort_eigenval)
        return distance


