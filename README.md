# global_intersection_analysis_sgi2024

An implementation of Baraff et al.'s "Untangling Cloth" algorithm, Global Intersection Analysis which focuse on resolving the intersection between two meshes which are forced to intersect each other at the beginning of a simulation. This was an outcome of a 1 week project by Sara Samy and Sachin Kishan under the mentorship of Zachary Ferguson under the MIT SGI 2024 program.

## To run
1. Clone the repository
2. Create a "build" subdirectory.
3. Run cmake
4. This should create a project build for you to run the code yourself depending on your OS with all dependencies.


## Screenshots

### Cloth "tangled" in the mesh
![Screenshot 2024-08-26 185107](https://github.com/user-attachments/assets/4482128d-85e8-46e4-9e28-7e976f94f66b)

### Cloth slowly being resolved to be removed from the mesh 
![Screenshot 2024-08-26 185040](https://github.com/user-attachments/assets/c9ed0f91-eb8a-4e05-b3e1-377231ba02e9)

### Cloth fully removed from the mesh, based off of the contour of intersection between both meshes
![Screenshot 2024-08-26 185132](https://github.com/user-attachments/assets/981d88c0-cde6-4087-b8a9-9dad3721f682)
