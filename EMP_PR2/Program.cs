using EMP_PR2;

const string spaceMeshInput = "Meshes/SpaceMesh.txt";
const string timeMeshInput = "Meshes/TimeMesh.txt";

FEM fem1 = new(new Mesh(spaceMeshInput), new Mesh(timeMeshInput), new Test2(), NonlinearMethod.SIMPLE_ITER, 1e-15, 1000);
fem1.Compute();

FEM fem2 = new(new Mesh(spaceMeshInput), new Mesh(timeMeshInput), new Test2(), NonlinearMethod.NEWTON, 1e-15, 1000);
fem2.Compute();