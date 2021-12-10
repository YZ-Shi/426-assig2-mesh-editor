var Filters = Filters || {};

// Space for your helper functions
// ----------- STUDENT CODE BEGIN ------------
// ----------- Our reference solution uses 105 lines of code.

// find the face: v0 is one of its vertices, vL and vR are on two edges
function faceBetweenThreeVerts(v0, vL, vR) {
  let he = vL.halfedge;
  while (he.vertex !== v0) {
    he = he.opposite.next;
  }
  let face;
  if (he.next.vertex !== vR) {
    face = he.opposite.face;
  } else {
    face = he.face;
  }
  return face;
}

// test if a vertex is on a crease (of selection)
function detectCrease(mesh, v) {
  let n_selected = mesh.getModifiableFaces().length;
  if (n_selected === mesh.faces.length) return false; // no selection
  let faces = mesh.facesOnVertex(v);
  for (let i = 0; i < faces.length; i++) {
    if (!faces[i].selected) return true; // one of the faces touching the vertex is not selected
  }
  return false;
}

// calculate the length of a halfedge
/*function edgeLength(mesh, he) {
let v1 = he.vertex;
let v2 = he.opposite.vertex;
let length = CopyVec(v1.position).sub(v2.position).length();
return;
}*/
// ----------- STUDENT CODE END ------------

// Translate all selected vertices in the mesh by the given x,y,z offsets.
Filters.translation = function(mesh, x, y, z) {
  const t = new THREE.Vector3(x, y, z);

  const verts = mesh.getModifiableVertices();

  const n_vertices = verts.length;
  for (let i = 0; i < n_vertices; ++i) {
    verts[i].position.add(t);
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Given x,y,z, the desired rotation around each axis, in radians,
// apply this rotation to all selected vertices in the mesh.
Filters.rotation = function(mesh, x, y, z) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 5 lines of code.
  const r = new THREE.Euler(x, y, z, 'XYZ');

  const n_vertices = verts.length;
  for (let i = 0; i < n_vertices; i++) {
    verts[i].position.applyEuler(r);
  }
  // ----------- STUDENT CODE END ------------

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Uniformly scale the position of all selected vertices in the mesh
// by the provided scale factor s
Filters.scale = function(mesh, s) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 4 lines of code.
  const n_vertices = verts.length;
  for (let i = 0; i < n_vertices; i++) {
    verts[i].position.multiplyScalar(s);
  }
  // ----------- STUDENT CODE END ------------

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// estimate the per-vertex gaussian vurvature of the mesh at each vertex.
// set that vertex's color to some value based on its curvature value.
// (the precise mapping of curvature to color is left to you)
Filters.curvature = function(mesh) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 102 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("Curvature is not implemented yet");
};

// Apply a random offset to each selected vertex in the direction of its normal
// scale the random offset by the provided factor and by
// the average length of edges at that vertex
Filters.noise = function(mesh, factor) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 13 lines of code.
  const n_vertices = verts.length;
  let offsets = []; // stores the offsets
  // calculate the positions for updating
  for (let i = 0; i < n_vertices; i++) {
    let normal = CopyVec(verts[i].normal);
    let random = Math.random() * 2 - 1; // random offset (-1, 1)
    normal.multiplyScalar(random);
    normal.multiplyScalar(mesh.averageEdgeLength(verts[i]));
    normal.multiplyScalar(factor);
    offsets[i] = normal;
  }
  // update the positions
  for (let i = 0; i < n_vertices; i++) {
    verts[i].position.add(offsets[i]);
  }
  // ----------- STUDENT CODE END ------------

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Smooth the mesh using the specified weighting scheme.
// In the standard case, this is done using uniform Laplacian smoothing,
// by moving each vertex towards the average position of its neighbors.
//
// Arguments:
//  - mesh: the mesh to smooth
//  - iter: the number of iterations of smoothing to apply
//  - delta: a scaling factor for the amount of smoothing
//  - curvFlow: a bool. if true, use cotangent weights instead of uniform (requires triangular mesh)
//  - scaleDep: a bool. if true, scale offsets differently for each vertex (see spec.)
//  - implicit: a bool. if true, perform implicit smoothing (see spec.)
//
// Note that the reference solution calls a giant utility function so the line
// count is not terribly representative of the true solution
//
// For implicit, you will want to compute I - M*L*delta, where I is the identity
// matrix, M is a diagonal "mass" matrix, and L is a Laplacian matrix. Then
// you will want to call math.lup() on your result in order to decompose the
// matrix. Finally, call math.lusolve() to compute the X,Y, and Z positions of
// vertices. Note that the decomposition step allows for fast solving multiple
// times. It would be possible to replace a few of these steps with simple matrix
// inversion; however, matrix inversion is far slower and less numerically stable
//
Filters.smooth = function(mesh, iter, delta, curvFlow, scaleDep, implicit) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 16 lines of code.
  const n_vertices = verts.length;
  if (curvFlow) Filters.triangulate(mesh); // curvFlow only: triangulate the mesh if using curvature-flow weighting scheme

  for (let n = 0; n < iter; n++) { // n iterations
    let new_positions = []; // stores the new positions
    let areas = []; // scaleDep only: stores the area of neighboring faces for each vertex
    let total_area = 0; // scaleDep only: average area of neighboring faces from all vertices

    // if scale-dependent, calculate the area of neighboring faces for each vertex and the average of all vertices
    if (scaleDep) {
      for (let i = 0; i < n_vertices; i++) {
        let neighb_f = mesh.facesOnVertex(verts[i]);
        let area = 0;
        for (let j = 0; j < neighb_f.length; j++) {
          area += mesh.calculateFaceArea(neighb_f[j]);
        }
        areas[i] = area;
        total_area += area;
      }
      total_area /= n_vertices;
    }
    // compute positions
    for (let i = 0; i < n_vertices; i++) {
      const vertex = verts[i]; // current vertex
      let v_pos = CopyVec(vertex.position); // current vertex position
      let avg_pos = new THREE.Vector3(0, 0, 0);
      // weighting schemes
      if (curvFlow) { // curvature-flow Laplacian Smoothing
        let total_w = 0;
        const original_he = vertex.halfedge;
        let he = original_he;
        do { // calculate weight for each neighboring vertex
          let angle1 = mesh.angleBetweenEdges(he.opposite.next.vertex, he.opposite.next.opposite, he.opposite.next.next);
          let angle2 = mesh.angleBetweenEdges(he.next.vertex, he.next.opposite, he.next.next);
          let weight = 0.5 * (1 / Math.tan(angle1) + 1 / Math.tan(angle2));
          avg_pos.add(CopyVec(he.vertex.position).multiplyScalar(weight));
          total_w += weight;
          he = he.opposite.next;
        } while (he != original_he);
        avg_pos.sub(v_pos.multiplyScalar(total_w));

      } else { // uniform Laplacian Smoothing
        const neighbors = mesh.verticesOnVertex(vertex);
        for (let j = 0; j < neighbors.length; j++) {
          avg_pos.add(neighbors[j].position);
        }
        avg_pos.sub(v_pos.multiplyScalar(neighbors.length));
        avg_pos.multiplyScalar(1 / neighbors.length); // optional - normalize the weights to keep the overall dimensions
      }

      avg_pos.multiplyScalar(delta);
      // scale-dependent smoothing
      if (scaleDep) {
        avg_pos.multiplyScalar(total_area / areas[i]);
      }
      // add offset to original position
      let new_pos = CopyVec(vertex.position);
      new_pos.add(avg_pos);
      new_positions[i] = new_pos;
    }
    // update positions
    for (let i = 0; i < n_vertices; i++) {
      verts[i].position = CopyVec(new_positions[i]);
    }
  }
  // ----------- STUDENT CODE END ------------
  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Sharpen the mesh by moving selected vertices away from the average position
// of their neighbors (i.e. Laplacian smoothing in the negative direction)
Filters.sharpen = function(mesh, iter, delta) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  const n_vertices = verts.length;

  for (let n = 0; n < iter; n++) { // for nth interation

    let new_positions = []; // stores the new positions
    for (let i = 0; i < n_vertices; i++) { // compute positions for each vertex
      const vertex = verts[i];
      let v_pos = CopyVec(vertex.position); // the original position
      const neighbors = mesh.verticesOnVertex(vertex);
      // compute the average offset of the neighboring vertices
      let avg_pos = new THREE.Vector3(0, 0, 0);
      for (let j = 0; j < neighbors.length; j++) {
        avg_pos.add(neighbors[j].position);
      }
      avg_pos.sub(v_pos.multiplyScalar(neighbors.length));
      avg_pos.multiplyScalar(delta);

      let new_pos = CopyVec(vertex.position);
      new_pos.sub(avg_pos); // subtract the average position instead of adding, moving the vertex away from the average position
      new_positions[i] = new_pos;
    }
    // update positions
    for (let i = 0; i < n_vertices; i++) {
      verts[i].position = CopyVec(new_positions[i]);
    }
  }
  // ----------- STUDENT CODE END ------------
  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Move every selected vertex along its normal direction
// Scale the amount by the provided factor and average edge length at that vertex
Filters.inflate = function(mesh, factor) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 16 lines of code.
  const n_vertices = verts.length;
  let scaleFactor = []; // stores the positions for updating
  // calculate the positions
  for (let i = 0; i < n_vertices; i++) {
    let normal = CopyVec(verts[i].normal);
    normal.multiplyScalar(factor);
    normal.multiplyScalar(mesh.averageEdgeLength(verts[i])); // optional augmentation with average edge length
    scaleFactor[i] = normal;
  }
  // update the positions
  for (let i = 0; i < n_vertices; i++) {
    verts[i].position.add(scaleFactor[i]);
  }
  // ----------- STUDENT CODE END ------------

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// rotate selected vertices around the Y axis by an amount
// proportional to its Y value times the scale factor.
Filters.twist = function(mesh, factor) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 8 lines of code.
  const n_vertices = verts.length;
  for (let i = 0; i < n_vertices; i++) {
    let rotate = verts[i].position.y * factor;
    let r = new THREE.Euler(0, rotate, 0, 'XYZ');
    verts[i].position.applyEuler(r);
  }
  // ----------- STUDENT CODE END ------------

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// warp a mesh using a nonlinear mapping of your choice
Filters.wacky = function(mesh, factor) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 3 lines of code.
  for (let i = 0; i < mesh.vertices.length; i++) {
    mesh.vertices[i].position.x *= Math.abs(mesh.vertices[i].position.y);
    mesh.vertices[i].position.x *= factor;
    mesh.vertices[i].position.z *= Math.abs(mesh.vertices[i].position.y);
    mesh.vertices[i].position.z *= factor;
  }
  // ----------- STUDENT CODE END ------------

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Convert the selected faces from arbitrary polygons into all triangles
Filters.triangulate = function(mesh) {
  const faces = mesh.getModifiableFaces();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 4 lines of code.
  const n_faces = faces.length;
  for (let i = 0; i < n_faces; i++) {
    mesh.triangulateFace(faces[i]);
  }
  // ----------- STUDENT CODE END ------------

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for splitEdgeMakeVert in mesh.js
Filters.splitEdge = function(mesh) {
  const verts = mesh.getSelectedVertices();

  if (verts.length === 2) {
    mesh.splitEdgeMakeVert(verts[0], verts[1], 0.5);
  } else {
    console.log("ERROR: to use split edge, select exactly 2 adjacent vertices");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for joinEdgeKillVert in mesh.js
Filters.joinEdges = function(mesh) {
  const verts = mesh.getSelectedVertices();

  if (verts.length === 3) {
    let v0 = verts[0],
    v1 = verts[1],
    v2 = verts[2];

    const he01 = mesh.edgeBetweenVertices(v0, v1);
    const he12 = mesh.edgeBetweenVertices(v1, v2);

    if (he01) {
      if (he12) {
        mesh.joinEdgeKillVert(verts[0], verts[1], verts[2]);
      } else {
        mesh.joinEdgeKillVert(verts[1], verts[0], verts[2]);
      }
    } else {
      if (he12) {
        mesh.joinEdgeKillVert(verts[0], verts[2], verts[1]);
      } else {
        console.log(
          "ERROR: to use join edge, select exactly 3 vertices such that one only has edges to the other two"
        );
      }
    }
  } else {
    console.log("ERROR: to use join edge, select exactly 3 vertices");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for splitFaceMakeEdge in mesh.js
Filters.splitFace = function(mesh) {
  const verts = mesh.getSelectedVertices();
  const faces = mesh.getModifiableFaces();

  if (verts.length === 2 && faces.length === 1) {
    mesh.splitFaceMakeEdge(faces[0], verts[0], verts[1]);
  } else {
    console.log("ERROR: to use split face, select exactly 1 face and 2 nonadjacent vertices on it");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for joinFaceKillEdge in mesh.js
Filters.joinFaces = function(mesh) {
  const verts = mesh.getSelectedVertices();
  const faces = mesh.getModifiableFaces();

  if (verts.length === 2 && faces.length === 2) {
    mesh.joinFaceKillEdge(faces[0], faces[1], verts[0], verts[1]);
  } else {
    console.log(
      "ERROR: to use split face, select exactly 2 adjacent faces the 2 vertices between them"
    );
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// extrude the selected faces from the mesh in the direction of their normal
// vector, scaled by the provided factor.
// See the spec for more detail.
Filters.extrude = function(mesh, factor) {
  const faces = mesh.getModifiableFaces();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 32 lines of code.
  const n_oldFaces = faces.length; // number of old faces
  for (let n = 0; n < n_oldFaces; n++) { // for each old face
    const oldVerts = mesh.verticesOnFace(faces[n]); // vertices on the old face
    const newVerts = [];
    const newFaces = [];
    // create new vertices
    for (let i = 0; i < oldVerts.length; i++) {
      if (i + 1 !== oldVerts.length) {
        newVerts.push(mesh.splitEdgeMakeVert(oldVerts[i], oldVerts[i + 1], 0));
      } else {
        newVerts.push(mesh.splitEdgeMakeVert(oldVerts[i], oldVerts[0], 0));
      }
    }
    // connect the old vertices, split the neighboring faces
    for (let i = 0; i < oldVerts.length; i++) {
      let f;
      let e = mesh.edgeBetweenVertices(oldVerts[i], newVerts[i]);
      if (e.face === faces[n]) {
        f = e.opposite.face;
      } else {
        f = e.face;
      }
      if (i + 1 !== oldVerts.length) {
        newFaces.push(mesh.splitFaceMakeEdge(f, oldVerts[i], oldVerts[i + 1], newVerts[i], true));
      } else {
        newFaces.push(mesh.splitFaceMakeEdge(f, oldVerts[i], oldVerts[0], newVerts[i], true));
      }
    }
    // connect the new vertices, split current face
    for (let i = 0; i < newVerts.length; i++) {
      if (i + 1 !== newVerts.length) {
        newFaces.push(mesh.splitFaceMakeEdge(faces[n], newVerts[i], newVerts[i + 1], oldVerts[i + 1], true));
      } else {
        newFaces.push(mesh.splitFaceMakeEdge(faces[n], newVerts[i], newVerts[0], oldVerts[0], true));
      }
    }
    // join each pair of two new faces
    for (let i = 0; i < newVerts.length; i++) {
      if (i + 1 !== newVerts.length) {
        mesh.joinFaceKillEdge(newFaces[i], newFaces[newVerts.length + i], newVerts[i], oldVerts[i + 1]);
      } else {
        mesh.joinFaceKillEdge(newFaces[i], newFaces[newVerts.length + i], newVerts[i], oldVerts[0]);
      }
    }
    // geometrical Move
    for (let i = 0; i < newVerts.length; i++) {
      let normal = CopyVec(faces[n].normal);
      newVerts[i].position.add(normal.multiplyScalar(factor));
    }
  }
  // ----------- STUDENT CODE END ------------

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Truncate the selected vertices of the mesh by "snipping off" corners
// and replacing them with faces. factor specifies the size of the truncation.
// See the spec for more detail.
Filters.truncate = function(mesh, factor) {
  const verts = mesh.getModifiableVertices();
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 64 lines of code.
  let n_oldVerts = verts.length;
  let new_pos = []; // new positions
  let newFaces = [];
  // topological operations
  for (let n = 0; n < n_oldVerts; n++) {
    const oldVert = verts[n]; // the current vertex
    let neighbors = mesh.verticesOnVertex(oldVert); // one-ring neighbors of the current vertex
    // create n-1 new vertices
    let newVerts = [];
    for (let i = 0; i < neighbors.length - 1; i++) {
      newVerts.push(mesh.splitEdgeMakeVert(oldVert, neighbors[i], 0));
    }
    // find face between every 3 vertices
    let oldFace = faceBetweenThreeVerts(oldVert, newVerts[0], newVerts[1]);
    // make a new edge
    for (let i = 0; i < newVerts.length - 1; i++) {
      newFaces.push(mesh.splitFaceMakeEdge(oldFace, newVerts[i], newVerts[i + 1]));
    }
    let faces = mesh.getModifiableFaces();
    if (faces !== mesh.faces) {
      for (let j = 0; j < newFaces.length; j++) {
        newFaces[j].selected = true;
      }
    }
    // calculate offsets - derived from meshUtils.splitEdgeMakeVert
    let offsets = [];
    for (let i = 0; i < neighbors.length; i++) {
      const offset = new THREE.Vector3();
      const p1 = CopyVec(oldVert.position);
      const p2 = CopyVec(neighbors[i].position);
      offset.add(p1.multiplyScalar(1 - factor));
      offset.add(p2.multiplyScalar(factor));
      offsets.push(offset);
    }
    new_pos[oldVert.id] = offsets[offsets.length - 1];
    for (let i = 0; i < offsets.length - 1; i++) {
      new_pos[newVerts[i].id] = offsets[i];
    }
  }
  // update the new positions
  for (let n = 0; n < new_pos.length; n++) {
    if (new_pos[n] === undefined) continue;
    mesh.vertices[n].position = new_pos[n];
  }
  // ----------- STUDENT CODE END ------------

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Apply the bevel operation to the mesh, scaling the degree of bevelling by factor
Filters.bevel = function ( mesh, factor ) {

  var verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 104 lines of code.
  let faces = mesh.getModifiableFaces();
  const n_faces = faces.length;
  Filters.truncate(mesh, 0); // truncate the vertices
  verts = mesh.getModifiableVertices();
  faces = mesh.getModifiableFaces();
  const n_verts = mesh.vertices.length; // number of vertices after truncate
  for (let i = n_faces; i < faces.length; i++) { // for faces created in truncate
    let oldVerts = mesh.verticesOnFace(faces[i]);
    for (let j = 0; j < oldVerts.length; j++) {
      let next = (j === oldVerts.length - 1) ? 0 : j + 1;
      mesh.splitEdgeMakeVert(oldVerts[j], oldVerts[next], 0.5);
    }
  }
  const n_facesTrunc = faces.length; // number of faces after truncate
  let newFaces = [];
  // connect the new vertices on surfaces created in truncate
  for (let i = 0; i < n_faces; i++) {
    let currentVerts = mesh.verticesOnFace(faces[i]);
    for (let j = 0; j < currentVerts.length; j++) {
      if (currentVerts[j].id >= n_verts) {
        let next = (j + 3 >= currentVerts.length) ? j + 3 - currentVerts.length : j + 3;
        let ref = (j + 1 === currentVerts.length) ? 0 : j + 1;
        newFaces.push(mesh.splitFaceMakeEdge(faces[i], currentVerts[j], currentVerts[next], currentVerts[ref], true));
      }
    }
  }
  // join new faces along the same edge
  for (let i = 0; i < newFaces.length; i++) {
    let edges = mesh.edgesOnFace(newFaces[i]);
    for (let j = 0; j < edges.length; j++) {
      if (edges[j].opposite.face.id > n_facesTrunc) {
        mesh.joinFaceKillEdgeSimple(edges[j]);
      }
    }
  }
  // kill extra vertices
  for (let i = n_faces; i < n_facesTrunc; i++) {
    let currentVerts = mesh.verticesOnFace(faces[i]);
    for (let j = 0; j < currentVerts.length; j++) {
      if (currentVerts[j].id >= n_verts) {
        let nextOld = (j + 1 === currentVerts.length) ? 0 : j + 1;
        let nextNew = (j + 2 >= currentVerts.length) ? j + 2 - currentVerts.length : j + 2;
        mesh.joinEdgeKillVert(currentVerts[j], currentVerts[nextOld], currentVerts[nextNew]);
      }
    }
  }
  // geometrically move vertices
  for (let i = 0; i < n_faces; i++) {
    const centroid = mesh.calculateFaceCentroid(faces[i]);
    let currentVerts = mesh.verticesOnFace(faces[i]);
    for (let j = 0; j < currentVerts.length; j++) {
      const newPos = new THREE.Vector3();
      const p1 = CopyVec(currentVerts[j].position);
      const p2 = CopyVec(centroid);
      newPos.add(p1.multiplyScalar(1 - factor));
      newPos.add(p2.multiplyScalar(factor));
      currentVerts[j].position = newPos;
    }
  }

  // ----------- STUDENT CODE END ------------

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Split the longest edges in the mesh into shorter edges.
// factor is a float in [0,1]. it tells the proportion
// of the total number of edges in the mesh that should be split.
Filters.splitLong = function(mesh, factor) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 35 lines of code.
  // determine number of splits
  const f = mesh.getModifiableFaces();
  let he_id = [];
  for (let n = 0; n < f.length; n++) {
    let he = mesh.edgesOnFace(f[n]);
    for (let k = 0; k < he.length; k++) {
      if (!he_id.includes(he[k].id)) he_id.push(he[k].id);
      if (!he_id.includes(he[k].opposite.id)) he_id.push(he[k].opposite.id);
    }
  }
  const n_edges = he_id.length / 2; // number of edges in the selection

  const iter = n_edges * factor;
  for (let n = 0; n < iter; n++) {
    const faces = mesh.getModifiableFaces();
    // find the longest edge
    let longest;
    let max_length = Number.NEGATIVE_INFINITY;
    for (let i = 0; i < faces.length; i++) {
      let edges = mesh.edgesOnFace(faces[i]);
      //debugger;
      for (let j = 0; j < edges.length; j++) {
        let length = edgeLength(edges[j]);
        if (length > max_length) {
          longest = edges[j];
          max_length = length;
        }
      }
    }
    // find the neighboring faces and vertices
    let f1 = longest.face;
    let f2 = longest.opposite.face;
    let v1 = longest.next.vertex;
    let v2 = longest.opposite.next.vertex;
    // split the edge
    let newVert = mesh.splitEdgeMakeVert(longest.vertex, longest.opposite.vertex, 0.5);
    let newFace1 = mesh.splitFaceMakeEdge(f1, v1, newVert);
    let newFace2 = mesh.splitFaceMakeEdge(f2, v2, newVert);
    if (faces.length !== mesh.faces.length) {
      if (f1.selected) newFace1.selected = true;
      if (f2.selected) newFace2.selected = true;
    }
  }
  // ----------- STUDENT CODE END ------------

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Triangulate a mesh, and apply triangular subdivision to its faces.
// Repeat for the specified number of levels.
Filters.triSubdiv = function(mesh, levels) {
  Filters.triangulate(mesh);

  for (let l = 0; l < levels; l++) {
    const faces = mesh.getModifiableFaces();
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 43 lines of code.
    const vertices = mesh.getModifiableVertices();
    const n_faces = faces.length;
    const n_verts = vertices.length;
    let v_id = []; // ids of modifiable vertices
    for (let i = 0; i < n_verts; i++) {
      v_id.push(vertices[i].id);
    }
    // split edges and add new vertices
    for (let i = 0; i < n_faces; i++) {
      let edges = mesh.edgesOnFace(faces[i]);
      for (let j = 0; j < edges.length; j++) {
        if (!v_id.includes(edges[j].vertex.id) || !v_id.includes(edges[j].opposite.vertex.id)) {
          continue; // avoid double spliting
        }
        mesh.splitEdgeMakeVert(edges[j].vertex, edges[j].opposite.vertex, 0.5);
      }
    }
    // join new vertices
    for (let i = 0; i < n_faces; i++) {
      let verts = mesh.verticesOnFace(faces[i]);
      let newFaces = [];
      for (let j = 0; j < verts.length; j++) {
        if (!v_id.includes(verts[j].id)) { // connect each new vertex with the next vertex
          let nextNew, oldVert;
          if (j + 2 === verts.length) {
            nextNew = 0;
            oldVert = verts.length - 1;
          } else if (j + 1 === verts.length) {
            nextNew = 1;
            oldVert = 0;
          } else {
            nextNew = j + 2;
            oldVert = j + 1;
          }
          newFaces.push(mesh.splitFaceMakeEdge(faces[i], verts[j], verts[nextNew], verts[oldVert], true));
        }
      }
      // propagate the selection if the mesh is partially selected
      if (faces !== mesh.faces) {
        for (let j = 0; j < newFaces.length; j++) {
          newFaces[j].selected = true;
        }
      }
    }
    // ----------- STUDENT CODE END ------------
  }
  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Triangulate the mesh and apply loop subdivision to the faces
// repeat for the specified number of levels.
Filters.loop = function(mesh, levels) {
  Filters.triangulate(mesh);

  for (let l = 0; l < levels; l++) {
    const faces = mesh.getModifiableFaces();
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 123 lines of code.
    let partial = false; // selected?
    if (faces.length < mesh.faces.length) partial = true;

    const vertices = mesh.getModifiableVertices();
    const n_meshVerts = mesh.vertices.length; // number of the complete list of vertices
    const n_faces = faces.length; // number of (selected) faces
    const n_verts = vertices.length; // number of (selected) vertices
    let v_id = []; // ids of (selected) (even) vertices
    for (let i = 0; i < n_verts; i++) {
      v_id.push(vertices[i].id);
    }

    if (partial) {
      // keep track of crease cases - even vertices
      var crease = new Array(mesh.vertices.length);
      crease.fill(false);
      for (let i = 0; i < n_verts; i++) {
        crease[vertices[i].id] = detectCrease(mesh, vertices[i]);
      }
      // find faces along the boundary
      var neighborFaces = new Array(mesh.faces.length);
      neighborFaces.fill(false);
      for (let i = 0; i < crease.length; i++) {
        if (crease[i]) {
          let f = mesh.facesOnVertex(mesh.vertices[i]);
          for (let j = 0; j < f.length; j++) {
            if (!faces.includes(f[j])) neighborFaces[f[j].id] = true;
          }
        }
      }
    }
    // find the even neighbors of even vertices
    let evenNeighbors = [];
    for (let i = 0; i < n_verts; i++) { // ids of even vertices
      let neighbors = mesh.verticesOnVertex(vertices[i]);
      if (partial && crease[vertices[i].id]) { // crease cases
        evenNeighbors[i] = [];
        for (let j = 0; j < neighbors.length; j++) {
          let he = mesh.edgeBetweenVertices(vertices[i], neighbors[j]);
          if (he.face.selected !== he.opposite.face.selected) {  // edge is on the boundary
            evenNeighbors[i].push(neighbors[j]);
          }
        }
      } else {
        evenNeighbors[i] = neighbors;
      }
    }
    // split the edge and create odd vertices
    let oddVerts = [];
    for (let i = 0; i < n_faces; i++) {
      let edges = mesh.edgesOnFace(faces[i]);
      for (let j = 0; j < edges.length; j++) { // split each edge at the midpoint
        if (!v_id.includes(edges[j].vertex.id) || !v_id.includes(edges[j].opposite.vertex.id)) continue;
        oddVerts.push(mesh.splitEdgeMakeVert(edges[j].vertex, edges[j].opposite.vertex, 0.5));
      }
    }
    if (partial) {
      // trisub neighboring faces
      for (let i = 0; i < neighborFaces.length; i++) {
        if (neighborFaces[i]) {
          let edges = mesh.edgesOnFace(mesh.faces[i]);
          for (let j = 0; j < edges.length; j++) { // split each edge at the midpoint
            if (edges[j].vertex.id >= n_meshVerts || edges[j].opposite.vertex.id >= n_meshVerts) continue;
            mesh.splitEdgeMakeVert(edges[j].vertex, edges[j].opposite.vertex, 0.5);
          }
        }
      }
      // keep track of crease cases - odd vertices
      for (let i = 0; i < oddVerts.length; i++) {
        crease[oddVerts[i].id] = detectCrease(mesh, oddVerts[i]);
      }
    }
    // masked weight for even vertices
    let evenPos = [];
    for (let i = 0; i < n_verts; i++) {
      evenPos[i] = new THREE.Vector3();
      let k = evenNeighbors[i].length;
      if (partial && crease[vertices[i].id]) { // boundary/crease
        for (let j = 0; j < k; j++) {
          if (crease[evenNeighbors[i][j].id]) {
            evenPos[i].add(CopyVec(evenNeighbors[i][j].position).multiplyScalar(0.125));
          }
        }
        evenPos[i].add(CopyVec(vertices[i].position).multiplyScalar(0.75));
        continue;
      }
      if (k === 3) {
        for (let j = 0; j < k; j++) {
          evenPos[i].add(CopyVec(evenNeighbors[i][j].position).multiplyScalar(3/16));
        }
        evenPos[i].add(CopyVec(vertices[i].position).multiplyScalar(1 - 3/16 * k));
      }
      if (k > 3) {
        for (let j = 0; j < k; j++) {
          evenPos[i].add(CopyVec(evenNeighbors[i][j].position).multiplyScalar(0.375 / k));
        }
        evenPos[i].add(CopyVec(vertices[i].position).multiplyScalar(0.625));
      }
    }
    // calculate new positions for odd verticesOnFace
    let oddPos = [];
    for (let i = 0; i < oddVerts.length; i++) {
      oddPos[i] = new THREE.Vector3();
      if (partial && crease[oddVerts[i].id]) {
        oddPos[i].add(CopyVec(oddVerts[i].halfedge.vertex.position).multiplyScalar(0.5));
        oddPos[i].add(CopyVec(oddVerts[i].halfedge.opposite.next.vertex.position).multiplyScalar(0.5));
        continue;
      }
      oddPos[i].add(CopyVec(oddVerts[i].halfedge.vertex.position).multiplyScalar(0.375));
      oddPos[i].add(CopyVec(oddVerts[i].halfedge.opposite.next.vertex.position).multiplyScalar(0.375));
      oddPos[i].add(CopyVec(oddVerts[i].halfedge.next.next.vertex.position).multiplyScalar(0.125));
      oddPos[i].add(CopyVec(oddVerts[i].halfedge.opposite.next.next.next.vertex.position).multiplyScalar(0.125));
    }
    // join new vertices
    for (let i = 0; i < n_faces; i++) {
      let verts = mesh.verticesOnFace(faces[i]);
      let newFaces = [];
      for (let j = 0; j < verts.length; j++) {
        if (!v_id.includes(verts[j].id)) {
          let nextNew, oldVert;
          if (j + 2 === verts.length) {
            nextNew = 0;
            oldVert = verts.length - 1;
          } else if (j + 1 === verts.length) {
            nextNew = 1;
            oldVert = 0;
          } else {
            nextNew = j + 2;
            oldVert = j + 1;
          }
          newFaces.push(mesh.splitFaceMakeEdge(faces[i], verts[j], verts[nextNew], verts[oldVert], true));
        }
      }
      if (partial) { // update selection
        for (let j = 0; j < newFaces.length; j++) {
          newFaces[j].selected = true;
        }
      }
    }
    // divide neighboring faces
    if (partial) {
      for (let i = 0; i < neighborFaces.length; i++) {
        if (neighborFaces[i]) {
          let verts = mesh.verticesOnFace(mesh.faces[i]);
          for (let j = 0; j < verts.length; j++) {
            if (verts[j].id >= n_meshVerts) { // a new vertice
              let nextNew, oldVert;
              if (j + 2 === verts.length) {
                nextNew = 0;
                oldVert = verts.length - 1;
              } else if (j + 1 === verts.length) {
                nextNew = 1;
                oldVert = 0;
              } else {
                nextNew = j + 2;
                oldVert = j + 1;
              }
              mesh.splitFaceMakeEdge(mesh.faces[i], verts[j], verts[nextNew], verts[oldVert], true);
            }
          }
        }
      }
    }
    // update vertice new_positions
    for (let i = 0; i < vertices.length; i++) {
      if (i < n_verts) { // even vertex
        vertices[i].position = evenPos[i];
      } else { // odd vertex
        vertices[i].position = oddPos[i - n_verts];
      }
    }
    // ----------- STUDENT CODE END ------------
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Requires a quad mesh. Apply quad subdivision to the faces of the mesh.
// Repeat for the specified number of levels.
Filters.quadSubdiv = function(mesh, levels) {
  for (let l = 0; l < levels; l++) {
    const faces = mesh.getModifiableFaces();
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 55 lines of code.
    const vertices = mesh.getModifiableVertices();
    const n_faces = faces.length;
    const n_verts = vertices.length;
    let v_id = []; // ids of modifiable vertices
    for (let i = 0; i < n_verts; i++) {
      v_id.push(vertices[i].id);
    }
    let f_centroids = []; // store existing face centroids
    for (let i = 0; i < n_faces; i++) {
      f_centroids[i] = mesh.calculateFaceCentroid(faces[i]);
    }
    let newV = []; // new vertices
    // split edges and add new vertices
    for (let i = 0; i < n_faces; i++) {
      let newVerts = [];
      let edges = mesh.edgesOnFace(faces[i]);
      for (let j = 0; j < edges.length; j++) {
        if (!v_id.includes(edges[j].vertex.id) || !v_id.includes(edges[j].opposite.vertex.id)) {
          continue; // avoid double spliting
        }
        newVerts.push(mesh.splitEdgeMakeVert(edges[j].vertex, edges[j].opposite.vertex, 0.5));
      }
      newV[i] = newVerts;
    }
    let centroids = []; // new face centroids
    // join new vertices
    for (let i = 0; i < n_faces; i++) {
      let verts = mesh.verticesOnFace(faces[i]);
      let newFaces = [];
      for (let j = 0; j < verts.length; j++) { // join two new vertices and split the new edge at midpoint
        if (!v_id.includes(verts[j].id)) {
          let nextNew = j + 2;
          let oldVert = j + 1;
          newFaces.push(mesh.splitFaceMakeEdge(faces[i], verts[j], verts[nextNew], verts[oldVert], true));
          centroids[i] = mesh.splitEdgeMakeVert(verts[j], verts[nextNew], 0.5); // midpoint of the new edge
          break;
        }
      }
      for (let j = 0; j < verts.length; j++) {
        if (!v_id.includes(verts[j].id) && mesh.edgeBetweenVertices(centroids[i], verts[j]) === undefined) { // connect each new vertex with the next vertex
          let oldVert;
          if (j + 1 === verts.length) {
            oldVert = 0;
          } else {
            oldVert = j - 1;
          }
          newFaces.push(mesh.splitFaceMakeEdge(faces[i], centroids[i], verts[j], verts[oldVert], true));
        }
      }
      // propagate the selection if the mesh is partially selected
      if (faces !== mesh.faces) {
        for (let j = 0; j < newFaces.length; j++) {
          newFaces[j].selected = true;
        }
      }
    }
    // geometric transformation
    for (let i = 0; i < centroids.length; i++) {
      centroids[i].position = f_centroids[i];
    }
    // ----------- STUDENT CODE END ------------
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Apply catmull clark subdivision to the faces of the mesh.
// Repeat for the specified number of levels.
Filters.catmullClark = function(mesh, levels) {
  for (let l = 0; l < levels; l++) {
    const faces = mesh.faces;
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 102 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce("Catmull-Clark subdivide is not implemented yet");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// ================= internal functions =======================

// internal function for selecting faces in the form of a loop
Filters.procFace = function(mesh, f) {
  const faceFlags = new Array(mesh.faces.length);
  for (let i = 0; i < mesh.faces.length; i++) {
    faceFlags[i] = 0;
  }
  let sum = f.area;
  const start_he = f.halfedge.opposite.next;
  let curr_he = start_he;
  do {
    if (faceFlags[curr_he.face.id] > 0) {
      break;
    }
    sum += curr_he.face.area;
    curr_he.face.selected = true;
    faceFlags[curr_he.face.id]++;
    const last_he = curr_he;
    curr_he = curr_he.opposite.next;
    if (curr_he.face == f) {
      curr_he = last_he.next.opposite.next;
    }
  } while (true);
};

Filters.parseSelected = function(sel) {
  if (sel === undefined || sel.replace === undefined) {
    return [];
  }
  if (typeof sel === "number") {
    return [sel];
  }
  // sel = sel.replace(/[\(\)]/g,'');
  sel = sel.split(",");
  const parsedSel = [];
  for (let i = 0; i < sel.length; i++) {
    const idx = parseInt(sel[i]);
    if (!isNaN(idx)) {
      parsedSel.push(idx);
    }
  }
  return parsedSel;
};

// internal filter for updating selection
Filters.selection = function(mesh, vertIdxs, faceIdxs) {
  mesh.setSelectedVertices(Filters.parseSelected(vertIdxs));
  mesh.setSelectedFaces(Filters.parseSelected(faceIdxs));
};

// internal filter for setting display settings
Filters.displaySettings = function(
  mesh,
  showLabels,
  showHalfedge,
  shading,
  showVN,
  showFN,
  showGrid,
  showVertDots,
  showAxes,
  showVC,
  meshColor
) {
  Main.displaySettings.showIdLabels = showLabels;
  Main.displaySettings.wireframe = showHalfedge;
  Main.displaySettings.shading = shading;
  Main.displaySettings.showVN = showVN;
  Main.displaySettings.showFN = showFN;
  Main.displaySettings.showGrid = showGrid;
  Main.displaySettings.showVertDots = showVertDots;

  Main.displaySettings.showAxes = showAxes;
  Main.displaySettings.showVC = showVC;
  // Main.displaySettings.meshColor = meshColor;

  // Main.refreshDisplaySettings();
};
