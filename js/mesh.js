// In this file you will implement traversal and analysis for your assignment.
// Make sure to familiarize yourself with the utility functions in meshUtils.js
// they might be useful for the second part of your assignment!

////////////////////////////////////////////////////////////////////////////////
// Traversal
////////////////////////////////////////////////////////////////////////////////

// Return all vertices on face f
Mesh.prototype.verticesOnFace = function(f) {
  const vertices = [];
  let he = f.halfedge;
  const first = he;
  while (true) {
    vertices.push(he.vertex);
    he = he.next;
    if (he === first) {
      break;
    }
  }
  return vertices;
};

// Return all halfedges on face f
Mesh.prototype.edgesOnFace = function(f) {
  const halfedges = [];

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  let he = f.halfedge;
  const first = he;
  while (true) {
      halfedges.push(he);
      he = he.next;
      if (he === first) {
          break;
      }
  }
  // ----------- STUDENT CODE END ------------

  return halfedges;
};

// Return all faces adjacent to input face, not
// including input face.
Mesh.prototype.facesOnFace = function(f) {
  const faces = [];

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  let he = f.halfedge;
  const first = he;
  while (true) {
      faces.push(he.opposite.face);
      he = he.next;
      if (he === first) {
          break;
      }
  }
  // ----------- STUDENT CODE END ------------

  return faces;
};

// Return one-ring neighbors of input vertex, not
// including the input vertex itself
Mesh.prototype.verticesOnVertex = function(v) {
  const vertices = [];

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  const original_he = v.halfedge;
  let he = original_he;
  do {
    vertices.push(he.vertex);
    he = he.opposite.next;
  } while (he != original_he);
  // ----------- STUDENT CODE END ------------

  return vertices;
};

// Return all halfedges that point away from v
Mesh.prototype.edgesOnVertex = function(v) {
  const halfedges = [];

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  const original_he = v.halfedge;
  let he = original_he;
  do {
    halfedges.push(he);
    he = he.opposite.next;
  } while (he != original_he);
  // ----------- STUDENT CODE END ------------

  return halfedges;
};

// Return all faces that include v as a vertex.
Mesh.prototype.facesOnVertex = function(v) {
  const faces = [];

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  const original_he = v.halfedge;
  let he = original_he;
  do {
    faces.push(he.face);
    he = he.opposite.next;
  } while (he != original_he);
  // ----------- STUDENT CODE END ------------

  return faces;
};

// Return the vertices that form the endpoints of a given edge
Mesh.prototype.verticesOnEdge = function(e) {
  const vertices = [];

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 2 lines of code.
  vertices.push(e.vertex);
  vertices.push(e.opposite.vertex);
  // ----------- STUDENT CODE END ------------

  return vertices;
};

// Return the faces that include a given edge
Mesh.prototype.facesOnEdge = function(e) {
  const faces = [];
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 2 lines of code.
  faces.push(e.face);
  faces.push(e.opposite.face);
  // ----------- STUDENT CODE END ------------
  return faces;
};

// Return the edge pointing from v1 to v2
Mesh.prototype.edgeBetweenVertices = function(v1, v2) {
  let out_he = undefined;
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 11 lines of code.
  const original_he = v1.halfedge;
  let he = original_he;
  do {
    if (he.vertex === v2) { // edges point from v1 to v2
      out_he = he;
      break;
    }
    he = he.opposite.next;
  } while (he != original_he);
  // ----------- STUDENT CODE END ------------
  return out_he;
};

////////////////////////////////////////////////////////////////////////////////
// Analysis
////////////////////////////////////////////////////////////////////////////////

// Return the surface area of a provided face f.
Mesh.prototype.calculateFaceArea = function(f) {
  let area = 0.0;
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 21 lines of code.
  const vertices = this.verticesOnFace(f);
  for (let i = 0; i < vertices.length - 2; i++) {
    area += triangleArea(vertices[0], vertices[i + 1], vertices[i + 2]);
  }
  // ----------- STUDENT CODE END ------------
  f.area = area;
  return area;
};

// helper function: calculate the area of a triangle with vertices va, vb, vc
function triangleArea(va, vb, vc) {
  let ab = new THREE.Vector3();
  let ac = new THREE.Vector3();
  ab = ab.subVectors(vb.position, va.position);
  ac = ac.subVectors(vc.position, va.position);
  let area = 0.5 * Math.sqrt((ab.y * ac.z - ab.z * ac.y) ** 2
  + (ab.z * ac.x - ab.x * ac.z) ** 2
  + (ab.x * ac.y - ab.y * ac.x) ** 2);
  return area;
}

// Update the area attributes of all faces in the mesh
Mesh.prototype.calculateFacesArea = function() {
  for (let i = 0; i < this.faces.length; ++i) {
    this.calculateFaceArea(this.faces[i]);
  }
};

// Calculate the vertex normal at a given vertex,
// using the face normals of bordering faces, weighted by face area
Mesh.prototype.calculateVertexNormal = function(v) {
  const v_normal = new THREE.Vector3();
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 11 lines of code.
  const original_he = v.halfedge;
  let he = original_he;
  let total_area = 0;
  do {
    let f_normal = CopyVec(he.face.normal);
    let area = he.face.area;
    v_normal.add(f_normal.multiplyScalar(area));
    total_area += area;
    he = he.opposite.next;
  } while (he != original_he);
  v_normal.multiplyScalar(1 / total_area);
  // ----------- STUDENT CODE END ------------
  v.normal = v_normal;
  return v_normal;
};

// update the vertex normals of every vertex in the mesh
Mesh.prototype.updateVertexNormals = function() {
  for (let i = 0; i < this.vertices.length; ++i) {
    this.calculateVertexNormal(this.vertices[i]);
  }
};

// compute the average length of edges touching v
Mesh.prototype.averageEdgeLength = function(v) {
  let avg = 0.0;

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  const original_he = v.halfedge;
  let he = original_he;
  let count = 0;
  do {
    avg += edgeLength(he);
    count++;
    he = he.opposite.next;
  } while (he != original_he);
  avg /= count;
  // ----------- STUDENT CODE END ------------
  return avg;
};

// helper function: calculating the length of an edge
function edgeLength(he) {
  let v1 = he.vertex;
  let v2 = he.opposite.vertex;
  let length = CopyVec(v1.position).sub(v2.position).length();
  return length;
}

////////////////////////////////////////////////////////////////////////////////
// Topology
////////////////////////////////////////////////////////////////////////////////

// Given a face in the shape of an arbitrary polygon,
// split that face so it consists only of several triangular faces.
Mesh.prototype.triangulateFace = function(f) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 8 lines of code.
  let vertices = this.verticesOnFace(f);
  if (vertices.length <= 3) return; // base case
  let newFace = this.splitFaceMakeEdge(f, vertices[0], vertices[2], vertices[1]);
  let m_faces = this.getModifiableFaces();
  if (m_faces.length !== this.faces.length) {
    newFace.selected = true;
  }
  this.triangulateFace(newFace);
  return;
  // ----------- STUDENT CODE END ------------
};
