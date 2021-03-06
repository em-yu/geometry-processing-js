let Vertex = require('./vertex.js');
let Edge = require('./edge.js');
let Face = require('./face.js');
let Halfedge = require('./halfedge.js');
let Corner = require('./corner.js');

class Mesh {
	/**
	 * This class represents a Mesh.
	 * @constructor module:Core.Mesh
	 * @property {module:Core.Vertex[]} vertices The vertices contained in this mesh.
	 * @property {module:Core.Edge[]} edges The edges contained in this mesh.
	 * @property {module:Core.Face[]} faces The faces contained in this mesh.
	 * @property {module:Core.Corner[]} corners The corners contained in this mesh.
	 * @property {module:Core.Halfedge[]} halfedges The halfedges contained in this mesh.
	 * @property {module:Core.Face[]} boundaries The boundary loops contained in this mesh.
	 * @property {Array.<module:Core.Halfedge[]>} generators An array of halfedge arrays, i.e.,
	 * [[h11, h21, ..., hn1], [h12, h22, ..., hm2], ...] representing this mesh's
	 * {@link https://en.wikipedia.org/wiki/Homology_(mathematics)#Surfaces homology generators}.
	 */
	constructor() {
		/** @type {Vertex[]} */
		this.vertices = [];
		/** @type {Edge[]} */
		this.edges = [];
		/** @type {Face[]} */
		this.faces = [];
		/** @type {Corner[]} */
		this.corners = [];
		/** @type {Halfedge[]} */
		this.halfedges = [];
		/** @type {Face[]} */
		this.boundaries = [];
		/** @type {Array.Halfedge[]} */
		this.generators = [];
	}

	/**
	 * Creates a vertex and adds it to the mesh
	 * @method
	 * @returns {Vertex}
	 */
	newVertex() {
		let v = new Vertex();
		const vi = this.vertices.length;
		v.index = vi;
		this.vertices.push(v);
		return v;
	}

	/**
	 * Creates an edge and adds it to the mesh
	 * @method
	 * @returns {Edge}
	 */
	newEdge() {
		let e = new Edge();
		const ei = this.edges.length;
		e.index = ei;
		this.edges.push(e);
		return e;
	}

	/**
	 * Creates a face and adds it to the mesh
	 * @method
	 * @returns {Face}
	 */
	newFace() {
		let f = new Face();
		const fi = this.faces.length;
		f.index = fi;
		this.faces.push(f);
		return f;
	}

	/**
	 * Creates a corner and adds it to the mesh
	 * @method
	 * @returns {Corner}
	 */
	newCorner() {
		let c = new Corner();
		const ci = this.corners.length;
		c.index = ci;
		this.corners.push(c);
		return c;
	}

	/**
	 * Creates a halfedge and adds it to the mesh
	 * @method
	 * @returns {Halfedge}
	 */
	newHalfedge() {
		let h = new Halfedge();
		const hi = this.halfedges.length;
		h.index = hi;
		this.halfedges.push(h);
		return h;
	}

	/**
	 * Links two halfedges as prev and next of each other
	 * @method
	 * @param {Halfedge} he1 the halfedge before he2
	 * @param {Halfedge} he2 the halfedge after he1
	 */
	linkHalfedges(he1, he2) {
		he1.next = he2;
		he2.prev = he1;
	}

	/**
	 * Glues two halfedges as twins of each other
	 * @method
	 * @param {Halfedge} he1 twin halfedge of he2
	 * @param {Halfedge} he2 twin halfedge of he1
	 */
	glueHalfedges(he1, he2) {
		he1.twin = he2;
		he2.twin = he1;
	} 

	/**
	 * Computes the euler characteristic of this mesh.
	 * @method module:Core.Mesh#eulerCharacteristic
	 * @returns {number}
	 */
	eulerCharacteristic() {
		return this.vertices.length - this.edges.length + this.faces.length;
	}

	/**
	 * Constructs this mesh.
	 * @method module:Core.Mesh#build
	 * @param {Object} polygonSoup A polygon soup mesh containing vertex positions and indices.
	 * @param {module:LinearAlgebra.Vector[]} polygonSoup.v The vertex positions of the polygon soup mesh.
	 * @param {number[]} polygonSoup.f The indices of the polygon soup mesh.
	 * @returns {boolean} True if this mesh is constructed successfully and false if not
	 * (when this mesh contains any one or a combination of the following - non-manifold vertices,
	 *  non-manifold edges, isolated vertices, isolated faces).
	 */
	build(polygonSoup) {
		// preallocate elements
		let positions = polygonSoup["v"];
		let indices = polygonSoup["f"];
		this.preallocateElements(positions, indices);

		// create and insert vertices
		let indexToVertex = new Map();
		for (let i = 0; i < positions.length; i++) {
			let v = new Vertex();
			this.vertices[i] = v;
			indexToVertex.set(i, v);
		}

		// create and insert halfedges, edges and non boundary loop faces
		let eIndex = 0;
		let edgeCount = new Map();
		let existingHalfedges = new Map();
		let hasTwinHalfedge = new Map();
		for (let I = 0; I < indices.length; I += 3) {
			// create new face
			let f = new Face();
			this.faces[I / 3] = f;

			// create a halfedge for each edge of the newly created face
			for (let J = 0; J < 3; J++) {
				let h = new Halfedge();
				this.halfedges[I + J] = h;
			}

			// initialize the newly created halfedges
			for (let J = 0; J < 3; J++) {
				// current halfedge goes from vertex i to vertex j
				let K = (J + 1) % 3;
				let i = indices[I + J];
				let j = indices[I + K];

				// set the current halfedge's attributes
				let h = this.halfedges[I + J];
				h.next = this.halfedges[I + K];
				h.prev = this.halfedges[I + (J + 3 - 1) % 3];
				h.onBoundary = false;
				hasTwinHalfedge.set(h, false);

				// point the new halfedge and vertex i to each other
				let v = indexToVertex.get(i);
				h.vertex = v;
				v.halfedge = h;

				// point the new halfedge and face to each other
				h.face = f;
				f.halfedge = h;

				// swap if i > j
				if (i > j) j = [i, i = j][0];

				let value = [i, j]
				let key = value.toString();
				if (existingHalfedges.has(key)) {
					// if a halfedge between vertex i and j has been created in the past, then it
					// is the twin halfedge of the current halfedge
					let twin = existingHalfedges.get(key);
					h.twin = twin;
					twin.twin = h;
					h.edge = twin.edge;

					hasTwinHalfedge.set(h, true);
					hasTwinHalfedge.set(twin, true);
					edgeCount.set(key, edgeCount.get(key) + 1);

				} else {
					// create an edge and set its halfedge
					let e = new Edge();
					this.edges[eIndex++] = e;
					h.edge = e;
					e.halfedge = h;

					// record the newly created edge and halfedge from vertex i to j
					existingHalfedges.set(key, h);
					edgeCount.set(key, 1);
				}

				// check for non-manifold edges
				if (edgeCount.get(key) > 2) {
					alert("Mesh has non-manifold edges!");
					return false;
				}
			}
		}

		// create and insert boundary halfedges and "imaginary" faces for boundary cycles
		// also create and insert corners
		let hIndex = indices.length;
		let cIndex = 0;
		for (let i = 0; i < indices.length; i++) {
			// if a halfedge has no twin halfedge, create a new face and
			// link it the corresponding boundary cycle
			let h = this.halfedges[i];
			if (!hasTwinHalfedge.get(h)) {
				// create new face
				let f = new Face();
				this.boundaries.push(f);

				// walk along boundary cycle
				let boundaryCycle = [];
				let he = h;
				do {
					// create a new halfedge
					let bH = new Halfedge();
					this.halfedges[hIndex++] = bH;
					boundaryCycle.push(bH);

					// grab the next halfedge along the boundary that does not have a twin halfedge
					let nextHe = he.next;
					while (hasTwinHalfedge.get(nextHe)) {
						nextHe = nextHe.twin.next;
					}

					// set the current halfedge's attributes
					bH.vertex = nextHe.vertex;
					bH.edge = he.edge;
					bH.onBoundary = true;

					// point the new halfedge and face to each other
					bH.face = f;
					f.halfedge = bH;

					// point the new halfedge and he to each other
					bH.twin = he;
					he.twin = bH;

					// continue walk
					he = nextHe;
				} while (he !== h);

				// link the cycle of boundary halfedges together
				let n = boundaryCycle.length;
				for (let j = 0; j < n; j++) {
					boundaryCycle[j].next = boundaryCycle[(j + n - 1) % n]; // boundary halfedges are linked in clockwise order
					boundaryCycle[j].prev = boundaryCycle[(j + 1) % n];
					hasTwinHalfedge.set(boundaryCycle[j], true);
					hasTwinHalfedge.set(boundaryCycle[j].twin, true);
				}
			}

			// point the newly created corner and its halfedge to each other
			if (!h.onBoundary) {
				let c = new Corner();
				c.halfedge = h;
				h.corner = c;

				this.corners[cIndex++] = c;
			}
		}

		// check if mesh has isolated vertices, isolated faces or
		// non-manifold vertices
		if (this.hasIsolatedVertices() ||
			this.hasIsolatedFaces() ||
			this.hasNonManifoldVertices()) {
			return false;
		}

		// index elements
		this.indexElements();

		return true;
	}

	/**
	 * Preallocates mesh elements.
	 * @private
	 * @method module:Core.Mesh#preallocateElements
	 * @param {module:LinearAlgebra.Vector[]} positions The vertex positions of a polygon soup mesh.
	 * @param {number[]} indices The indices of a polygon soup mesh.
	 */
	preallocateElements(positions, indices) {
		let nBoundaryHalfedges = 0;
		let sortedEdges = new Map();
		for (let I = 0; I < indices.length; I += 3) {
			for (let J = 0; J < 3; J++) {
				let K = (J + 1) % 3;
				let i = indices[I + J];
				let j = indices[I + K];

				// swap if i > j
				if (i > j) j = [i, i = j][0];

				let value = [i, j]
				let key = value.toString();
				if (sortedEdges.has(key)) {
					nBoundaryHalfedges--;

				} else {
					sortedEdges.set(key, value);
					nBoundaryHalfedges++;
				}
			}
		}

		let nVertices = positions.length;
		let nEdges = sortedEdges.size;
		let nFaces = indices.length / 3;
		let nHalfedges = 2 * nEdges;
		let nInteriorHalfedges = nHalfedges - nBoundaryHalfedges;

		// clear arrays
		this.vertices.length = 0;
		this.edges.length = 0;
		this.faces.length = 0;
		this.halfedges.length = 0;
		this.corners.length = 0;
		this.boundaries.length = 0;
		this.generators.length = 0;

		// allocate space
		this.vertices = new Array(nVertices);
		this.edges = new Array(nEdges);
		this.faces = new Array(nFaces);
		this.halfedges = new Array(nHalfedges);
		this.corners = new Array(nInteriorHalfedges);
	}

	/**
	 * Checks whether this mesh has isolated vertices.
	 * @private
	 * @method module:Core.Mesh#hasIsolatedVertices
	 * @returns {boolean}
	 */
	hasIsolatedVertices() {
		for (let v of this.vertices) {
			if (v.isIsolated()) {
				alert("Mesh has isolated vertices!");
				return true;
			}
		}

		return false;
	}

	/**
	 * Checks whether this mesh has isolated faces.
	 * @private
	 * @method module:Core.Mesh#hasIsolatedFaces
	 * @returns {boolean}
	 */
	hasIsolatedFaces() {
		for (let f of this.faces) {
			let boundaryEdges = 0;
			for (let h of f.adjacentHalfedges()) {
				if (h.twin.onBoundary) boundaryEdges++;
			}

			if (boundaryEdges === 3) {
				alert("Mesh has isolated faces!");
				return true;
			}
		}

		return false;
	}

	/**
	 * Checks whether this mesh has non-manifold vertices.
	 * @private
	 * @method module:Core.Mesh#hasNonManifoldVertices
	 * @returns {boolean}
	 */
	hasNonManifoldVertices() {
		let adjacentFaces = new Map();
		for (let v of this.vertices) {
			adjacentFaces.set(v, 0);
		}

		for (let f of this.faces) {
			for (let v of f.adjacentVertices()) {
				adjacentFaces.set(v, adjacentFaces.get(v) + 1);
			}
		}

		for (let b of this.boundaries) {
			for (let v of b.adjacentVertices()) {
				adjacentFaces.set(v, adjacentFaces.get(v) + 1);
			}
		}

		for (let v of this.vertices) {
			if (adjacentFaces.get(v) !== v.degree()) {
				return true;
			}
		}

		return false;
	}

	/**
	 * Assigns indices to this mesh's elements.
	 * @private
	 * @method module:Core.Mesh#indexElements
	 */
	indexElements() {
		let index = 0;
		for (let v of this.vertices) {
			v.index = index++;
		}

		index = 0;
		for (let e of this.edges) {
			e.index = index++;
		}

		index = 0;
		for (let f of this.faces) {
			f.index = index++;
		}

		index = 0;
		for (let h of this.halfedges) {
			h.index = index++;
		}

		index = 0;
		for (let c of this.corners) {
			c.index = index++;
		}

		index = 0;
		for (let b of this.boundaries) {
			b.index = index++;
		}
	}
}


/**
 * Assigns an index to each element in elementList. Indices can be accessed by using
 * elements as keys in the returned dictionary.
 * @global
 * @function module:Core.indexElements
 * @param {Object[]} elementList An array of any one of the following mesh elements -
 * vertices, edges, faces, corners, halfedges, boundaries.
 * @returns {Object} A dictionary mapping each element in elementList to a unique index
 * between 0 and |elementList|-1.
 * @example
 * let vertexIndex = indexElements(mesh.vertices);
 * let v = mesh.vertices[0];
 * let i = vertexIndex[v];
 * console.log(i); // prints 0
 */
function indexElements(elementList) {
	let i = 0;
	let index = {};
	for (let element of elementList) {
		index[element] = i++;
	}

	return index;
}

// module.exports = [Mesh, indexElements]
module.exports = {
	Mesh: Mesh,
	indexElements: indexElements,
}