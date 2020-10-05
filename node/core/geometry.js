let LinearAlgebra = require('../linear-algebra/linear-algebra.js');
let Vector = LinearAlgebra.Vector;
let Complex = LinearAlgebra.Complex;
let SparseMatrix = LinearAlgebra.SparseMatrix;
let Triplet = LinearAlgebra.Triplet;
let ComplexSparseMatrix = LinearAlgebra.ComplexSparseMatrix;
let ComplexTriplet = LinearAlgebra.ComplexTriplet;

/**
 * @typedef {import('./edge.js')} Edge
 */


class Geometry {
	/**
	 * This class represents the geometry of a {@link module:Core.Mesh Mesh}. This includes information such
	 * as the position of vertices as well as methods to compute edge lengths, corner
	 * angles, face area, normals, discrete curvatures etc.
	 * @constructor module:Core.Geometry
	 * @param {module:Core.Mesh} mesh The mesh this class describes the geometry of.
	 * @param {module:LinearAlgebra.Vector[]} positions An array containing the position of each vertex in a mesh.
	 * @param {boolean} normalizePositions flag to indicate whether positions should be normalized. Default value is false.
	 * @param {boolean} normalizeEdges flag to indicate whether positions should be updated st edges of the mesh have approximately unit length. Default value is true.
	 * @property {module:Core.Mesh} mesh The mesh this class describes the geometry of.
	 * @property {Object} positions A dictionary mapping each vertex to a normalized position.
	 */
	constructor(mesh, positions, maxIndices, normalizePositions = false, normalizeEdges = true) {
		this.mesh = mesh;
		/** @type {Vector[]} */
		this.positions = {};
		this.maxIndices = maxIndices;
		this.indices = undefined;
		for (let i = 0; i < positions.length; i++) {
			let v = this.mesh.vertices[i];
			let p = positions[i];

			this.positions[v] = p;
		}

		// Build indices list (for rendering)
		let F = this.mesh.faces.length;
		let indices = new Uint32Array(F * 3);
		for (let f of mesh.faces) {
			let i = 0;
			for (let v of f.adjacentVertices()) {
				indices[3 * f.index + i++] = v.index;
			}
		};
		this.indices = indices;

		if (normalizePositions) {
			normalize(this.positions, mesh.vertices);
		}

		if (normalizeEdges) {
			// Center mesh on origin
			normalize(this.positions, mesh.vertices, false);

			// Compute avg edge length
			let edgeLength = 0;
			let n = 0;
			for (let e of mesh.edges) {
				edgeLength += this.length(e);
				n++;
			}
			edgeLength /= n;

			// Rescale mesh to have ~ unit length edges
			for (let v of mesh.vertices) {
				let p = this.positions[v];
				p.divideBy(edgeLength);
			}
		}
	}

	/**
	 * Check whether an edge can be flipped and whether the flip miminizes the max angle
	 * @param {Edge} edge 
	 */
	isFlippable(edge) {
		// Check boundary
		if (edge.onBoundary())
			return false;
		let he = edge.halfedge;
		let twin = edge.halfedge.twin;
		let x2 = he.prev.vertex;
		let x3 = twin.prev.vertex;

		// Check angles
		let alpha0a = this.angle(he.next.corner);
		let alpha0b = this.angle(twin.prev.corner);
		let alpha1a = this.angle(he.prev.corner);
		let alpha1b = this.angle(twin.next.corner);
		if (alpha0a + alpha0b >= Math.PI || alpha1a + alpha1b >= Math.PI)
			return false;
		
		// Check dihedral angle
		if (Math.abs(this.dihedralAngle(he)) > 0.3)
			return false;

		// Check that flip maximizes min angle
		let alpha2 = this.angle(he.corner);
		let alpha3 = this.angle(twin.corner);
		let oldMaxAngle = Math.max(alpha0a, alpha0b, alpha1a, alpha1b, alpha2, alpha3);

		// New angles in case of flip
		let alpha0 = alpha0a + alpha0b;
		let alpha1 = alpha1a + alpha1b;
		let v23 = this.positions[x3.index].minus(this.positions[x2.index]).unit();
		let alpha2a = v23.angle(this.vector(he.prev));
		let alpha2b = v23.angle(this.vector(he.next).negated());
		let alpha3a = v23.negated().angle(this.vector(twin.next).negated());
		let alpha3b = v23.negated().angle(this.vector(twin.prev));
		let newMaxAngle = Math.max(alpha0, alpha1, alpha2a, alpha2b, alpha3a, alpha3b);

		return newMaxAngle < oldMaxAngle;
	}

	/**
	 * Flips an edge and update mesh connectivity and indices list
	 * @param {Edge} edge 
	 */
	flip(edge) {
		let he = edge.halfedge;
		let twin = edge.halfedge.twin;

		let nhe = he.next;
		let phe = he.prev;
		let ntwin = twin.next;
		let ptwin = twin.prev;

		let fhe = he.face;
		let ftwin = twin.face;

		let x0 = he.vertex;
		let x1 = twin.vertex;
		let x2 = phe.vertex;
		let x3 = ptwin.vertex;

		// Connectivity updates
		x0.halfedge = ntwin;
		x1.halfedge = nhe;

		ntwin.face = fhe;
		nhe.face = ftwin;

		fhe.halfedge = he;
		ftwin.halfedge = twin;

		he.vertex = x3;
		twin.vertex = x2;

		this.mesh.linkHalfedges(ntwin, he);
		this.mesh.linkHalfedges(he, phe);
		this.mesh.linkHalfedges(phe, ntwin);

		this.mesh.linkHalfedges(nhe, twin);
		this.mesh.linkHalfedges(twin, ptwin);
		this.mesh.linkHalfedges(ptwin, nhe);

		// Update indices list
		let newIndices = [...this.indices];
		newIndices.splice(3 * fhe.index, 3, ...[x0.index, x3.index, x2.index]);
		newIndices.splice(3 * ftwin.index, 3, ...[x1.index, x2.index, x3.index]);
		this.indices = newIndices;

	}

	/**
	 * Splits an edge, add new elements to mesh, update mesh connectivity and update indices list
	 * @param {Edge} edge 
	 * @return {number} The index of the new vertex
	 */
	split(edge) {
		if (this.mesh.vertices.length >= this.maxPoints - 1) {
			console.error("Max number of vertices reached");
			return null;
		}

		const onBoundary = edge.onBoundary();

		// Halfedge and twin of the original edge
		const ohe = edge.halfedge;
		const otwin = edge.halfedge.twin;

		// Create new vertex in the middle of the edge
		let vNewPos = this.midpoint(edge);
		let vNew = this.mesh.newVertex();
		this.positions[vNew.index] = vNewPos;

		// Create a new edge
		let e0New = this.mesh.newEdge();

		// Create halfedges
		let he01 = this.mesh.newHalfedge();
		let he02 = this.mesh.newHalfedge();

		// he01 is halfedge of new edge
		e0New.halfedge = he01;
		he01.edge = e0New;

		// he02 is halfedge of old edge
		he02.edge = edge;
		
		// Associate new vertex with new halfedges
		he01.vertex = vNew;
		he02.vertex = vNew;
		vNew.halfedge = he01;

		// Create a new face for each new halfedge not on a boundary
		let newFaces = [];
		let newHalfedges;
		let originalHalfedges;
		if (!onBoundary) {
			newHalfedges = [he01, he02];
			originalHalfedges = [ohe, otwin];
		}
		else {
			newHalfedges = [he01];
			originalHalfedges = [ohe];

			// he02 is on boundary
			he02.face = this.mesh.boundaries[0];
			he02.onBoundary = true;
			// Relink boundary halfedge loop
			this.mesh.linkHalfedges(he02, otwin.next);
			this.mesh.linkHalfedges(otwin, he02);
			
		}
		for (let i = 0; i < newHalfedges.length; i++) {
			let he = originalHalfedges[i];
			let nhe = newHalfedges[i];

			// Create a new face
			let nface = this.mesh.newFace();
			newFaces.push(nface);

			// Link new face with new halfedge
			nface.halfedge = nhe;
			nhe.face = nface;
			// Link old face with old halfedge
			he.face.halfedge = he;

			// Set face for halfedge opposite new vertex
			he.next.face = nface;

			// Create a new edge
			let ne = this.mesh.newEdge();

			// Create halfedges around new edge
			let hea = this.mesh.newHalfedge();
			let heb = this.mesh.newHalfedge();
			// Set edge
			ne.halfedge = hea;
			hea.edge = ne;
			heb.edge = ne;
			// Set twins
			this.mesh.glueHalfedges(hea, heb);
			// Set face
			hea.face = nface;
			heb.face = he.face;
			// Set vertex
			hea.vertex = he.prev.vertex;
			heb.vertex = vNew;

			// Link halfedges together
			this.mesh.linkHalfedges(nhe, he.next);
			this.mesh.linkHalfedges(hea, nhe);
			this.mesh.linkHalfedges(nhe.next, hea);
			this.mesh.linkHalfedges(heb, he.prev);
			this.mesh.linkHalfedges(he, heb);

			// Create corners for the halfedges
			let interiorHalfedges = [nhe, hea, heb];
			for (let k = 0; k < interiorHalfedges.length; k++) {
				let hek = interiorHalfedges[k];
				let c = this.mesh.newCorner();
				c.halfedge = hek;
				hek.corner = c;
				hek.onBoundary = false;
			}

		}

		// Glue old and new halfedges together
		this.mesh.glueHalfedges(he01, otwin);
		this.mesh.glueHalfedges(he02, ohe);
		
		// Update indices list
		let newIndices = [...this.indices];
		let oldFaces;
		if (!onBoundary) {
			oldFaces = [ohe.face, otwin.face];
		}
		else {
			oldFaces = [ohe.face];
		}

		// Replace old faces indices
		for (let oldFace of oldFaces) {
			// Get face indices
			let indices = [];
			for (let v of oldFace.adjacentVertices()) {
				indices.push(v.index);
			}
			// Replace
			newIndices.splice(3 * oldFace.index, 3, ...indices);
		}
		// Add new faces indices
		for(let newFace of newFaces) {
			// Get face indices
			let indices = [];
			for (let v of newFace.adjacentVertices()) {
				indices.push(v.index);
			}
			// Add
			newIndices.splice(3 * newFace.index, 0, ...indices);
		}
		this.indices = newIndices;

		return vNew.index;

	}

	/**
	 * Computes the vector along a halfedge.
	 * @method module:Core.Geometry#vector
	 * @param {module:Core.Halfedge} h The halfedge along which the vector needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vector(h) {
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];

		return b.minus(a);
	}

	/**
	 * Computes the length of an edge.
	 * @method module:Core.Geometry#length
	 * @param {module:Core.Edge} e The edge whose length needs to be computed.
	 * @returns {number}
	 */
	length(e) {
		return this.vector(e.halfedge).norm();
	}

	/**
	 * Computes the midpoint of an edge.
	 * @method module:Core.Geometry#midpoint
	 * @param {module:Core.Edge} e The edge whose midpoint needs to be computed.
	 * @returns {number}
	 */
	midpoint(e) {
		let h = e.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.twin.vertex];

		return (a.plus(b)).over(2);
	}

	/**
	 * Computes the mean edge length of all the edges in a mesh.
	 * @method module:Core.Geometry#meanEdgeLength
	 * @returns {number}
	 */
	meanEdgeLength() {
		let sum = 0;
		let edges = this.mesh.edges;
		for (let e of edges) {
			sum += this.length(e);
		}

		return sum / edges.length;
	}

	/**
	 * Computes the area of a face.
	 * @method module:Core.Geometry#area
	 * @param {module:Core.Face} f The face whose area needs to be computed.
	 * @returns {number}
	 */
	area(f) {
		if (f.isBoundaryLoop()) return 0.0;

		let u = this.vector(f.halfedge);
		let v = this.vector(f.halfedge.prev).negated();

		return 0.5 * u.cross(v).norm();
	}

	/**
	 * Computes the total surface area of a mesh.
	 * @method module:Core.Geometry#totalArea
	 * @returns {number}
	 */
	totalArea() {
		let sum = 0.0;
		for (let f of this.mesh.faces) {
			sum += this.area(f);
		}

		return sum;
	}

	/**
	 * Computes the normal of a face.
	 * @method module:Core.Geometry#faceNormal
	 * @param {module:Core.Face} f The face whose normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	faceNormal(f) {
		if (f.isBoundaryLoop()) return undefined;

		let u = this.vector(f.halfedge);
		let v = this.vector(f.halfedge.prev).negated();

		return u.cross(v).unit();
	}

	/**
	 * Computes the centroid of a face.
	 * @method module:Core.Geometry#centroid
	 * @param {module:Core.Face} f The face whose centroid needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	centroid(f) {
		let h = f.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];
		let c = this.positions[h.prev.vertex];

		if (f.isBoundaryLoop()) return a.plus(b).over(2);

		return a.plus(b).plus(c).over(3);
	}

	/**
	 * Computes the circumcenter of a face.
	 * @method module:Core.Geometry#circumcenter
	 * @param {module:Core.Face} f The face whose circumcenter needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	circumcenter(f) {
		let h = f.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];
		let c = this.positions[h.prev.vertex];

		if (f.isBoundaryLoop()) return a.plus(b).over(2);

		let ac = c.minus(a);
		let ab = b.minus(a);
		let w = ab.cross(ac);

		let u = (w.cross(ab)).times(ac.norm2());
		let v = (ac.cross(w)).times(ab.norm2());
		let x = (u.plus(v)).over(2 * w.norm2());

		return x.plus(a);
	}

	/**
	 * Computes an orthonormal bases for a face.
	 * @method module:Core.Geometry#orthonormalBases
	 * @param {module:Core.Face} f The face on which the orthonormal bases needs to be computed.
	 * @returns {module:LinearAlgebra.Vector[]} An array containing two orthonormal vectors tangent to the face.
	 */
	orthonormalBases(f) {
		let e1 = this.vector(f.halfedge).unit();

		let normal = this.faceNormal(f);
		let e2 = normal.cross(e1);

		return [e1, e2];
	}

	/**
	 * Computes the angle (in radians) at a corner.
	 * @method module:Core.Geometry#angle
	 * @param {module:Core.Corner} c The corner at which the angle needs to be computed.
	 * @returns {number} The angle clamped between 0 and π.
	 */
	angle(c) {
		let u = this.vector(c.halfedge.prev).unit();
		let v = this.vector(c.halfedge.next).negated().unit();

		return Math.acos(Math.max(-1.0, Math.min(1.0, u.dot(v))));
	}

	/**
	 * Computes the cotangent of the angle opposite to a halfedge.
	 * @method module:Core.Geometry#cotan
	 * @param {module:Core.Halfedge} h The halfedge opposite to the angle whose cotangent needs to be computed.
	 * @returns {number}
	 */
	cotan(h) {
		if (h.onBoundary) return 0.0;

		let u = this.vector(h.prev);
		let v = this.vector(h.next).negated();

		return u.dot(v) / u.cross(v).norm();
	}

	/**
	 * Computes the signed angle (in radians) between two adjacent faces.
	 * @method module:Core.Geometry#dihedralAngle
	 * @param {module:Core.Halfedge} h The halfedge (shared by the two adjacent faces) on which
	 * the dihedral angle is computed.
	 * @returns {number} The dihedral angle.
	 */
	dihedralAngle(h) {
		if (h.onBoundary || h.twin.onBoundary) return 0.0;

		let n1 = this.faceNormal(h.face);
		let n2 = this.faceNormal(h.twin.face);
		let w = this.vector(h).unit();

		let cosTheta = n1.dot(n2);
		let sinTheta = n1.cross(n2).dot(w);

		return Math.atan2(sinTheta, cosTheta);
	}

	/**
	 * Computes the barycentric dual area of a vertex.
	 * @method module:Core.Geometry#barycentricDualArea
	 * @param {module:Core.Vertex} v The vertex whose barycentric dual area needs to be computed.
	 * @returns {number}
	 */
	barycentricDualArea(v) {
		let area = 0.0;
		for (let f of v.adjacentFaces()) {
			area += this.area(f) / 3;
		}

		return area;
	}

	/**
	 * Computes the circumcentric dual area of a vertex.
	 * @see {@link http://www.cs.cmu.edu/~kmcrane/Projects/Other/TriangleAreasCheatSheet.pdf}
	 * @method module:Core.Geometry#circumcentricDualArea
	 * @param {module:Core.Vertex} v The vertex whose circumcentric dual area needs to be computed.
	 * @returns {number}
	 */
	circumcentricDualArea(v) {
		let area = 0.0;
		for (let h of v.adjacentHalfedges()) {
			let u2 = this.vector(h.prev).norm2();
			let v2 = this.vector(h).norm2();
			let cotAlpha = this.cotan(h.prev);
			let cotBeta = this.cotan(h);

			area += (u2 * cotAlpha + v2 * cotBeta) / 8;
		}

		return area;
	}

	/**
	 * Computes the normal at a vertex using the "equally weighted" method.
	 * @method module:Core.Geometry#vertexNormalEquallyWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalEquallyWeighted(v) {
		let n = new Vector();
		for (let f of v.adjacentFaces()) {
			let normal = this.faceNormal(f);

			n.incrementBy(normal);
		}

		n.normalize();

		return n;
	}

	/**
	 * Computes the normal at a vertex using the "face area weights" method.
	 * @method module:Core.Geometry#vertexNormalAreaWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalAreaWeighted(v) {
		let n = new Vector();
		for (let f of v.adjacentFaces()) {
			let normal = this.faceNormal(f);
			let area = this.area(f);

			n.incrementBy(normal.times(area));
		}

		n.normalize();

		return n;
	}

	/**
	 * Computes the normal at a vertex using the "tip angle weights" method.
	 * @method module:Core.Geometry#vertexNormalAngleWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalAngleWeighted(v) {
		let n = new Vector();
		for (let c of v.adjacentCorners()) {
			let normal = this.faceNormal(c.halfedge.face);
			let angle = this.angle(c);

			n.incrementBy(normal.times(angle));
		}

		n.normalize();

		return n;
	}

	/**
	 * Computes the normal at a vertex using the "gauss curvature" method.
	 * @method module:Core.Geometry#vertexNormalGaussCurvature
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalGaussCurvature(v) {
		let n = new Vector();
		for (let h of v.adjacentHalfedges()) {
			let weight = 0.5 * this.dihedralAngle(h) / this.length(h.edge);

			n.decrementBy(this.vector(h).times(weight));
		}

		n.normalize();

		return n;
	}

	/**
	 * Computes the normal at a vertex using the "mean curvature" method (same as the "area gradient" method).
	 * @method module:Core.Geometry#vertexNormalMeanCurvature
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalMeanCurvature(v) {
		let n = new Vector();
		for (let h of v.adjacentHalfedges()) {
			let weight = 0.5 * (this.cotan(h) + this.cotan(h.twin));

			n.decrementBy(this.vector(h).times(weight));
		}

		n.normalize();

		return n;
	}

	/**
	 * Computes the normal at a vertex using the "inscribed sphere" method.
	 * @method module:Core.Geometry#vertexNormalSphereInscribed
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalSphereInscribed(v) {
		let n = new Vector();
		for (let c of v.adjacentCorners()) {
			let u = this.vector(c.halfedge.prev);
			let v = this.vector(c.halfedge.next).negated();

			n.incrementBy(u.cross(v).over(u.norm2() * v.norm2()));
		}

		n.normalize();

		return n;
	}

	/**
	 * Computes the angle defect at a vertex (= 2π minus the sum of incident angles
	 * at an interior vertex or π minus the sum of incident angles at a boundary vertex).
	 * @method module:Core.Geometry#angleDefect
	 * @param {module:Core.Vertex} v The vertex whose angle defect needs to be computed.
	 * @returns {number}
	 */
	angleDefect(v) {
		let angleSum = 0.0;
		for (let c of v.adjacentCorners()) {
			angleSum += this.angle(c);
		}

		return v.onBoundary() ? Math.PI - angleSum : 2 * Math.PI - angleSum;
	}

	/**
	 * Computes the (integrated) scalar gauss curvature at a vertex.
	 * @method module:Core.Geometry#scalarGaussCurvature
	 * @param {module:Core.Vertex} v The vertex whose gauss curvature needs to be computed.
	 * @returns {number}
	 */
	scalarGaussCurvature(v) {
		return this.angleDefect(v);
	}

	/**
	 * Computes the (integrated) scalar mean curvature at a vertex.
	 * @method module:Core.Geometry#scalarMeanCurvature
	 * @param {module:Core.Vertex} v The vertex whose mean curvature needs to be computed.
	 * @returns {number}
	 */
	scalarMeanCurvature(v) {
		let sum = 0.0;
		for (let h of v.adjacentHalfedges()) {
			sum += 0.5 * this.length(h.edge) * this.dihedralAngle(h);
		}

		return sum;
	}

	/**
	 * Computes the total angle defect (= 2π times the euler characteristic of the mesh).
	 * @method module:Core.Geometry#totalAngleDefect
	 * @returns {number}
	 */
	totalAngleDefect() {
		let totalDefect = 0.0;
		for (let v of this.mesh.vertices) {
			totalDefect += this.angleDefect(v);
		}

		return totalDefect;
	}

	/**
	 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
	 * @method module:Core.Geometry#principalCurvatures
	 * @param {module:Core.Vertex} v The vertex on which the principal curvatures need to be computed.
	 * @returns {number[]} An array containing the minimum and maximum principal curvature values at a vertex.
	 */
	principalCurvatures(v) {
		let A = this.circumcentricDualArea(v);
		let H = this.scalarMeanCurvature(v) / A;
		let K = this.angleDefect(v) / A;

		let discriminant = H * H - K;
		if (discriminant > 0) discriminant = Math.sqrt(discriminant);
		else discriminant = 0;

		let k1 = H - discriminant;
		let k2 = H + discriminant;

		return [k1, k2];
	}

	/**
	 * Builds a sparse laplace matrix. The laplace operator is negative semidefinite;
	 * instead we build a positive definite matrix by multiplying the entries of the
	 * laplace matrix by -1 and shifting the diagonal elements by a small constant (e.g. 1e-8).
	 * @method module:Core.Geometry#laplaceMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	laplaceMatrix(vertexIndex) {
		let V = this.mesh.vertices.length;
		let T = new Triplet(V, V);
		for (let v of this.mesh.vertices) {
			let i = vertexIndex[v];
			let sum = 1e-8;

			for (let h of v.adjacentHalfedges()) {
				let j = vertexIndex[h.twin.vertex];
				let weight = (this.cotan(h) + this.cotan(h.twin)) / 2;
				sum += weight;

				T.addEntry(-weight, i, j);
			}

			T.addEntry(sum, i, i);
		}

		return SparseMatrix.fromTriplet(T);
	}

	/**
	 * Builds a sparse diagonal mass matrix containing the barycentric dual area of each vertex
	 * of a mesh.
	 * @method module:Core.Geometry#massMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	massMatrix(vertexIndex) {
		let V = this.mesh.vertices.length;
		let T = new Triplet(V, V);
		for (let v of this.mesh.vertices) {
			let i = vertexIndex[v];

			T.addEntry(this.barycentricDualArea(v), i, i);
		}

		return SparseMatrix.fromTriplet(T);
	}

	/**
	 * Builds a sparse complex laplace matrix. The laplace operator is negative semidefinite;
	 * instead we build a positive definite matrix by multiplying the entries of the
	 * laplace matrix by -1 and shifting the diagonal elements by a small constant (e.g. 1e-8).
	 * @method module:Core.Geometry#complexLaplaceMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.ComplexSparseMatrix}
	 */
	complexLaplaceMatrix(vertexIndex) {
		let V = this.mesh.vertices.length;
		let T = new ComplexTriplet(V, V);
		for (let v of this.mesh.vertices) {
			let i = vertexIndex[v];
			let sum = 1e-8;

			for (let h of v.adjacentHalfedges()) {
				let j = vertexIndex[h.twin.vertex];
				let weight = (this.cotan(h) + this.cotan(h.twin)) / 2;
				sum += weight;

				T.addEntry(new Complex(-weight), i, j);
			}

			T.addEntry(new Complex(sum), i, i);
		}

		return ComplexSparseMatrix.fromTriplet(T);
	}
}

/**
 * Centers a mesh about the origin and rescales it to unit radius.
 * @global
 * @function module:Core.normalize
 * @param {module:LinearAlgebra.Vector[]} positions The position of each vertex in the vertices array.
 * @param {module:Core.Vertex[]} vertices The vertices of a mesh.
 * @param {boolean} rescale A flag indicating whether mesh positions should be scaled to a unit radius.
 */
function normalize(positions, vertices, rescale = true) {
	// compute center of mass
	let N = vertices.length;
	let cm = new Vector();
	for (let v of vertices) {
		let p = positions[v];

		cm.incrementBy(p);
	}
	cm.divideBy(N);

	// translate to origin and determine radius
	let radius = -1;
	for (let v of vertices) {
		let p = positions[v];

		p.decrementBy(cm);
		radius = Math.max(radius, p.norm());
	}

	// rescale to unit radius
	if (rescale) {
		for (let v of vertices) {
			let p = positions[v];

			p.divideBy(radius);
		}
	}
}

// module.exports = [Geometry, normalize]
module.exports = {
	Geometry: Geometry,
	normalize: normalize,
}