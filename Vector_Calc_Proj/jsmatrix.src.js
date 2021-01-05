/*
*
*             JSMatrix : JavaScript matrix and vector implementation          
*                                                                             
*                            version 0.8 (2014-08-20)
*  
*                    Copyright (C)  2011 - 2014  Jan Stransky                 
*  
*           Czech Technical University, Faculty of Civil Engineering,         
*       Department of Structural Mechanics, 166 29 Prague, Czech Republic     
*  
*  JSMatrix is free software: you can redistribute it and/or modify it under
*  the terms of the GNU Lesser General Public License as published by the Free
*  Software Foundation, either version 3 of the License, or (at your option)
*  any later version.
*  
*  JSMatrix is distributed in the hope that it will be useful, but WITHOUT ANY
*  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
*  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
*  more details.
*  
*  You should have received a copy of the GNU Lesser General Public License
*  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


/**
* @fileOverview JavaScript vector and matrix implementation.
* This file contains JavaScript implementation of vectors and matrices, i.e 1d and 2d arrays of real (floating point) numbers and corresponding mathematical operations. The objects creation is inspired by Numerical Python (numpy) package ( v = JSMatrix.vec([1,2,3]), m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]) ). The implementation is such that coefficient access can be done directly by a=v[1], v[1]=3., a=m[1][0], m[1][0] = 4. (it is simple and easy, but no dimension is checked), or (in safer but longer way) by get/set method, where dimensions are checked ( a=v.get(1), v.set(1,3.), a=m.get(1,0), m.set(1,0,4.) ). Access function numbering is 0-based(!) (first element has index 0).
* <br /><br />
* When a method returns a new object (e.g. ret = m1.add(m2); ), if operation is successfull (e.g. dimensions agree), the object is returned, if operation fails (e.g. dimensions are mismatched), null value is returned (inspired by <a href="http://sylvester.jcoglan.com/">Sylvester project</a>).
* <br /><br />
* In vector-matrix context, vectors are considered in their natural representation (row x column vectors). For example in matrix*vector multiplication the vector is considered as 1d column matrix, while in vector*matrix the vector is considered as 1d row matrix (similar approach is used in Mathematica program). In documantation, vectors are considered as 1d column matrices (e.g. vec^T indiceates transposition of vec, i.e. 1d row matrix).
* <br /><br />
* Method naming is usualy following (example is on matrix transposition): verb infinitive (e.g. transpose) is used for function changing receiver - receiver is changed (to its transposition for example) and nothing is returned. Verb past simple (e.g. transposed) is used for function returning new object (transposed copy of receiver for example), receiver is not changed. Names beSomethingOf (e.g. beTranspositionOf) indicates that receiver is modified to become something (transposition of given matrix for example), receiver is changed and given parameter is not. Names isSomething or canSomething (e.g. isTranspositionOf) returns usually bools (true if receiver is transposition of given matrix for example).
* <br /><br />
* Current implementation contains: basic operations (coefficient access, +, -, *, /, +=, -=, *=, vector norms), testing operations, subvector/submatrix operations, solving linear system of equations, eigenvalues and eigenvector solving, matrix decmpositions (LU, LDU, Cholesky, QR, SVD, Schur, spectral, polar) and many others.
* <br /><br />
* For more information see <a href="http://mech.fsv.cvut.cz/~stransky/software/jsmatrix/">project homepage</a> or further documentation.
* <br /><br />
* JSMatrix is a free software distributed under <a href='http://www.gnu.org/licenses/lgpl.html'>GNU LGPL license</a>.
* <br /><br />
* Example:<br />
* v = new JSMatrix.Vector(); // creates empty vector object<br />
* v = new JSMatrix.Vector(4); // creates (zeroed) vector with 4 elements<br />
* v = JSMatrix.Vector.create([2,4,3]); // creates vector corresponding to given array (length 3, elements=2,4,3)<br />
* a = v[1]; // a = 4.<br />
* a[1] = 5.; // sets 5. to a[1]<br />
* a[9] = 3.; // sets 3. to a[9] (out of bounds, no error!!)<br />
* a.set(9,3.); // sets 3. to a[9], checks bounds, finds that index 9 is out of bounds and does nothing<br />
* <br />
* References:<br />
* [1] Quarteroni, A., Sacco, R. and Saleri, F. 2010. Numerical Mathematics, 2nd edition. Springer-Verlag New York<br />
* [2] http://en.wikipedia.org/wiki/Matrix_decomposition<br />
* [3] Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. 2007. The Numerical Recipes: The Art of Scientific Computing, Third Edition<br />
* [4] inspiration from <a href="http://www.oofem.org>OOFEM project</a><br />
* [5] inspiration from <a href="http://sylvester.jcoglan.com/">Sylvester project</a><br />
* @author Jan Stránský <http://mech.fsv.cvut.cz/~stransky>
* @version 0.8 (2014-08-20)
*/








/** Global namespace object
* @namespace
*/
var JSMatrix = {
	/** ( = 0 )
	* @constant
	*/
	X: 0,
	/** ( = 1 )
	* @constant
	*/
	Y: 1,
	/** ( = 2 )
	* @constant
	*/
	Z: 2,
	/** ( = 1e-8)
	* @constant
	*/
	TOL: 1e-8,
	/** ( = 1e-32)
	* @constant
	*/
	TOL2: 1e-32
};




/** Vector implementation
* @class Represents a vector (1D matrix) of real (floating point) numbers
* @param {number=} nElements (default=0) number of elements of newly created Vector object
* @property {number} nElements Number of vector elements
* @this {JSMatrix.Vector}
* @constructor
*/
JSMatrix.Vector = function(nElements) {
	this.nElements = nElements || 0;
	for (var e=0; e<this.nElements; e++) { this[e] = 0.; }
};

/** Matrix implementation
* @class Represents 2D matrix of real (floating point) numbers
* @param {number=} nRows (default=0) number or rows of newly created Matrix
* @param {number=} nCols (default=0) number or columns of newly created Matrix
* @property {number} nRows number of rows of receiver
* @property {number} nCols number of columns of receiver
* @this {JSMatrix.Matrix}
* @constructor
*/
JSMatrix.Matrix = function(nRows,nCols) {
	this.nRows = nRows || 0;
	this.nCols = nCols || 0;
	var r,c;
	for (r=0; r<this.nRows; r++) {
		this[r] = {};
		for (c=0; c<this.nCols; c++) { this[r][c] = 0.; }
	}
};





/** Empties receiver (resizes it to size 0)
* @example
* v = JSMatrix.vec([1,2,3,4]);
* v.empty(); // v = Vector([])
*/
JSMatrix.Vector.prototype.empty = function() {
	this.resize(0);
};

/** Tests if receiver is empty (has size 0)
* @returns true if receiver is empty, false otherwise
* @example
* v = new JSMatrix.Vector(3);
* b1 = v.isEmpty(); // b1 = false
* v.empty();
* b2 = v.isEmpty(); // b2 = true
*/
JSMatrix.Vector.prototype.isEmpty = function() {
	return this.nElements==0;
};

/** Zeroes all elements of receiver
* @example
* v = JSMatrix.vec([1,2,3,4]);
* v.zero(); // v = Vector([0,0,0,0])
*/
JSMatrix.Vector.prototype.zero = function() {
	for (var e=0; e<this.nElements; e++) { this[e] = 0.; }
};

/** Sets all elements of receiver to 1.0
* @example
* v = JSMatrix.vec([1,2,3,4]);
* v.beFullOfOnes(); // v = Vector([1,1,1,1])
*/
JSMatrix.Vector.prototype.beFullOfOnes = function() {
	for (var e=0; e<this.nElements; e++) { this[e] = 1.; }
};

/** Alias for beFullOfOnes(), see {@link JSMatrix.Vector#beFullOfOnes}
* @function
*/
JSMatrix.Vector.prototype.one = JSMatrix.Vector.prototype.beFullOfOnes;

/** Returns size (numer of elements) of receiver
* @returns {number} length (number of elements) of receiver
* @example
* v = JSMatrix.vec([1,5,8,2]);
* s = v.size(); // s = 4
*/
JSMatrix.Vector.prototype.size = function() {
	return this.nElements;
};

/** Returns sum of all elements of receiver
* @returns {number} sum of all elements of receiver
* @example
* v = JSMatrix.vec([1,2,3]);
* s = v.sum(); // s = 6
*/
JSMatrix.Vector.prototype.sum = function() {
	var ret = 0.;
	for (var e=0; e<this.nElements; e++) { ret += this[e]; }
	return ret;
};

/** Returns maximal value of receiver elements
* @returns {number} maximal value of receiver elements
* @example
* v = JSMatrix.vec([1,-4,3,2]);
* m = v.max(); // m = 3
*/
JSMatrix.Vector.prototype.max = function() {
	var ret = -1e64;
	for (var e=0; e<this.nElements; e++) {
		if (this[e] > ret) { ret = this[e]; }
	}
	return ret;
};

/** Returns minimal value of receiver elements
* @returns {number} minimal value of receiver elements
* @example
* v = JSMatrix.vec([1,-4,3,2]);
* m = v.min(); // m = -4
*/
JSMatrix.Vector.prototype.min = function() {
	var ret = 1e64;
	for (var e=0; e<this.nElements; e++) {
		if (this[e] < ret) { ret = this[e]; }
	}
	return ret;
};

/** Returns Vector with elementwise abs value of receiver elements ( ret[i] = abs(this[i]) )
* @returns {JSMatrix.Vector} elementwise absolute value of receiver
* @example
* v = JSMatrix.vec([1,-4,3,-2]);
* m = v.abs(); // m = Vector([1,4,3,2])
*/
JSMatrix.Vector.prototype.abs = function() {
	var ret = new JSMatrix.Vector(this.nElements);
	for (var e=0; e<this.nElements; e++) {
		ret[e] = Math.abs(this[e]);
	}
	return ret;
};

/** Returns index of first occurrence of given value
* @param {number} val value whose index will be returned
* @returns {number|null} index of first occurrence of val
* @example
* v = JSMatrix.vec([5,6,1,2,6,8,6,2]);
* i1 = v.indexOf(2); // i1 = 3
* i2 = v.indexOf(6); // i2 = 1
* i3 = v.indexOf(8); // i3 = 5
*/
JSMatrix.Vector.prototype.indexOf = function(val) {
	for (var e=0; e<this.nElements; e++) {
		if (Math.abs(this[e]-val) < JSMatrix.TOL) { return e; }
	}
	return null;
};

/** Testing vectors equality
* @param {JSMatrix.Vector} vec vector to be compared with receiver
* @returns {boolean} true if this==vec, false otherwise
* @example
* v1 = JSMatrix.vec([1,2,3]);
* v2 = JSMatrix.vec([1,2,2.999999999]);
* v3 = JSMatrix.vec([1,2,1]);
* t1 = v1.isEqualTo(v2); // t1 = true
* t2 = v1.isEqualTo(v3); // t2 = false
*/
JSMatrix.Vector.prototype.isEqualTo = function(vec) {
	if (!this.isSameSizeAs(vec)) { return false; }
	for (var e=0; e<this.nElements; e++) {
		if (Math.abs(this[e]-vec[e]) > JSMatrix.TOL) { return false; }
	}
	return true;
};

/** Coefficient access function, returns one element of receiver ( ret = this[r] ). Vector.get(e) is almost equivalent to direct Vector[e] access method
* @param {number} e index of element to be returned
* @returns {number|null} value of the e-th element
* @example
* v = JSMatrix.vec([4,7,2,4]);
* a = v.get(1); // a = 7
* a = v[1]; // a = 7
* a = v.get(9); // a = null
* a = v[9]; // a = undefined
*/
JSMatrix.Vector.prototype.get = function(e) {
	if (e<0 || e>=this.nElements) { /*jsmLogger.warn('Vector.get: out fo bounds');*/ return null; }
	return this[e];
};

/** Coefficient access function, sets one element of receiver ( this[r] = val ). Vector.set(e) is almost equivalent to direct Vector[e] access method
* @param {number} e index of element to be set
* @param {number} val value to be set
* @example
* v = JSMatrix.vec([4,7,2,4]);
* v.set(1,68.); // v = Vector([4,68,2,4])
* v[1] = 43.; // v = Vector([4,43,2,4])
* v.set(9,68.); // nothing is done - out of bounds
* // v = Vector([4,43,2,4])
* v[9] = 43.;
* /// new object this[9] is created, but no reflection to actual Vector is made
* // v = Vector([4,43,2,4])
*/
JSMatrix.Vector.prototype.set = function(e,val) {
	if (e<0 || e>=this.nElements) { /*jsmLogger.warn('Vector.set: out of bounds');*/ return; }
	this[e] = val;
};

/** Coefficient access function, increment one element of receiver ( this[r] += val ). Vector.incr(e,val) is much safer than direct Vector[e] += val (see example)
* @param {number} e index of element to be incremented
* @param {number} val value to be added
* @example
* v = JSMatrix.vec([4,7,2,4]);
* v.incr(1,4.); // v = Vector([4,11,2,4])
* v[1] += 4.; // v = Vector([4,15,2,4])
* v.incr(9,4.); // nothing is done, out of bounds
* // v = Vector([4,15,2,4])
* /// v[9] += 4. would result into error (v[9] is not defined)
*/
JSMatrix.Vector.prototype.incr = function(e,val) {
	if (e<0 || e>=this.nElements) { /*jsmLogger.warn('Vector.incr: out of bounds');*/ return; }
	if (val) { this[e] += val; }
};

/** Testing vectors size equality
* @param {JSMatrix.Vector} vec vector to be tested
* @returns {boolean} true if vec and receiver has same size, false otherwise
* @example
* a = JSMatrix.vec([1,2,3,4]);
* b = JSMatrix.vec([5,6,7,8,9,10,11,12]);
* c = JSMatrix.vec([13,14,15,16]);
* t1 = a.isSameSizeAs(b); // t1 = false
* t2 = a.isSameSizeAs(c); // t2 = true
*/
JSMatrix.Vector.prototype.isSameSizeAs = function(vec) {
	return this.nElements==vec.nElements;
};

/** Testing if receiver can multiply given matrix from left ( this^T * mat )
* @param {JSMatrix.Matrix} mat matrix to be tested
* @returns {boolean} true if the multiplication is possible, false otherwise
* @example
* v = JSMatrix.vec([1,2,3]);
* m1 = JSMatrix.mat([[11,12,13],[21,22,23]]);
* m2 = JSMatrix.mat([[11,12],[21,22],[31,32]]);
* t1 = v.canMultiplyMat(m1); // t1 = true
* t2 = v.canMultiplyMat(m2); // t2 = false
*/
JSMatrix.Vector.prototype.canMultiplyMat = function(mat) {
return this.nElements ==  mat.nCols;
};

/** Alias for canMultiplyMat(), see {@link JSMatrix.Vector#canMultiplyMat}
* @function
*/
JSMatrix.Vector.canMultiplyMatFromLeft = JSMatrix.Vector.prototype.canMultiplyMat;

/** Sets elements of receiver form given Array object
* @param {Array.<number>} arry Array object containing new receiver elements
* @example
* v.fromArray([3,6,5,2]); // v = Vector([3,6,5,2])
*/
JSMatrix.Vector.prototype.fromArray = function(arry) {
	var n = arry.length;
	var e;
	for (e=0; e<n; e++) { this[e] = arry[e]; }
	for (e=n; e<this.nElements; e++) { delete this[e]; }
	this.nElements = n;
};

/** Returns elements of receiver as an Array object
* @returns {Array.<number>} Array object containing receiver's elements
* @example
* v = JSMatrix.vec([3,4,5,6]);
* s = v.toArray(); // s = [3,4,5,6]
*/
JSMatrix.Vector.prototype.toArray = function() {
	var ret = [];
	for (var i=0; i<this.nElements; i++) { ret[i] = this[i]; }
	return ret;
};

/** Returns sum of receiver and vector ( ret = this + vec )
* @param {JSMatrix.Vector} vec 2nd vector of the sum
* @returns {JSMatrix.Vector} sum of receiver and vec
* @example
* v1 = JSMatrix.vec([1,2,3]);
* v2 = JSMatrix.vec([6,1,5]);
* v4 = v1.add(v2); // v4 = Vector([7,3,8])
*/
JSMatrix.Vector.prototype.add = function(vec) {
	if (!this.isSameSizeAs(vec)) { /*jsmLogger.warn('Vector.add: dimensions mismatched');*/ return null; }
	var ret = new JSMatrix.Vector();
	for (var e=0; e<this.nElements; e++) { ret[e] = this[e] + vec[e]; }
	ret.nElements = this.nElements;
	return ret;
};

/** Add given vector to receiver ( this += vec )
* @param {JSMatrix.Vector} vec vector to be added
* @example
* v1 = JSMatrix.vec([1,2,3]);
* v2 = JSMatrix.vec([6,1,5]);
* v1.iadd(v2); // v1 = Vector([7,3,8])
*/
JSMatrix.Vector.prototype.iadd = function(vec) {
	if (!this.isSameSizeAs(vec)) { /*jsmLogger.warn('Vector.iadd: dimensions mismatched');*/ return; }
	for (var e=0; e<this.nElements; e++) { this[e] += vec[e]; }
};

/** Modifies receiver to become sum of two vectors ( this = vec1 + vec2 ). Receiver's size is adjusted
* @param {JSMatrix.Vector} vec1 1st vector of the sum
* @param {JSMatrix.Vector} vec2 2nd vector of the sum
* @example
* v = JSMatrix.vec([3,4]);
* v1 = JSMatrix.vec([1,2,3]);
* v2 = JSMatrix.vec([6,1,5]);
* v.beSumOf(v1,v2); // v = Vector([7,3,8])
*/
JSMatrix.Vector.prototype.beSumOf = function(vec1,vec2) {
	if (!vec1.isSameSizeAs(vec2)) { /*jsmLogger.warn('Vector.beSumOf: dimensions mismatched');*/ return; }
	var e, n = vec1.nElements;
	for (e=0; e<n; e++) { this[e] = vec1[e] + vec2[e]; }
	for (e=n; e<this.nElements; e++) { delete this[e]; }
	this.nElements = n;
};

/** Returns difference of receiver and vector ( ret = this - vec )
* @param {JSMatrix.Vector} vec 2nd vector of subtraction
* @returns {JSMatrix.Vector} difference of receiver and vec
* @example
* v1 = JSMatrix.vec([1,2,3]);
* v2 = JSMatrix.vec([6,1,5]);
* v4 = v1.sub(v2); // v4 = Vector([-5,1,-2])
*/
JSMatrix.Vector.prototype.sub = function(vec) {
	if (!this.isSameSizeAs(vec)) { /*jsmLogger.warn('Vector.sub: dimensions mismatched');*/ return null; }
	var ret = new JSMatrix.Vector();
	for (var e=0; e<this.nElements; e++) { ret[e] = this[e] - vec[e]; }
	ret.nElements = this.nElements;
	return ret;
};

/** Subtract given vector to receiver ( this -= vec )
* @param {JSMatrix.Vector} vec vector to be subtracted
* @example
* v1 = JSMatrix.vec([1,2,3]);
* v1.isub(v2); // v1 = Vector([-5,1,-2])
*/
JSMatrix.Vector.prototype.isub = function(vec) {
	if (!this.isSameSizeAs(vec)) { /*jsmLogger.warn('Vector.isub: dimensions mismatched');*/ return; }
	for (var e=0; e<this.nElements; e++) { this[e] -= vec[e]; }
};

/** Modifies receiver to become difference of two vectors ( this = vec1 - vec2 ). Receiver's size is adjusted
* @param {JSMatrix.Vector} vec1 1st vector of the difference
* @param {JSMatrix.Vector} vec2 2nd vector of the difference
* @example
* v = JSMatrix.vec([4,5]);
* v1 = JSMatrix.vec([1,2,3]);
* v2 = JSMatrix.vec([6,1,5]);
* v.beDifferenceOf(v1,v2); // v = Vector([-5,1,-2])
*/
JSMatrix.Vector.prototype.beDifferenceOf = function(vec1,vec2) {
	if (!vec1.isSameSizeAs(vec2)) { /*jsmLogger.warn('Vector.beDifferenceOf: dimensions mismatched');*/ return; }
	var e, n = vec1.nElements;
	for (e=0; e<n; e++) { this[e] = vec1[e] - vec2[e]; }
	for (e=n; e<this.nElements; e++) { delete this[e]; }
	this.nElements = n;
};

/** Returns dot (inner) product of receiver and given vector
* @param {JSMatrix.Vector} vec vector for dot product
* @param {number=} n =0] dot product will be made from first n elements. If not specified, zero or negative, all elements are used
* @returns {number|null} dot product of receiver and vec from first n elements
* @example
* v1 = JSMatrix.vec([1,2,3]);
* v2 = JSMatrix.vec([6,1,5]);
* d2 = v1.dot(v2,2); // d2 = 8
*/
JSMatrix.Vector.prototype.dot = function(vec,n) {
	n = n || -1;
	if (n<=0 || n>this.nElements) { n = this.nElements; }
	if (n > vec.nElements) { /*jsmLogger.warn('Vector.dot: wrong number of elements');*/ return null; }
	var ret = 0.;
	for (var e=0; e<n; e++) { ret += this[e]*vec[e]; }
	return ret;
};

/** Returns squared (Euclidean) norm of receiver (ret = this.^T * this )
* @returns {number} squared Euclidean norm of receiver
* @example
* v = JSMatrix.vec([4,2,4]);
* n = v.squaredNorm(); // n = 36
*/
JSMatrix.Vector.prototype.squaredNorm = function() {
	var ret = 0., n=this.nElements;
	for (var i=0; i<n; i++) {
		ret += this[i]*this[i];
	}
	return ret;
};

/** Returns (Euclidean) norm of receiver (ret = sqrt(this^T * this) )
* @returns {number} Euclidean norm of receiver
* @example
* v = JSMatrix.vec([4,2,4]);
* n = v.norm(); // n = 6
*/
JSMatrix.Vector.prototype.norm = function() {
	return Math.sqrt(this.dot(this));
};

/** Normlize receiver ( this.norm()=1. )
* @example
* v = JSMatrix.vec([4,2,4]);
* v.normalize(); // v = Vector([0.6666,0.3333,0.6666])
*/
JSMatrix.Vector.prototype.normalize = function() {
	var n = this.norm();
	if (n==0.) { /*jsmLogger.warn('Vector.normalize: norm == 0.');*/ return; }
	this.imulf(1/n);
};

/** Returns normalized copy of receiver ( ret.norm()=1. )
* @returns {JSMatrix.Vector} normalized copy of receiver
* @example
* v1 = JSMatrix.vec([4,2,4]);
* v2 = v1.normalized(); // v2 = Vector([0.6666,0.3333,0.6666])
*/
JSMatrix.Vector.prototype.normalized = function() {
	var n = this.norm();
	if (n==0.) { /*jsmLogger.warn('Vector.normalized: norm == 0.');*/ return null; }
	var e, ret = new JSMatrix.Vector(), nn=1/n;
	for (e=0; e<this.nElements; e++) { ret[e] = this[e]*nn; }
	ret.nElements = this.nElements;
	return ret;
};

/** Returns squared energy norm of receiver (ret = this^T*mat*this)
* @param {JSMatrix.Matrix} mat marix of norm
* @returns {number|null} squared energy norm of receiver
* @example
* v = JSMatrix.vec([1,2,3]);
* m = JSMatrix.mat([[8,0,0],[0,2,1],[0,1,3]]);
* n = v.squaredEnergyNorm(m); // n = 55
*/
JSMatrix.Vector.prototype.squaredEnergyNorm = function(mat) {
	return this.mulm(mat).dot(this);
};

/** Returns energy norm of receiver (ret = sqrt(this*mat*this))
* @param {JSMatrix.Matrix} mat marix of norm
* @returns {number} energy norm of receiver
* @example
* v = JSMatrix.vec([1,2,3]);
* m = JSMatrix.mat([[8,0,0],[0,2,1],[0,1,3]]);
* n = v.energyNorm(m); // n = 7.4162
*/
JSMatrix.Vector.prototype.energyNorm = function(mat) {
	return Math.sqrt(this.squaredEnergyNorm(mat));
};

/** Returns dyadic (tensor, outer, direct,..) product of receiver and given vector
* @param {JSMatrix.Vector} vec vector for multiplication
* @returns {JSMatrix.Matrix} dyadic product of receiver and vec
* @example
* v1 = JSMatrix.vec([1,2,3]);
* v2 = JSMatrix.vec([4,5,6]);
* m = v1.dyadic(v2); // m = Matrix([[4,5,6],[8,10,12],[12,15,18]])
*/
JSMatrix.Vector.prototype.dyadic = function(vec) {
	var ret = new JSMatrix.Matrix();
	for (var r=0; r<this.nElements; r++) {
		ret[r] = {};
		for (var c=0; c<vec.nElements; c++) { ret[r][c] = this[r] * vec[c]; }
	}
	ret.nRows = this.nElements;
	ret.nCols = vec.nElements;
	return ret;
};

/** Alias for dyadic() (inspired by LaTeX symbol for dyadic product operation), see {@link JSMatrix.Vector#dyadic}
* @function
*/
JSMatrix.Vector.prototype.otimes = JSMatrix.Vector.prototype.dyadic;

/** Alias for dyadic(), see {@link JSMatrix.Vector#dyadic}
* @function
*/
JSMatrix.Vector.prototype.outer = JSMatrix.Vector.prototype.dyadic;

/** Alias for dyadic(), see {@link JSMatrix.Vector#dyadic}
* @function
*/
JSMatrix.Vector.prototype.tensor = JSMatrix.Vector.prototype.dyadic;

/** Returns vector (cross) product of receiver and given vector. Both receiver and given vector have to be of length 3.
* @param {JSMatrix.Vector} vec 2nd vector of multiplication
* @returns {JSMatrix.Vector} cross product of receiver and vec
* @example
* v1 = JSMatrix.vec([1,2,3]);
* v2 = JSMatrix.vec([4,5,6]);
* v3 = v1.cross(v2); // v3 = Vector([-3,6,-3])
*/
JSMatrix.Vector.prototype.cross = function(vec) {
	if (this.nElements!=3 || vec.nElements!=3) { /*jsmLogger.warn('Vector.cross: wrong dimensions');*/ return null; }
	return JSMatrix.Vector.create([
		this[1]*vec[2]-this[2]*vec[1],
		this[2]*vec[0]-this[0]*vec[2],
		this[0]*vec[1]-this[1]*vec[0]]);
};

/** Modifies receiver to become cross product of two vectors. Both given vectors have to be of length 3. Receiver's size is adjusted
* @param {JSMatrix.Vector} vec1 1st vector for multiplication
* @param {JSMatrix.Vector} vec2 2nd vector for multiplication
* @example
* v = JSMatrix.vec([5,6]);
* v1 = JSMatrix.vec([1,2,3]);
* v2 = JSMatrix.vec([4,5,6]);
* v.beCrossProductOf(v1,v2); // v = Vector([-3,6,-3])
*/
JSMatrix.Vector.prototype.beCrossProductOf = function(vec1,vec2) {
	if (vec1.nElements!=3 || vec2.nElements!=3) { /*jsmLogger.warn('Vector.beCrossProductOf: wrong dimensions');*/ return; }
	this[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1];
	this[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2];
	this[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0];
	for (var e=3; e<this.nElements; e++) { delete this[e]; }
	this.nElements = 3;
};

/** Returns negative of receiver ( ret = -this )
* @returns {JSMatrix.Vector} negative of receiver
* @example
* v1 = JSMatrix.vec([1,2,3]);
* v2 = v1.negated(); // v2 = Vector([-1,-2,-3])
*/
JSMatrix.Vector.prototype.negated = function() {
	return this.mulf(-1.);
};

/** Alias for negated(), see {@link JSMatrix.Vector#negated}
* @function
*/
JSMatrix.Vector.prototype.neg = JSMatrix.Vector.prototype.negated;

/** Negate receiver ( ret *= -1., ret = -ret )
* @example
* v1 = JSMatrix.vec([1,2,3]);
* v1.negate(); // v1 = Vector([-1,-2,-3])
*/
JSMatrix.Vector.prototype.negate = function() {
	for (var e=0; e<this.nElements; e++) { this[e] *= -1.; }
};

/** Alias for neagete(), see {@link JSMatrix.Vector#negate}
* @function
*/
JSMatrix.Vector.prototype.ineg = JSMatrix.Vector.prototype.negate;

/** Modifies receiver to become negative of given vector (this = -vec). Receiver's size is adjusted
* @param {JSMatrix.Vector} vec given vector
* @example
* v1 = JSMatrix.vec([4,5])
* v2 = JSMatrix.vec([1,2,3]);
* v1.beNegativeOf(v2); // v1 = Vector([-1,-2,-3])
*/
JSMatrix.Vector.prototype.beNegativeOf = function(vec) {
	var e, n = vec.nElements;
	for (e=0; e<n; e++) { this[e] = -vec[e]; }
	for (e=n; e<this.nElements; e++) { delete this[e]; }
	this.nElements = n;
};

/** Returns receiver multiplied by float f ( ret = this * f )
* @param {number} f float multiplier
* @returns {JSMatrix.Vector} copy of receiver multiplied by f
* @example
* v1 = JSMatrix.vec([1,2,3]);
* v2 = v1.mulf(3); // v2 = Vector([3,6,9])
*/
JSMatrix.Vector.prototype.mulf = function(f) {
	var ret = new JSMatrix.Vector();
	for (var e=0; e<this.nElements; e++) { ret[e] = f*this[e]; }
	ret.nElements = this.nElements;
	return ret;
};

/** Multiply receiver by float f ( this *= f )
* @param {number} f float multiplier
* @example
* v = JSMatrix.vec([1,2,3]);
* v.imulf(3); // v = Vector([3,6,9])
*/
JSMatrix.Vector.prototype.imulf = function(f) {
	for (var e=0; e<this.nElements; e++) { this[e] *= f; }
};

/** Returns product of receiver and given matrix ( ret = this^T * mat )
* @param {JSMatrix.Matrix} mat matrix to multiply
* @returns {JSMatrix.Vector} copy of receiver multiplied by mat
* @example
* v1 = JSMatrix.vec([1,2,3]);
* m = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]])
* v2 = v1.mulm(m); // v2 = Vector([146,152,158])
*/
JSMatrix.Vector.prototype.mulm = function(mat) {
	if (!this.canMultiplyMat(mat)) { /*jsmLogger.warn('Vector.mulm: dimensions mismatched');*/ return null; }
	var ret = new JSMatrix.Vector(), temp, nc = mat.nCols, nr = this.nElements;
	for (var c=0; c<nc; c++) {
		temp = 0.;
		for (var r=0; r<nr; r++) {
			temp += this[r]*mat[r][c];
		}
		ret[c] = temp;
	}
	ret.nElements = nc;
	return ret;
};

/** Returns product of receiver and given multiplier (float or Matrix)
* @param {JSMatrix.Matrix|number} what matrix or float to multiply
* @returns {JSMatrix.Vector} copy of receiver multiplied by what
* @example
* v1 = JSMatrix.vec([1,2,3]);
* v2 = v1.mul(JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]])); // v2 = Vector([146,152,158])
* v3 = v1.mul(3); // v3 = Vector([3,6,9])
*/
JSMatrix.Vector.prototype.mul = function(what) {
	if (what instanceof JSMatrix.Matrix) { return this.mulm(what); }
	if ((typeof what).toLowerCase() == 'number') { return this.mulf(what); }
	/*jsmLogger.warn('Vector.mul: wrong argument');*/
	return null;
};

/** Alias for mul(), see {@link JSMatrix.Vector#mul}
* @function
*/
JSMatrix.Vector.prototype.x = JSMatrix.Vector.prototype.mul;

/** Modifies receiver to become product of given matrix and vector ( this =  mat * vec ). Receiver's size is adjusted
* @param {JSMatrix.Matrix} mat matrix of multiplication
* @param {JSMatrix.Vector} vec vector of multiplication
* @example
* v = JSMatrix.vec([4,5]);
* m = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* v1 = JSMatrix.vec([1,2,3]);
* v.beProductOf(m,v1); // v = Vector([74,134,194])
*/
JSMatrix.Vector.prototype.beProductOf = function(mat,vec) {
	if (!mat.canMultiplyVec(vec)) { /*jsmLogger.warn('Vector.beProductOf: dimensions mismatched');*/ return; }
	var nr = mat.nRows;
	var nc = mat.nCols;
	var temp, r, c;
	for (r=0; r<nr; r++) {
		temp = 0.;
		for (c=0; c<nc; c++) {
			temp += mat[r][c]*vec[c];
		}
		this[r] = temp;
	}
	for (r=nr; r<this.nElements; r++) { delete this[r]; }
	this.nElements = nr;
};

/** Modifies receiver to become product of given vector and matrix ( this =  vec^T * mat )
* @param {JSMatrix.Vector} JSMatrix.vec vector of multiplication
* @param {JSMatrix.Matrix} mat matrix of multiplication
* @example
* v= JSMatrix.vec([5,6]);
* v1 = JSMatrix.vec([1,2,3]);
* m = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]])
* v.beTProductOf(v1,m); // v = Vector([146,152,158])
*/
JSMatrix.Vector.prototype.beTProductOf = function(vec,mat) {
if (!vec.canMultiplyMat(mat)) { /*jsmLogger.warn('Vector.beTProductOf: dimensions mismatched');*/ return; }
var nr = mat.nRows;
var nc = mat.nCols;
var temp, r, c;
for (c=0; c<nc; c++) {
	temp = 0.;
	for (r=0; r<nr; r++) {
		temp += vec[r]*mat[r][c];
	}
	this[c] = temp;
};
for (c=nc; c<this.nElements; c++) { delete this[c]; }
this.nElements = nc;
};

/** Returns subvector of receiver. ret = this[coeffs]
* @param {Array.<number>|JSMatrix.Vector} coeffs array containing coefficients of desired subvector
* @returns {JSMatrix.Vector} desired subvector
* @example
* v1 = JSMatrix.vec([4,7,9,1,8,3,2,4,6,7,5,1,6,9]);
* v2 = v1.getSubVector([2,3,6,8,9]); // v2 = Vector([9,1,2,6,7])
*/
JSMatrix.Vector.prototype.getSubVector = function(coeffs) {
	var n = coeffs.nElements==undefined? coeffs.length : coeffs.nElements;
	if (n>this.nElements) { /*jsmLogger.warn('Vector.getSubVector: wrong argument dimensions');*/ return null; }
	var e, ce;
	var ret = new JSMatrix.Vector();
	for (e=0; e<n; e++) {
		ce = coeffs[e];
		if (ce<0 || ce>=this.nElements) { /*jsmLogger.warn('Vector.getSubVector: wrong argument values');*/ return null; }
		ret[e] = this[coeffs[e]];
	}
	ret.nElements = n;
	return ret;
};

/** Alias for getSubVector, see {@link JSMatrix.Vector#getSubVector}
* @function
*/
JSMatrix.Vector.prototype.getsv = JSMatrix.Vector.prototype.getSubVector;

/** Modifies receiver to become subvector of given vector (this = vec[coeffs]). Receiver's size is adjusted
* @param {JSMatrix.Vector} vec vector to get subvector from
* @param {Array.<number>|JSMatrix.Vector} coeffs array containing coefficients of desired subvector
* @example
* v1 = JSMatrix.vec([4,7,9,1,8,3,2,4,6,7,5,1,6,9]);
* v2 = JSMatrix.vec([3,4]);
* v2.beSubVectorOf(v1,[2,3,6,8,9]); // v2 = Vector([9,1,2,6,7])
*/
JSMatrix.Vector.prototype.beSubVectorOf = function(vec,coeffs) {
	var n = coeffs.nElements==undefined? coeffs.length : coeffs.nElements;
	if (n>vec.nElements) { /*jsmLogger.warn('Vector.beSubVectorOf: wrong argument dimensions');*/ return; }
	var e;
	for (e=0; e<n; e++) { if (coeffs[e]<0 || coeffs[e]>=vec.nElements) { /*jsmLogger.warn('Vector.beSubVectorOf: wrong argument values');*/ return; } }
	for (e=0; e<n; e++) { this[e] = vec[coeffs[e]]; }
	for (e=n; e<this.nElements; e++) { delete this[e]; }
	this.nElements = n;
};

/** Sets subvector at defined positions of receiver ( this[coeffs]=vec )
* @param {Array.<number>|JSMatrix.Vector} coeffs array containing coefficients to be set
* @param {JSMatrix.Vector} vec vector to be set
* @example
* v = JSMatrix.vec([4,7,9,1,8,3,2,4,6,7,5,1,6,9]);
* v.setSubVector([2,3,6,8,9],JSMatrix.vec([1.1,2.2,3.3,4.4,5.5])); // v = Vector([4,7,1.1,2.2,8,3,3.3,4,4.4,5.5,5,1,6,9])
*/
JSMatrix.Vector.prototype.setSubVector = function(coeffs,vec) {
	var n = coeffs.nElements==undefined? coeffs.length : coeffs.nElements;
	var vi = vec.nElements;
	if (vi != n) { /*jsmLogger.warn('Vector.setSubVector: wrong argument dimensions');*/ return; }
	if (n>this.nElements) { /*jsmLogger.warn('Vector.setSubVector: wrong argument dimensions');*/ return; }
	var e;
	for (e=0; e<n; e++) { if (coeffs[e]<0 || coeffs[e]>=this.nElements) { /*jsmLogger.warn('Vector.setSubVector: wrong argument values');*/ return; } }
	for (e=0; e<n; e++) { this[coeffs[e]] = vec[e]; }
};

/** Alias for setSubVector, see {@link JSMatrix.Vector#setSubVector}
* @function
*/
JSMatrix.Vector.prototype.setsv = JSMatrix.Vector.prototype.setSubVector;

/** Increment subvector at defined positions of receiver (this[coeffs] += vec)
* @param {Array.<number>|JSMatrix.Vector} coeffs array containing coefficients to be incremented
* @param {JSMatrix.Vector} vec vector to be added
* @example
* v = JSMatrix.vec([4,7,9,1,8,3,2,4,6,7,5,1,6,9]);
* v.incrSubVector([2,3,6,8,9],JSMatrix.vec([1.1,2.2,3.3,4.4,5.5])); // v = Vector([4,7,10.1,3.2,8,3,5.3,4,10.4,12.5,5,1,6,9])
*/
JSMatrix.Vector.prototype.incrSubVector = function(coeffs,vec) {
	var n = coeffs.nElements==undefined? coeffs.length : coeffs.nElements;
	var vi = vec.nElements;
	if (vi != n) { /*jsmLogger.warn('Vector.incrSubVector: wrong argument dimensions');*/ return; }
	if (n>this.nElements) { /*jsmLogger.warn('Vector.incrSubVector: wrong argument dimensions');*/ return; }
	var e;
	for (e=0; e<n; e++) { if (coeffs[e]<0 || coeffs[e]>=this.nElements) { /*jsmLogger.warn('Vector.incrSubVector: wrong argument values');*/ return; } }
	for (e=0; e<n; e++) { this[coeffs[e]] += vec[e]; }
};

/** Alias for incrSubVector, see {@link JSMatrix.Vector#incrSubVector}
* @function
*/
JSMatrix.Vector.prototype.incrsv = JSMatrix.Vector.prototype.incrSubVector;

/** Decrement subvector at defined positions of receiver (this[coeffs] -= vec)
* @param {Array.<number>|JSMatrix.Vector} coeffs array containing coefficients to be decremented
* @param {JSMatrix.Vector} vec vector to be subtracted
* @example
* v = JSMatrix.vec([4,7,9,1,8,3,2,4,6,7,5,1,6,9]);
* v.decrSubVector([2,3,6,8,9],JSMatrix.vec([1.1,2.2,3.3,4.4,5.5])); // v = Vector([4,7,7.9,-1.2,8,3,-1.3,4,1.6,1.5,5,1,6,9])
*/
JSMatrix.Vector.prototype.decrSubVector = function(coeffs,vec) {
	var n = coeffs.nElements==undefined? coeffs.length : coeffs.nElements;
	var vi = vec.nElements;
	if (vi != n) { /*jsmLogger.warn('Vector.decrSubVector: wrong argument dimensions');*/ return; }
	if (n>this.nElements) { /*jsmLogger.warn('Vector.decrSubVector: wrong argument dimensions');*/ return; }
	var e;
	for (e=0; e<n; e++) { if (coeffs[e]<0 || coeffs[e]>=this.nElements) { /*jsmLogger.warn('Vector.decrSubVector: wrong argument values');*/ return; } }
	for (e=0; e<n; e++) { this[coeffs[e]] -= vec[e]; }
};

/** Alias for decrSubVector, see {@link JSMatrix.Vector#decrSubVector}
* @function
*/
JSMatrix.Vector.prototype.decrsv = JSMatrix.Vector.prototype.decrSubVector;

/** assemble receiver to given vector ( vec[coeffs] += this )
* @param {JSMatrix.Vector} vec vector where receiver is assembled
* @param {Array.<number>|JSMatrix.Vector} coeffs array containing coefficients of vec to be incremented
* @example
* v1 = JSMatrix.vec([5,8,2,4,3,5]);
* v2 = JSMatrix.vec([1,2,3]);
* v2.assemble(v1,[5,1,2]); // v1 = Vector([5,10,5,4,3,6])
*/
JSMatrix.Vector.prototype.assemble = function(vec,coeffs) {
	var n = coeffs.nElements==undefined? coeffs.length : coeffs.nElements;
	var vi = this.nElements;
	if (vi != n) { /*jsmLogger.warn('Vector.assemble: wrong argument dimensions');*/ return; }
	if (n>vec.nElements) { /*jsmLogger.warn('Vector.assemble: wrong argument dimensions');*/ return; }
	var e;
	for (e=0; e<n; e++) { if (coeffs[e]<0 || coeffs[e]>=vec.nElements) { /*jsmLogger.warn('Vector.assemble: wrong argument values');*/ return; } }
	for (e=0; e<n; e++) { vec[coeffs[e]] += this[e]; }
};

/** Modifies receiver to become row of given matrix - indexing from 0. Receiver's size is adjusted
* @param {JSMatrix.Matrix} mat given matrix
* @param {number} r desired row (indexing from 0)
* @example
* v = JSMatrix.vec([1,2]);
* m = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* v.beRowOf(m,1) // v = Vector([21,22,23])
*/
JSMatrix.Vector.prototype.beRowOf = function(mat,r) {
	var c, nc = mat.nCols;
	for (c=0; c<nc; c++) { this[c] = mat[r][c]; }
	for (c=nc; c<this.nElements; c++) { delete this[c]; }
	this.nElements = nc;
};

/** Modifies receiver to become column of given matrix - indexing from 0. Receiver's size is adjusted
* @param {JSMatrix.Matrix} mat given matrix
* @param {number} c desired column (indexing from 0)
* @example
* v = JSMatrix.vec([1,2]);
* m = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* v.beColOf(m,1) // v = Vector([12,22,32])
*/
JSMatrix.Vector.prototype.beColOf = function(mat,c) {
	var r, nr = mat.nRows;
	for (r=0; r<nr; r++) { this[r] = mat[r][c]; }
	for (r=nr; r<this.nElements; r++) { delete this[r]; }
	this.nElements = nr;
};

/** Alias for beColOf(), see {@link JSMatrix.Vector#beColOf}
* @function
*/
JSMatrix.Vector.prototype.beColumnOf = JSMatrix.Vector.prototype.beColOf;

/** Returns copy of receiver
* @returns {JSMatrix.Vector} copy of receiver
* @example
* v1 = JSMatrix.vec([1,2,3]);
* v2 = v1;
* v3 = v1.copied();
* v1.set(1,6);
* // v1; = Vector([1,6,3])
* // v2; = Vector([1,6,3])
* // v3; = Vector([1,2,3])
*/
JSMatrix.Vector.prototype.copied = function() {
	var ret = new JSMatrix.Vector(), e, n = this.nElements;
	for (e=0; e<n; e++) { ret[e] = this[e]; }
	ret.nElements = n;
	return ret;
};

/** Alias for copied(), see {@link JSMatrix.Vector#copied}
* @function
*/
JSMatrix.Vector.prototype.copy = JSMatrix.Vector.prototype.copied;

/** Modifies receiver to become copy of given vector (this = vec). Receiver's size is adjusted
* @param {JSMatrix.Vector} vec vector to be copied to receiver
* @example
* v1 = JSMatrix.vec([4,5]);
* v2 = JSMatrix.vec([1,2,3]);
* v1.beCopyOf(v2) // v2 = Vector([1,2,3])
*/
JSMatrix.Vector.prototype.beCopyOf = function(vec) {
	var s = vec.nElements, e;
	for (e=0; e<s; e++) { this[e] = vec[e]; }
	for (e=s; e<this.nElements; e++) { delete this[e]; }
	this.nElements = s;
};

/** Returns string representation of receiver
* @returns {string} string representation of receiver
* @example
* v1 = JSMatrix.vec([0.5999999999999999,2,3]);
* s = v1.toString(); // s = 'Vector([ 1, 2, 3 ])'
*/
JSMatrix.Vector.prototype.toString = function() {
	var nn = this.nElements, e, val;
	if (nn == 0) { return "Vector([])\n"; }
	var ret = "Vector([ ";
	for (e=0; e<nn-1; e++) {
		val = this[e];
		ret += val.toExponential(3)+", ";
	}
	val = this[nn-1];
	ret += val.toExponential(3)+" ])\n";
	return ret;
};

/** Swaps two elements of receiver
* @param {number} e1 index of first element to swap
* @param {number} e2 index of second element to swap
* @example
* v = JSMatrix.vec([1,2,3,4,5,6,7]);
* v.swapElements(1,4); // v = Vector([1,5,3,4,2,6,7])
*/
JSMatrix.Vector.prototype.swapElements = function(e1,e2) {
	if (e1==e2) { return; }
	if (e1<0 || e2<0 || e1>=this.nElements || e2>=this.nElements) { /*jsmLogger.warn('Vector.swapElements: wrong coefficients');*/ return; }
	var temp = this[e1];
	this[e1] = this[e2];
	this[e2] = temp;
};

/** Appends given float/vector to the end of receiver
* @param {number|JSMatrix.Vector|Array.<JSMatrix.Vector>} newElems new element(s) to be appended
* @example
* v = JSMatrix.vec([1,2]);
* v.appendElements(4); // v = Vector([1,2,4]);
* v.appendElements(JSMatrix.vec([7,6])); // v = Vector([1,2,4,7,6]);
* v.appendElements([JSMatrix.vec([1,2]),JSMatrix.vec([3,4])]); // v = Vector([1,2,4,6,7,1,2,3,4])
*/
JSMatrix.Vector.prototype.appendElements = function(newElems) {
	if (newElems instanceof JSMatrix.Vector) {
		var n = newElems.nElements, e;
		for (e=0; e<n; e++) {
			this[this.nElements] = newElems[e];
			this.nElements++;
		}
		return;
	}
	if (newElems instanceof Array) {
		var n = newElems.length, e;
		for (e=0; e<n; e++) {
			this.appendElements(newElems[e]);
		}
		return;
	}
	if ((typeof newElems).toLowerCase() == 'number') {
		this[this.nElements] = newElems;
		this.nElements++;
		return;
	}
};

/** Permute elements of receiver according to given coefficients
* @param {Array.<number>|JSMatrix.Vector} coeffs coefficients of permutation
* @param {boolean=} backward =false] if false, receiver is permutated to given coefficients (forward). If true, from given coefficients (backward)
* @example
* v = JSMatrix.vec([7,9,6]);
* v.permuteElements([1,2,0]); // v = Vector([9,6,7])
* v.permuteElements([1,2,0],true); // v = Vector([7,9,6])
*/
JSMatrix.Vector.prototype.permuteElements = function(coeffs,backward) {
	var n;
	if (coeffs.nElements==undefined) { // coeffs == Array
		n = coeffs.length;
		coeffs = coeffs.slice();
	} else { // coeffs == Vector
		n = coeffs.nElements;
		coeffs = coeffs.toArray();
	}
	if (n != this.nElements) { /*jsmLogger.warn('Vector.permuteElements: wrong argument size');*/ return; }
	for (var e=0; e<n; e++) { if (e<0 || e>=this.nElements) { /*jsmLogger.warn('Vector.permuteElements: wrong argument values');*/ return; } }
	backward = backward || false;
	var ee;
	if (backward) {
		for (var e=0; e<n; e++) {
			ee = coeffs[e];
			while (ee!=e) {
				this.swapElements(e,ee);
				coeffs[e] = coeffs[ee];
				coeffs[ee] = ee;
				ee = coeffs[e];
			}
		}
	} else {
		for (var e=n-1; e>0; e--) {
			ee = coeffs[e];
			this.swapElements(e,ee);
			coeffs[coeffs.indexOf(e)] = ee;
		}
	}
	return;
};

/** Returns permutated copy of receiver according to given coefficients
* @param {Array.<number>|JSMatrix.Vector} coeffs coefficients of permutation. The content will be changed(!), use a copy if you need to preserve the content
* @param {boolean=} backward =false] if false, receiver is permutated to given coefficients (forward). If true, from given coefficients (backward)
* @returns {JSMatrix.Vector} permutated copy of receiver
* @example
* v1 = JSMatrix.vec([7,9,6]);
* v2 = v1.elementsPermutation([1,2,0]); // v2 = Vector([9,6,7])
* v3 = v1.elementsPermutation([1,2,0],true); // v3 = Vector([6,7,9])
*/
JSMatrix.Vector.prototype.elementsPermutation = function(coeffs,backward) {
	var n = coeffs.nElements==undefined? coeffs.length : coeffs.nElements;
	if (n != this.nElements) { /*jsmLogger.warn('Vector.elementsPermutation: wrong argument size');*/ return null; }
	for (var e=0; e<n; e++) { if (e<0 || e>=this.nElements) { /*jsmLogger.warn('Vector.elementsPermutation: wrong argument values');*/ return null; } }
	var ret = new JSMatrix.Vector();
	if (backward) {
		for (e=0; e<this.nElements; e++) {
			ret[coeffs[e]] = this[e];
		}
	} else {
		for (e=0; e<this.nElements; e++) {
			ret[e] = this[coeffs[e]];
		}
	}
	ret.nElements = this.nElements;
	return ret;
};

/** Resize receiver according to given size (delete extra elements or add zero elements)
* @param {number} nElements new number of elements
* @example
* v1 = JSMatrix.vec([4,7,9,1,7,3])
* v2 = JSMatrix.vec([4,6]);
* v1.resize(4); // v1 = Vector([4,7,9,1])
* v2.resize(4); // v2 = Vector([4,6,0,0])
*/
JSMatrix.Vector.prototype.resize = function(nElements) {
	if (nElements==this.nElements) { return; }
	if (nElements > this.nElements) {
		for (var e=this.nElements; e<nElements; e++) { this[e] = 0.; }
		this.nElements = nElements;
		return;
	}
	if (nElements < this.nElements) {
		for (var e=nElements; e<this.nElements; e++) { delete this[e]; }
		this.nElements = nElements;
		return;
	}
};

/** Returns resized copy of receiver according to given size (delete extra elements or add zero elements)
* @param {number} nElements new number of elements
* @returns {JSMatrix.Vector} resized copy of receiver
* @example
* v1 = JSMatrix.vec([4,7,9,1]);
* v2 = v1.resized(2); // v2 = Vector([4,7])
* v3 = v1.resized(6); // v3 = Vector([4,7,9,1,0,0])
*/
JSMatrix.Vector.prototype.resized = function(nElements) {
	var e, ret = new JSMatrix.Vector();
	if (nElements<=this.nElements) {
		for (e=0; e<nElements; e++) { ret[e] = this[e]; }
	} else {
		for (e=0; e<this.nElements; e++) { ret[e] = this[e]; }
		for (e=this.nElements; e<nElements; e++) { ret[e] = 0.; }
	}
	ret.nElements = nElements;
	return ret;
};

/** Returns matrix with receiver's elements on its diagonal
* @returns {JSMatrix.Matrix} diagonal matrix with receiver's elements on its diagonal
* @example
* v = JSMatrix.vec([1,2,3]);
* m = v.toDiagonalMatrix(); // m = Matrix([[1,0,0],[0,2,0],[0,0,3]])
*/
JSMatrix.Vector.prototype.toDiagonalMatrix = function() {
	var ret = JSMatrix.Matrix.Zeros(this.nElements,this.nElements);
	for (var e=0; e<this.nElements; e++) { ret[e][e] = this[e]; }
	return ret;
};

/** Modifies receiver to contain diagonal elements of mat ( this = diag(mat) )
* @param {JSMatrix.Matrix} mat matrix whose diagonal is copied to receiver
* @example
* v = JSMatrix.vec([6,7]);
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* v.beDiagonalFrom(m); // v = Vector([1,5,9])
*/
JSMatrix.Vector.prototype.beDiagonalFrom = function(mat) {
	if (!mat.isSquare()) { /*jsmLogger.warn('Vector.beDiagonalFrom: non-square matrix');*/ return; }
	var e, s = mat.nRows;
	for (e=0; e<s; e++) { this[e] = mat[e][e]; }
	for (e=s; s<this.nElements; e++) { delete this[e]; }
	this.nElements = s;
};

/** Returns matrix with receiver's elements on diagonal (if no argument is specified, see{@link JSMatrix.Vector#toDiagonalMatrix}) or modifies receiver to contain diagonal elements of given matrix (if given matrix is specified, see{@link JSMatrix.Vector#beDiagonalFrom})
* @param {JSMatrix.Matrix=} mat =undefined] matrix whose diagonal is copied to receiver. If not specified, diagnal matrix is returned
* @returns {JSMatrix.Matrix|undefined} if no input parameter, returns diagonal matrix, otherwise nothing
* @example
* v = JSMatrix.vec([6,7]);
* m1 = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* v.diag(m1); // v = Vector([1,5,9])
* m2 = v.diag(); // m2 = Matrix([[1,0,0],[0,5,0],[0,0,9]])
*/
JSMatrix.Vector.prototype.diag = function(mat) {
	if (mat==undefined) { return this.toDiagonalMatrix(); }
	this.beDiagonalFrom(mat);
};

/** Modifies receiver to be a solution of linear system of equations mat*this = rhs. Receiver's size is adjusted
* @param {JSMatrix.Matrix} mat matrix of the system
* @param {JSMatrix.Vector} rhs vector of right hand side
* @param {string=} method ="default"] see {@JSMatrix.Matrix#linSolve}
* @param {boolean=} saveOrig =true] see {@JSMatrix.Matrix#linSolve}
* @param {Array.<JSMatrix.Matrix>=} precompDecomps =undefined] see {@JSMatrix.Matrix#linSolve}
* @example
* v = JSMatrix.vec([2,3]);
* a = JSMatrix.mat([[1,2,9],[8,3,2],[3,7,3]]);
* b = JSMatrix.vec([32,20,26]);
* v.beSolutionOf(a,b); // v = Vector([1,2,3])
*/
JSMatrix.Vector.prototype.beSolutionOf = function(mat,rhs,method,saveOrig,precompDecomps) {
	saveOrig = saveOrig==undefined? true : saveOrig;
	mat = saveOrig? mat.copy() : mat;
	this.beCopyOf(rhs);
	mat.linSolve(this,method,false,precompDecomps);
};

/** Constructs new Vector from given array
* @param {Array.<number>} arry (default=[]) array containing elements of new vector
* @returns {JSMatrix.Vector} new Vector object 
* @example
* v = JSMatrix.Vector.create([1,2,3,4]); // v = Vector([1,2,3,4])
*/
JSMatrix.Vector.create = function(arry) {
	var a = arry || [];
	var n = a.length;
	var ret = new JSMatrix.Vector();
	for (var i=0; i<n; i++) { ret[i] = a[i]; }
	ret.nElements = n;
	return ret;
};

/** Creates unit vector in x direction
* @returns {JSMatrix.Vector} unit x vector
* @example
* v = JSMatrix.Vector.UnitX(); // v = Vector([1,0,0])
*/
JSMatrix.Vector.UnitX = function() { return JSMatrix.Vector.create([1,0,0]); }

/** Creates unit vector in y direction
* @returns {JSMatrix.Vector} unit y vector
* @example
* v = JSMatrix.Vector.UnitY(); // v = Vector([0,1,0])
*/
JSMatrix.Vector.UnitY = function() { return JSMatrix.Vector.create([0,1,0]); }

/** Creates unit vector in z direction
* @returns {JSMatrix.Vector} unit z vector
* @example
* v = JSMatrix.Vector.UnitZ(); // v = Vector([0,0,1])
*/
JSMatrix.Vector.UnitZ = function() { return JSMatrix.Vector.create([0,0,1]); }

/** Creates a vector full of zeros
* @param {number=} nElements =0] number of elements
* @returns {JSMatrix.Vector} new vector full of zeros
* @example
* v = JSMatrix.Vector.Zeros(4); // v = Vector([0,0,0,0])
*/
JSMatrix.Vector.Zeros = function(nElements) {
	nElements = nElements || 0;
	return new JSMatrix.Vector(nElements)
};

/** Creates a vector full of ones
* @param {number=} nElements =0] number of elements
* @returns {JSMatrix.Vector} new vector full of zeros
* @example
* v = JSMatrix.Vector.Ones(6); // v = Vector([1,1,1,1,1,1])
*/
JSMatrix.Vector.Ones = function(nElements) {
	nElements = nElements || 0;
	var r, ret = new JSMatrix.Vector();
	for (r=0; r<nElements; r++) { ret[r] = 1.; }
	ret.nElements = nElements;
	return ret;
};









/** Empties receiver (resizes it to size 0)
* @example
* m = JSMatrix.mat([[1,2],[3,4]]);
* m.empty(); // m = Matrix([[]])
*/
JSMatrix.Matrix.prototype.empty = function() {
	this.resize(0,0);
};

/** Tests if receiver is empty (has size 0,0)
* @returns true if receiver is empty, false otherwise
* @example
* v1 = new JSMatrix.Matrix(2,3);
* v2 = new JSMatrix.Matrix();
* b1 = v1.isEmpty(); // b1 = false
* b2 = v2.isEmpty(); // b2 = true
*/
JSMatrix.Matrix.prototype.isEmpty = function() {
	return this.nRows==0 && this.nCols==0;
};

/** Zeroes all elements of receiver
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* m.zero(); // m = Matrix([[0,0,0],[0,0,0],[0,0,0]])
*/
JSMatrix.Matrix.prototype.zero = function() {
	for (var r=0; r<this.nRows; r++) {
		for (var c=0; c<this.nCols; c++) { this[r][c] = 0.; }
	}
};

/** Sets all elements of receiver to 1.0
* @example
* m = JSMatrix.mat([[2,3,4],[5,6,7]]);
* m.beFullOfOnes(); // m = Matrix([[1,1,1],[1,1,1]])
*/
JSMatrix.Matrix.prototype.beFullOfOnes = function() {
	for (var r=0; r<this.nRows; r++) {
		for (var c=0; c<this.nCols; c++) { this[r][c] = 1.; }
	}
};

/** Alias for beFullOfOnes(), see {@JSMatrix.Matrix#beFullOfOnes}
* @function
*/
JSMatrix.Matrix.prototype.one = JSMatrix.Matrix.prototype.beFullOfOnes;

/** Modifies receiver to be unit matrix
* @param {number=} newSize (default=0) new size or receiver. If omitted or 0, current size is considered (then required square matrix), else resized to newSize
* @example
* m = JSMatrix.mat([[1],[2]]);
* m.beUnitMatrix(); // m = Matrix([[1],[2]])
* /// nothing happened, m is not square
* m.beUnitMatrix(2); // m = Matrix([[1,0],[0,1]])
* m = JSMatrix.mat([[2,3],[4,5]]);
* m.beUnitMatrix(); // m = Matrix([[1,0],[0,1]])
*/
JSMatrix.Matrix.prototype.beUnitMatrix = function(newSize) {
	newSize = newSize===undefined? 0 : newSize;
	if (newSize<=0) {
		if (!this.isSquare()) { /*jsmLogger.warn('Matrix.beUnitMatrix: matrix is not square');*/ return; }
		this.zero();
		for (var i=0; i<this.nRows; i++) { this[i][i] = 1.; }
		return;
	}
	var r,c;
	this.resize(newSize,newSize);
	for (r=0; r<newSize; r++) { this[r][r] = 1.; }
};

/** Returns size of receiver as [nRows,nCols]
* @returns {{nRows:number,nCols:number}} {nRows:number,nCols:number} size of receiver
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6]]);
* s = m.size(); // s = {nRows:2,nCols:3}
* nRows = s.nRows; // nRows = 2
* nCols = s.nCols; // nCols = 3
*/
JSMatrix.Matrix.prototype.size = function() {
	return {nRows:this.nRows,nCols:this.nCols};
};

/** Returns trace (sum of diagonal elements) of the receiver
* @returns {number|null} trace of the receiver
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]])
* tr = m.trace() // tr = 15
*/
JSMatrix.Matrix.prototype.trace = function() {
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.trace: matrix is not square');*/ return null; }
	var ret = 0.;
	for (var e=0; e<this.nRows; e++) {
		ret += this[e][e];
	}
	return ret;
};

/** Alias for trace, see {@JSMatrix.Matrix#trace}
* @function
*/
JSMatrix.Matrix.prototype.diagSum = JSMatrix.Matrix.prototype.trace;

/** Returns product of diagonal elements of the receiver
* @returns {number|null} product of diagonal elements of the receiver
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* p = m.diagProduct(); // p = 45
*/
JSMatrix.Matrix.prototype.diagProduct = function() {
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.trace: matrix is not square');*/ return null; }
	var ret = 1.;
	for (var e=0; e<this.nRows; e++) {
		ret *= this[e][e];
		if (ret==0.) { return 0.; }
	}
	return ret;
};

/** Alias for diagProduct, see {@JSMatrix.Matrix#diagProduct}
* @function
*/
JSMatrix.Matrix.prototype.dprod = JSMatrix.Matrix.prototype.diagProduct;

/** Coefficient access function, returns one element of receiver ( ret = this[r][c] ). Matrix.get(r,c) is much safer than direct Matrix[r][c] access method
* @param {number} r index of row of element to be returned
* @param {number} c index of column of element to be returned
* @returns {number|null} value of the element at r-th row and c-th column
* @example
* m = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* g = m.get(1,2); // g = 23
* g = m[1][2]; // g = 23
* /// c out of bounds
* g = m.get(1,5); // g = null
* g = m[1][5]; // g = undefined
* /// r out of bounds
* g = m.get(6,5) // g = null
* /// g = m[6][5]; would be error (m[6] is not defined)
*/
JSMatrix.Matrix.prototype.get = function(r,c) {
	if (r<0 || c<0) { /*jsmLogger.warn('Matrix.get: coefficients out of bounds');*/ return null; }
	if (r>=this.nRows || c>=this.nCols) { return null; }
	return this[r][c];
};

/** Coefficient access function, sets one element of receiver ( this[r][c] = val ). Matrix.set(r,c,val) is much safer than direct Matrix[r][c]=val access method
* @param {number} r index of row of element to be set
* @param {number} c index of column of element to be set
* @param {number} val value to be set
* @example
* m = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* m.set(1,2,3.45); // m = Matrix([[11,12,13],[21,22,3.45],[31,32,33]])
* m[1][2] = 3.45; // m = Matrix([[11,12,13],[21,22,3.45],[31,32,33]])
* m = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* /// c out of bounds (nothing is done by both appraches)
* m.set(1,5,8.45); // m = Matrix([[11,12,13],[21,22,3.45],[31,32,33]])
* m[1][5] = 8.45; // m = Matrix([[11,12,13],[21,22,3.45],[31,32,33]])
* /// r out of bounds
* m.set(6,5,8.45); // m = Matrix([[11,12,13],[21,22,3.45],[31,32,33]])
* /// m[6][5] = 8.45; would be error (m[6] is not defined)
*/
JSMatrix.Matrix.prototype.set = function(r,c,val) {
	if (r<0 || c<0) { /*jsmLogger.warn('JSMatrix.Matrix.set: coefficients out of bounds');*/ return; }
	if (r>=this.nRows || c>=this.nCols) { return; }
	this[r][c] = val;
};

/** Coefficient access function, increment one element of receiver ( this[r][c] += val ). Matrix.incr(r,c,val) is much safer than direct Matrix[r][c]+=val access method
* @param {number} r index of row of element to be incremented
* @param {number} c index of column of element to be incremented
* @param {number} val value to be add
* @example
* m = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* m.incr(1,2,3.45); // m = Matrix([[11,12,13],[21,22,26.45],[31,32,33]])
* m[1][2] += 3.45; // m = Matrix([[11,12,13],[21,22,29.9],[31,32,33]])
* /// c out of bounds
* m.incr(1,5,3.45); // m = Matrix([[11,12,13],[21,22,29.9],[31,32,33]])
* /// m[1][5] += 3.45); would be error (m[1][5] is not defined)
* /// r out of bounds
* m.incr(6,5,3.45); // m = Matrix([[11,12,13],[21,22,29.9],[31,32,33]])
* /// m[6][5] += 3.45; would be error (m[6] is not defined)
*/
JSMatrix.Matrix.prototype.incr = function(r,c,val) {
	if (r<0 || c<0) { /*jsmLogger.warn('Matrix.incr: coefficients out of bounds');*/ return; }
	if (r>=this.nRows || c>=this.nCols) { return; }
	if (val) { this[r][c] += val; }
};

/** Sets elements of receiver form given Array object
* @param {Array.<Array.<number>>} arry 2D array containing new receiver elements
* @example
* m = JSMatrix.mat([[3,2],[1,1]])
m.fromArray([[1,2,3],[4,5,6],[7,8,9]]); // m = Matrix([[1,2,3],[4,5,6],[7,8,9]])
*/
JSMatrix.Matrix.prototype.fromArray = function(arry) {
	var r, c, nr = arry.length, nc = arry[0].length;
	for (r=0; r<nr; r++) {
		if (r >= this.nRows) { this[r] = {}; }
		for (c=0; c<nc; c++) {
			this[r][c] = arry[r][c];
		}
		for (c=nc; c<this.nCols; c++) { delete this[r][c]; }
	}
	for (r=nr; r<this.nRows; r++) { delete this[r]; }
	this.nRows = nr;
	this.nCols = nc;
};

/** Returns elements of receiver as an Array object
* @returns {Array.<Array.<number>>} 2D array containing receiver's elements
* @example
* m =JSMatrix.mat([[11,12],[21,22]]);
* a = m.toArray(); // a = [[11,12],[21,22]]
*/
JSMatrix.Matrix.prototype.toArray = function() {
	var r, c, ret = []
	for (r=0; r<this.nRows; r++) {
		ret[r] = [];
		for (c=0; c<this.nCols; c++) {
			ret[r][c] = this[r][c];
		}
	}
	return ret;
};

/** Returns sum of receiver and matrix ( ret = this + mat )
* @param {JSMatrix.Matrix} mat matrix to be added
* @returns {JSMatrix.Matrix} sum of receiver and matrix
* @example
* m1 = JSMatrix.mat([[1,2],[3,4]]);
* m2 = JSMatrix.mat([[2,5],[3,2]]);
* m3 = m1.add(m2); // m3 = Matrix([[3,7],[6,6]])
*/
JSMatrix.Matrix.prototype.add = function(mat) {
	if (!this.isSameSizeAs(mat)) { /*jsmLogger.warn('Matrix.add: dimensions mismatched');*/ return null; }
	var r, c, nr = this.nRows, nc = this.nCols, ret = new JSMatrix.Matrix();
	for (r=0; r<nr; r++) {
		ret[r] = {};
		for (c=0; c<nc; c++) { ret[r][c] = this[r][c] + mat[r][c]; }
	}
	ret.nRows = nr;
	ret.nCols = nc;
	return ret;
};

/** Add given matrix to receiver ( this += mat )
* @param {JSMatrix.Matrix} mat matrix to be added
* @example
* m1 = JSMatrix.mat([[1,2],[3,4]]);
* m2 = JSMatrix.mat([[2,5],[3,2]]);
* m1.iadd(m2); // m1 = Matrix([[3,7],[6,6]])
*/
JSMatrix.Matrix.prototype.iadd = function(mat) {
	if (!this.isSameSizeAs(mat)) { /*jsmLogger.warn('Matrix.iadd: dimensions mismatched');*/ return; }
	for (var r=0; r<this.nRows; r++) {
		for (var c=0; c<this.nCols; c++) { this[r][c] += mat[r][c]; }
	}
};

/** Modifies receiver to become sum of two matrices ( this = mat1 + mat2 ). Receiver's size is adjusted
* @param {JSMatrix.Matrix} mat1 1st matrix of the sum
* @param {JSMatrix.Matrix} mat2 2nd matrix of the sum
* @example
* m = JSMatrix.mat([[5],[9]]);
* m1 = JSMatrix.mat([[1,2],[3,4]]);
* m2 = JSMatrix.mat([[2,5],[3,2]]);
* m.beSumOf(m1,m2); // m = Matrix([[3,7],[6,6]])
*/
JSMatrix.Matrix.prototype.beSumOf = function(mat1,mat2) {
	if (!mat1.isSameSizeAs(mat2)) { /*jsmLogger.warn('Matrix.beSumOf: dimensions mismatched');*/ return; }
	var r, c , nr = mat1.nRows, nc = mat1.nCols;
	for (r=0; r<nr; r++) {
		if (r >= this.nRows) { this[r] = {}; }
		for (c=0; c<nc; c++) { this[r][c] = mat1[r][c] + mat2[r][c]; }
		for (c=nc; c<this.nCols; c++) { delete this[r][c]; }
	}
	for (r=nr; r<this.nRows; r++) { delete this[r]; }
	this.nRows = nr;
	this.nCols = nc;
};

/** Returns difference of receiver and matrix ( ret = this - mat )
* @param {JSMatrix.Matrix} mat matrix to be added
* @returns {JSMatrix.Matrix} sum of receiver and matrix
* @example
* m1 = JSMatrix.mat([[1,2],[3,4]]);
* m2 = JSMatrix.mat([[2,5],[3,2]]);
* m3 = m1.sub(m2); // m3 = Matrix([[-1,-3],[0,2]])
*/
JSMatrix.Matrix.prototype.sub = function(mat) {
	if (!this.isSameSizeAs(mat)) { /*jsmLogger.warn('Matrix.sub: dimensions mismatched');*/ return null; }
	var r, c, nr = this.nRows, nc = this.nCols, ret = new JSMatrix.Matrix();
	for (r=0; r<nr; r++) {
		ret[r] = {};
		for (c=0; c<nc; c++) { ret[r][c] = this[r][c] - mat[r][c]; }
	}
	ret.nRows = nr;
	ret.nCols = nc;
	return ret;
};

/** Subtract given matrix to receiver ( this -= mat )
* @param {JSMatrix.Matrix} mat matrix to be added
* @example
* m1 = JSMatrix.mat([[1,2],[3,4]]);
* m2 = JSMatrix.mat([[2,5],[3,2]]);
* m1.isub(m2); // m1 = Matrix([[-1,-3],[0,2]])
*/
JSMatrix.Matrix.prototype.isub = function(mat) {
	if (!this.isSameSizeAs(mat)) { /*jsmLogger.warn('Matrix.isub: dimensions mismatched');*/ return; }
	for (var r=0; r<this.nRows; r++) {
		for (var c=0; c<this.nCols; c++) { this[r][c] -= mat[r][c]; }
	}
};

/** Modifies receiver to become difference of two matrices ( this = mat1 - mat2 ). Receiver's size is adjusted
* @param {JSMatrix.Matrix} mat1 1st matrix of the difference
* @param {JSMatrix.Matrix} mat2 2nd matrix of the difference
* @example
* m = JSMatrix.mat([[5],[9]]);
* m1 = JSMatrix.mat([[1,2],[3,4]]);
* m2 = JSMatrix.mat([[2,5],[3,2]]);
* m.beDifferenceOf(m1,m2); // m = Matrix([[-1,-3],[0,2]])
*/
JSMatrix.Matrix.prototype.beDifferenceOf = function(mat1,mat2) {
	if (!mat1.isSameSizeAs(mat2)) { /*jsmLogger.warn('Matrix.beDifferenceOf: dimensions mismatched');*/ return; }
	var r, c , nr = mat1.nRows, nc = mat1.nCols;
	for (r=0; r<nr; r++) {
		if (r >= this.nRows) { this[r] = {}; }
		for (c=0; c<nc; c++) { this[r][c] = mat1[r][c] - mat2[r][c]; }
		for (c=nc; c<this.nCols; c++) { delete this[r][c]; }
	}
	for (r=nr; r<this.nRows; r++) { delete this[r]; }
	this.nRows = nr;
	this.nCols = nc;
};

/** Returns negative of receiver ( ret = -this )
* @returns {JSMatrix.Matrix} negative of receiver
* @example
* m1 = JSMatrix.mat([[1,2,3],[4,5,6]]);
* m2 = m1.negated(); // m2 = Matrix([[-1,-2,-3],[-4,-5,-6]])
*/
JSMatrix.Matrix.prototype.negated = function() {
	return this.mulf(-1.);
};

/** Alias for neageted(), see {@JSMatrix.Matrix#negated}
* @function
*/
JSMatrix.Matrix.prototype.neg = JSMatrix.Matrix.prototype.negated;

/** Negate receiver ( ret *= -1., ret = -ret )
* @example
* m1 = JSMatrix.mat([[1,2,3],[4,5,6]]);
* m1.negate(); // m1 = Matrix([[-1,-2,-3],[-4,-5,-6]])
*/
JSMatrix.Matrix.prototype.negate = function() {
	for (var r=0; r<this.nRows; r++) {
		for (var c=0; c<this.nCols; c++) { this[r][c] *= -1.; }
	}
};

/** Alias for neagete(), see {@JSMatrix.Matrix#negate}
* @function
*/
JSMatrix.Matrix.prototype.ineg = JSMatrix.Matrix.prototype.negate;

/** Modifies receiver to become negative of given matrix (this = -vec). Receiver's size is adjusted
* @param {JSMatrix.Matrix} mat given matrix
* @example
* m1 = JSMatrix.mat([[2],[3]]);
* m2 = JSMatrix.mat([[1,2,3],[4,5,6]]);
* m1.beNegativeOf(m2); // m1 = Matrix([[-1,-2,-3],[-4,-5,-6]])
*/
JSMatrix.Matrix.prototype.beNegativeOf = function(mat) {
	var r, c , nr = mat.nRows, nc = mat.nCols;
	for (r=0; r<nr; r++) {
		if (r >= this.nRows) { this[r] = {}; }
		for (c=0; c<nc; c++) { this[r][c] = -mat[r][c]; }
		for (c=nc; c<this.nCols; c++) { delete this[r][c]; }
	}
	for (r=nr; r<this.nRows; r++) { delete this[r]; }
	this.nRows = nr;
	this.nCols = nc;
};

/** Returns one row of receiver
* @param {number} r index of row to return
* @returns {JSMatrix.Vector} r-th row of receiver as vector
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* v = m.getRow(1); // v = Vector([4,5,6])
*/
JSMatrix.Matrix.prototype.getRow = function(r) {
	var ret = new JSMatrix.Vector();
	for (var c=0; c<this.nCols; c++) { ret[c] = this[r][c]; }
	ret.nElements = this.nCols;
	return ret;
};

/** Returns one column of receiver
* @param {number} c index of column to return
* @returns {JSMatrix.Vector} c-th column of receiver as vector
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* v = m.getCol(1); // v = Vector([2,5,8])
*/
JSMatrix.Matrix.prototype.getCol = function(c) {
	var ret = new JSMatrix.Vector()
	for (var r=0; r<this.nRows; r++) { ret[r] = this[r][c]; }
	ret.nElements = this.nRows;
	return ret;
};

/** Sets one row of receiver
* @param {number} r index of row to set
* @param {JSMatrix.Vector} vec vector to be set as new r-th row
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* m.setRow(1,JSMatrix.vec([11,12,13])); // m = Matrix([[1,2,3],[11,12,13],[7,8,9]])
*/
JSMatrix.Matrix.prototype.setRow = function(r,vec) {
	if (r<0 || r>=this.nRows) { /*jsmLogger.warn('Matrix.setRow: wrong row');*/ return; }
	var n = vec.nElements;
	if (this.nCols != n) { /*jsmLogger.warn('Matrix.setRow: wrong argument');*/ return; }
	for (var c=0; c<this.nCols; c++) { this[r][c] = vec[c]; }
};

/** Sets one column of receiver
* @param {number} c index of column to set
* @param {JSMatrix.Vector} vec vector to be set as new c-th column
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* m.setCol(1,JSMatrix.vec([11,12,13])); // m = Matrix([[1,11,3],[4,12,6],[7,13,9]])
*/
JSMatrix.Matrix.prototype.setCol = function(c,vec) {
	if (c<0 || c>=this.nCols) { /*jsmLogger.warn('Matrix.setCol: wrong column');*/ return; }
	var n = vec.nElements;
	if (this.nRows != n) { /*jsmLogger.warn('Matrix.setCol: wrong argument');*/ return; }
	for (var r=0; r<this.nRows; r++) { this[r][c] = vec[r]; }
};

/** Incerements one row of receiver
* @param {number} r index of row to set
* @param {JSMatrix.Vector} vec vector to be incremented to r-th row
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* m.incrRow(1,JSMatrix.vec([11,12,13])); // m = Matrix([[1,2,3],[15,17,19],[7,8,9]])
*/
JSMatrix.Matrix.prototype.incrRow = function(r,vec) {
	if (r<0 || r>=this.nRows) { /*jsmLogger.warn('Matrix.incrRow: wrong column');*/ return; }
	var n = vec.nElements;
	if (this.nCols != n) { /*jsmLogger.warn('Matrix.incrRow: wrong argument');*/ return; }
	for (var c=0; c<this.nCols; c++) { this[r][c] += vec[c]; }
};

/** Incerements one column of receiver
* @param {number} c index of column to set
* @param {JSMatrix.Vector} vec vector to be incremented to c-th column
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* m.incrCol(1,JSMatrix.vec([11,12,13])); // m = Matrix([[1,13,3],[4,17,6],[7,21,9]])
*/
JSMatrix.Matrix.prototype.incrCol = function(c,vec) {
	if (c<0 || c>=this.nCols) { /*jsmLogger.warn('Matrix.incrCol: wrong column');*/ return; }
	var n = vec.nElements;
	if (this.nRows != n) { /*jsmLogger.warn('Matrix.incrCol: wrong argument');*/ return; }
	for (var r=0; r<this.nRows; r++) { this[r][c] += vec[r]; }
};

/** Decerements one row of receiver
* @param {number} r index of row to set
* @param {JSMatrix.Vector} vec vector to be decremented from r-th row
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* m.decrRow(1,JSMatrix.vec([11,12,13])); // m = Matrix([[1,2,3],[-7,-7,-7],[7,8,9]])
*/
JSMatrix.Matrix.prototype.decrRow = function(r,vec) {
	if (r<0 || r>=this.nRows) { /*jsmLogger.warn('Matrix.decrRow: wrong column');*/ return; }
	var n = vec.nElements;
	if (this.nCols != n) { /*jsmLogger.warn('Matrix.decrRow: wrong argument');*/ return; }
	for (var c=0; c<this.nCols; c++) { this[r][c] -= vec[c]; }
};

/** Decerements one column of receiver
* @param {number} c index of column to set
* @param {JSMatrix.Vector} vec vector to be decremented from c-th column
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* m.decrCol(1,JSMatrix.vec([11,12,13])); // m = Matrix([[1,-9,3],[4,-7,6],[7,-5,9]])
*/
JSMatrix.Matrix.prototype.decrCol = function(c,vec) {
	if (c<0 || c>=this.nCols) { /*jsmLogger.warn('Matrix.decrCol: wrong column');*/ return; }
	var n = vec.nElements;
	if (this.nRows != n) { /*jsmLogger.warn('Matrix.decrCol: wrong argument');*/ return; }
	for (var r=0; r<this.nRows; r++) { this[r][c] -= vec[r]; }
};

/** Swaps two rows of receiver
* @param {number} r1 index of first row to swap
* @param {number} r2 index of second row to swap
* @example
* m = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* m.swapRows(0,2); // m = Matrix([[31,32,33],[21,22,23],[11,12,13]])
*/
JSMatrix.Matrix.prototype.swapRows = function(r1,r2) {
	if (r1==r2) { return; }
	if (r1<0 || r1>=this.nRows || r2<0 || r2>=this.nRows) { /*jsmLogger.warn('Matrix.swapRows: wrong coefficients');*/ return; }
	var temp = this[r1];
	this[r1] = this[r2];
	this[r2] = temp;
};

/** Swaps two columns of receiver
* @param {number} c1 index of first row to swap
* @param {number} c2 index of second row to swap
* @example
* m = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* m.swapCols(0,2); // m = Matrix([[13,12,11],[23,22,21],[33,32,31]])
*/
JSMatrix.Matrix.prototype.swapCols = function(c1,c2) {
	if (c1==c2) { return; }
	if (c1<0 || c1>=this.nCols || c2<0 || c2>=this.nCols) { /*jsmLogger.warn('Matrix.swapCols: wrong coefficients');*/ return; }
	var temp;
	for (var r=0; r<this.nRows; r++) {
		temp = this[r][c1];
		this[r][c1] = this[r][c2];
		this[r][c2] = temp;
	}
};

/** Permute rows of receiver according to given coefficients
* @param {Array.<number>|JSMatrix.Vector} coeffs coefficients of row permutation
* @param {boolean=} backward =false] if false, receiver's rows is permutated to given coefficients (forward). If true, from given coefficients (backward)
* @example
* m = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* m.permuteRows([1,2,0]); // m = Matrix([[21,22,23],[31,32,33],[11,12,13]])
* m.permuteRows([1,2,0],true); // m = mat([[11,12,13],[21,22,23],[31,32,33]])
*/
JSMatrix.Matrix.prototype.permuteRows = function(coeffs,backward) {
	var n, r;
	if (coeffs.nElements==undefined) { // coeffs == Array
		n = coeffs.length;
		coeffs = coeffs.slice();
	} else { // coeffs == Vector
		n = coeffs.nElements;
		coeffs = coeffs.toArray();
	}
	if (n != this.nRows) { /*jsmLogger.warn('Matrix.permuteRows: wrong argument size');*/ return; }
	for (r=0; r<n; r++) { if (r<0 || r>=this.nRows) { /*jsmLogger.warn('Matrix.permuteRows: wrong argument values');*/ return; } }
	backward = backward || false;
	var rr;
	if (backward) {
		for (r=0; r<n; r++) {
			rr = coeffs[r];
			while (rr!=r) {
				this.swapRows(r,rr);
				coeffs[r] = coeffs[rr];
				coeffs[rr] = rr;
				rr = coeffs[r];
			}
		}
	} else {
		for (r=n-1; r>0; r--) {
			rr = coeffs[r];
			this.swapRows(r,rr);
			coeffs[coeffs.indexOf(r)] = rr;
		}
	}
	return;
};


/** Returns copy of receiver with rows permutated according to given coefficients
* @param {Array.<number>|JSMatrix.Vector} coeffs coefficients of permutation.
* @param {boolean=} backward =false] if false, receiver is permutated to given coefficients (forward). If true, from given coefficients (backward)
* @returns {JSMatrix.Matrix} copy of receiver with permutated rows
* @example
* m1 = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* m2 = m1.rowPermutation([1,2,0]); // m2 = Matrix([[21,22,23],[31,32,33],[11,12,13]])
* m3 = m1.rowPermutation([1,2,0],true); // m3 = mat([[31,32,33],[11,12,13],[21,22,23]])
*/
JSMatrix.Matrix.prototype.rowPermutation = function(coeffs,backward) {
	var r, c, n = coeffs.nElements==undefined? coeffs.length : coeffs.nElements;
	if (n != this.nRows) { /*jsmLogger.warn('Matrix.rowPermutation: wrong argument size');*/ return null; }
	for (r=0; r<n; r++) { if (r<0 || r>=this.nRows) { /*jsmLogger.warn('Matrix.rowPermutation: wrong argument values');*/ return null; } }
	var ret = new JSMatrix.Matrix();
	if (backward) {
		for (r=0; r<this.nRows; r++) { ret[r] = {}; }
		for (r=0; r<this.nRows; r++) {
			for (c=0; c<this.nCols; c++) {
				ret[coeffs[r]][c] = this[r][c];
			}
		}
	} else {
		for (r=0; r<this.nRows; r++) {
			ret[r] = {};
			for (c=0; c<this.nCols; c++) {
				ret[r][c] = this[coeffs[r]][c];
			}
		}
	}
	ret.nRows = this.nRows;
	ret.nCols = this.nCols;
	return ret;
};

/** Returns receiver multiplied by float f ( ret = this * f )
* @param {number} f float multiplier
* @returns {JSMatrix.Matrix} copy of receiver multiplied by f
* @example
* m1 = JSMatrix.mat([[1,2,3],[4,5,6]]);
* m2 = m1.mulf(3); // m2 = Matrix([[3,6,9],[12,15,18]])
*/
JSMatrix.Matrix.prototype.mulf = function(f) {
	var ret = new JSMatrix.Matrix(), r, c, nr = this.nRows, nc = this.nCols;
	for (r=0; r<nr; r++) {
		ret[r] = {};
		for (c=0; c<nc; c++) {
			ret[r][c] = f*this[r][c];
		}
	}
	ret.nRows = nr;
	ret.nCols = nc;
	return ret;
};

/** Multiply receiver by float f ( this *= f )
* @param {number} f float multiplier
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6]]);
* m.imulf(3); // m = Matrix([[3,6,9],[12,15,18]])
*/
JSMatrix.Matrix.prototype.imulf = function(f) {
	for (var r=0; r<this.nRows; r++) {
		for (var c=0; c<this.nCols; c++) { this[r][c] *= f; }
	}
};

/** Returns product of receiver and given vector ( ret = this * vec )
* @param {JSMatrix.Vector} vec vector to multiply
* @returns {JSMatrix.Vector} copy of receiver multiplied by vec
* @example
* m = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* v1 = JSMatrix.vec([1,2,3]);
* v2 = m.mulv(v1); // v2 = Vector([74,134,194])
*/
JSMatrix.Matrix.prototype.mulv = function(vec) {
	if (!this.canMultiplyVec(vec)) { /*jsmLogger.warn('Matrix.mulv: dimensions mismatched');*/ return null; }
	var ret = new JSMatrix.Vector();
	var temp;
	for (var r=0; r<this.nRows; r++) {
		temp = 0.;
		for (var c=0; c<this.nCols; c++) {
			temp += this[r][c]*vec[c]
		}
		ret[r] = temp;
	}
	ret.nElements = this.nRows
	return ret;
};

/** Returns product of receiver and given matrix ( ret = this * mat )
* @param {JSMatrix.Matrix} mat matrix to multiply
* @returns {JSMatrix.Matrix} copy of receiver multiplied by mat
* @example
* m1 = JSMatrix.mat([[11,12],[21,22]]);
* m2 = JSMatrix.mat([[1,2],[3,4]]);
* m3 = m1.mulm(m2); // m3 = Matrix([[47,70],[87,130]])
*/
JSMatrix.Matrix.prototype.mulm = function(mat) {
	if (!this.canMultiplyMat(mat)) { /*jsmLogger.warn('Matrix.mulm: dimensions mismatched');*/ return null; }
	var ret = new JSMatrix.Matrix(), r, c, k;
	var temp;
	for (r=0; r<this.nRows; r++) {
		ret[r] = {};
		for (c=0; c<mat.nCols; c++) {
			temp = 0.;
			for (k=0; k<this.nCols; k++) {
				temp += this[r][k]*mat[k][c];
			}
			ret[r][c] = temp;
		}
	}
	ret.nRows = this.nRows;
	ret.nCols = mat.nCols;
	return ret;
};

/** Returns product of receiver and given multiplier (Matrix, Vector or float)
* @param {JSMatrix.Matrix|JSMatrix.Vector|number} what matrix, vector or float to multiply
* @returns {JSMatrix.Matrix|JSMatrix.Vector} copy of receiver multiplied by what
* @example
* m1 = JSMatrix.mat([[1,2,3],[4,5,6]]);
* m2 = m1.mul(3); // m2 = Matrix([[3,6,9],[12,15,18]])
* //
* m3 = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* v1 = JSMatrix.vec([1,2,3]);
* v2 = m3.mul(v1) // v2 = Vector([74,134,194]);
* //
* m4 = JSMatrix.mat([[11,12],[21,22]]);
* m5 = JSMatrix.mat([[1,2],[3,4]]);
* m6 = m4.mul(m5); // m6 = Matrix([[47,70],[87,130]])
*/
JSMatrix.Matrix.prototype.mul = function(what) {
	if (what instanceof JSMatrix.Vector) { return this.mulv(what); }
	if (what instanceof JSMatrix.Matrix) { return this.mulm(what); }
	if ((typeof what).toLowerCase() == 'number') { return this.mulf(what); }
	/*jsmLogger.warn('Matrix.mul: wrong argument');*/
	return null;
};

/** Alias for mul(), see {@JSMatrix.Matrix#mul}
* @function
*/
JSMatrix.Matrix.prototype.x = JSMatrix.Matrix.prototype.mul;

/** Modifies receiver to become dyadic product of two vectors (this = vec1*vec2^T ). Receiver's size is adjusted
* @param {JSMatrix.Vector} vec1 1st vector of product
* @param {JSMatrix.Vector} vec2 2nd vector of product
* @example
* m = JSMatrix.mat([[4],[5]]);
* v1 = JSMatrix.vec([1,2,3]);
* v2 = JSMatrix.vec([4,5,6]);
* m.beDyadicProductOf(v1,v2); // m = Matrix([[4,5,6],[8,10,12],[12,15,18]])
*/
JSMatrix.Matrix.prototype.beDyadicProductOf = function(vec1,vec2) {
	var r, c, nr = vec1.nElements, nc = vec2.nElements;
	for (r=0; r<nr; r++) {
		if (r >= this.nRows) { this[r] = {}; }
		for (c=0; c<nc; c++) { this[r][c] = vec1[r]*vec2[c]; }
		for (c=nc; c<this.nCols; c++) { delete this[r][c]; }
	}
	for (r=nr; r<this.nRows; r++) { delete this[r]; }
	this.nRows = nr;
	this.nCols = nc;
};

/** Modifies receiver to become product of two matrices (this = mat1*mat2 ). Receiver's size is adjusted
* @param {JSMatrix.Matrix} mat1 1st matrix of product
* @param {JSMatrix.Matrix} mat2 2nd matrix of product
* @example
* m1 = JSMatrix.mat([[1],[2]]);
* m2 = JSMatrix.mat([[11,12],[21,22]]);
* m3 = JSMatrix.mat([[1,2],[3,4]]);
* m1.beProductOf(m2,m3); // m1 = Matrix([[47,70],[87,130]])
*/
JSMatrix.Matrix.prototype.beProductOf = function(mat1,mat2) {
	if (!mat1.canMultiplyMat(mat2)) { /*jsmLogger.warn('Matrix.beProductOf: dimensions mismatched');*/ return; }
	var r, c, k, temp, nr = mat1.nRows, nc = mat2.nCols;
	for (r=0; r<nr; r++) {
		if (r >= this.nRows) { this[r] = {}; }
		for (c=0; c<nc; c++) {
			temp = 0.;
			for (k=0; k<mat1.nCols; k++) { temp += mat1[r][k]*mat2[k][c]; }
			this[r][c] = temp;
		}
		for (c=nc; c<this.nCols; c++) { delete this[r][c]; }
	}
	for (r=nr; r<this.nRows; r++) { delete this[r]; }
	this.nRows = nr;
	this.nCols = nc;
};

/** Modifies receiver to become product of two matrices (this = mat1^T*mat2 )
* @param {JSMatrix.Matrix} mat1 1st matrix of product
* @param {JSMatrix.Matrix} mat2 2nd matrix of product
* @example
* m1 = JSMatrix.mat([[1],[2]]);
* m2 = JSMatrix.mat([[11,21],[12,22]]);
* m3 = JSMatrix.mat([[1,2],[3,4]]);
* m1.beTProductOf(m2,m3); // m1 = Matrix([[47,70],[87,130]])
*/
JSMatrix.Matrix.prototype.beTProductOf = function(mat1,mat2) {
	var r, c, k, temp, nr = mat1.nCols, nc = mat2.nCols;
	if (mat1.nRows != mat2.nRows) { /*jsmLogger.warn('Matrix.beTProductOf: dimensions mismatched');*/ return; }
	for (r=0; r<nr; r++) {
		if (r >= this.nRows) { this[r] = {}; }
		for (c=0; c<nc; c++) {
			temp = 0.;
			for (k=0; k<mat1.nRows; k++) { temp += mat1[k][r]*mat2[k][c]; }
			this[r][c] = temp;
		}
		for (c=nc; c<this.nCols; c++) { delete this[r][c]; }
	}
	for (r=nr; r<this.nRows; r++) { delete this[r]; }
	this.nRows = nr;
	this.nCols = nc;
};

/** Modifies receiver to become product of two matrices (this = mat1*mat2^T )
* @param {JSMatrix.Matrix} mat1 1st matrix of product
* @param {JSMatrix.Matrix} mat2 2nd matrix of product
* @example
* m1 = JSMatrix.mat([[1],[2]]);
* m2 = JSMatrix.mat([[11,12],[21,22]]);
* m3 = JSMatrix.mat([[1,3],[2,4]]);
* m1.beProductTOf(m2,m3); // m1 = Matrix([[47,70],[87,130]])
*/
JSMatrix.Matrix.prototype.beProductTOf = function(mat1,mat2) {
	var r, c, k, temp, nr = mat1.nRows, nc = mat2.nRows;
	if (mat1.nCols != mat2.nCols) { /*jsmLogger.warn('Matrix.beProductTOf: dimensions mismatched');*/ return; }
	for (r=0; r<nr; r++) {
		if (r >= this.nRows) { this[r] = {}; }
		for (c=0; c<nc; c++) {
			temp = 0.;
			for (k=0; k<mat1.nCols; k++) { temp += mat1[r][k]*mat2[c][k]; }
			this[r][c] = temp;
		}
		for (c=nc; c<this.nCols; c++) { delete this[r][c]; }
	}
	for (r=nr; r<this.nRows; r++) { delete this[r]; }
	this.nRows = nr;
	this.nCols = nc;
};

/** Modifies receiver to become product of two matrices (this = mat1^T*mat2^T )
* @param {JSMatrix.Matrix} mat1 1st matrix of product
* @param {JSMatrix.Matrix} mat2 2nd matrix of product
* @example
* m1 = JSMatrix.mat([[1],[2]]);
* m2 = JSMatrix.mat([[11,21],[12,22]]);
* m3 = JSMatrix.mat([[1,3],[2,4]]);
* m1.beTProductTOf(m2,m3); // m1 = Matrix([[47,70],[87,130]])
*/
JSMatrix.Matrix.prototype.beTProductTOf = function(mat1,mat2) {
	var r, c, k, temp, nr = mat1.nRows, nc = mat2.nCols;
	if (mat1.nCols != mat2.nCols) { /*jsmLogger.warn('Matrix.beTProductTOf: dimensions mismatched');*/ return; }
	for (r=0; r<nr; r++) {
		if (r >= this.nRows) { this[r] = {}; }
		for (c=0; c<nc; c++) {
			temp = 0.;
			for (k=0; k<mat1.nRows; k++) { temp += mat1[k][r]*mat2[c][k]; }
			this[r][c] = temp;
		}
		for (c=nc; c<this.nCols; c++) { delete this[r][c]; }
	}
	for (r=nr; r<this.nRows; r++) { delete this[r]; }
	this.nRows = nr;
	this.nCols = nc;
};

/** Returns transposition of receiver ( ret = this^T )
* @returns {JSMatrix.Matrix} transposed copy of receiver
* @example
* m1 = JSMatrix.mat([[1,2,3],[4,5,6]]);
* m2 = m1.transposed(); // m2 = Matrix([[1,4],[2,5],[3,6]])
*/
JSMatrix.Matrix.prototype.transposed = function() {
	var ret = new JSMatrix.Matrix(this.nCols,this.nRows);
	for (var r=0; r<this.nCols; r++) {
		ret[r] = {};
		for (var c=0; c<this.nRows; c++) {
			ret[r][c] = this[c][r];
		}
	}
	return ret;
};

/** Alias for transposed(), see {@JSMatrix.Matrix#transposed}
* @function
*/
JSMatrix.Matrix.prototype.T = JSMatrix.Matrix.prototype.transposed;

/** Modifies receiver to become transposition of itself ( this = this^T )
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6]]);
* m.transpose(); // m = Matrix([[1,4],[2,5],[3,6]])
* m.transpose(); // m = Matrix([[1,2,3],[4,5,6]])
*/
JSMatrix.Matrix.prototype.transpose = function() {
	var temp, r, c, nr = this.nRows, nc = this.nCols, n=Math.min(nr,nc);
	for (r=0; r<n; r++) {
		for (c=r+1; c<n; c++) {
			temp = this[r][c];
			this[r][c] = this[c][r];
			this[c][r] = temp;
		}
	}
	if (nr>nc) {
		for (r=nc; r<nr; r++) {
			for (c=0; c<nc; c++) {
				this[c][r] = this[r][c];
			}
			delete this[r];
		}
	} else if (nr<nc) {
		for (c=nr; c<nc; c++) {
			this[c] = {};
			for (r=0; r<nr; r++) {
				this[c][r] = this[r][c];
				delete this[r][c];
			}
		}
	}
	this.nRows = nc;
	this.nCols = nr;
};

/** Modifies receiver to become transposition of given matrix (this = mat^T ). Receiver's size is adjusted
* @param {JSMatrix.Matrix} mat matrix of transposition
* @example
* m1 = JSMatrix.mat([[4],[5]]);
* m2 = JSMatrix.mat([[1,2,3],[4,5,6]]);
* m1.beTranspositionOf(m2) // m1 = mat([[1,4],[2,5],[3,6]])
*/
JSMatrix.Matrix.prototype.beTranspositionOf = function(mat) {
	var r, c, nr = mat.nCols, nc = mat.nRows;
	for (r=0; r<nr; r++) {
		if (r >= this.nRows) { this[r] = {}; }
		for (c=0; c<nc; c++) { this[r][c] = mat[c][r]; }
		for (c=nc; c<this.nCols; c++) { delete this[r][c]; }
	}
	for (r=nr; r<this.nRows; r++) { delete this[r]; }
	this.nRows = nr;
	this.nCols = nc;
};

/** Symmetrize receiver ( this = .5*(this+this^T) ). Receiver must be square
* @example
* m = JSMatrix.mat([[1,2,3],[0,-1,5],[-1,1,7]])
* m.symmetrize(); // m = Matrix([[1,1,1],[1,-1,3],[1,3,7]])
* b = m.isSymmetric(); // b = true
*/
JSMatrix.Matrix.prototype.symmetrize = function() {
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.symmetrize: matrix is not square');*/ return; }
	var r, c, temp, n = this.nRows;
	for (r=0; r<n-1; r++) {
		for (c=r+1; c<n; c++) {
			temp = .5*(this[r][c]+this[c][r]);
			this[r][c] = temp;
			this[c][r] = temp;
		}
	}
};

/** Returns symmetric part of receiver ( ret = .5*(this+this^T) ). Receiver must be square
* @example
* m1 = JSMatrix.mat([[1,2,3],[0,-1,5],[-1,1,7]])
* m2 = m1.symmetrized(); // m2 = Matrix([[1,1,1],[1,-1,3],[1,3,7]])
* b = m2.isSymmetric(); // b = true
*/
JSMatrix.Matrix.prototype.symmetrized = function() {
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.symmetrized: matrix is not square');*/ return null; }
	var ret = new JSMatrix.Matrix();
	var r, c, temp, n = this.nRows;
	for (r=0; r<n; r++) { ret[r] = {}; }
	for (r=0; r<n; r++) { ret[r][r] = this[r][r]; }
	for (r=0; r<n-1; r++) {
		for (c=r+1; c<n; c++) {
			temp = .5*(this[r][c]+this[c][r]);
			ret[r][c] = temp;
			ret[c][r] = temp;
		}
	}
	ret.nRows = n;
	ret.nCols = n;
	return ret;
};

/** Alias for symmetrized(), see {@JSMatrix.Matrix#symmetrized}
* @function
*/
JSMatrix.Matrix.prototype.giveSymmetricPart = JSMatrix.Matrix.prototype.symmetrized;

/** Modifies receiver to become symmetric part of given matrix (this = .5*(mat+mat^T) ). mat has to be square
* @param {JSMatrix.Matrix} mat given matrix. Receiver's size is adjusted
* @example
* m1 = JSMatrix.mat([[1],[3]]);
* m2 = JSMatrix.mat([[1,2,3],[0,-1,5],[-1,1,7]])
* m1.beSymmetricPartOf(m2); // m1 = Matrix([[1,1,1],[1,-1,3],[1,3,7]])
* b = m1.isSymmetric(); // b = true
*/
JSMatrix.Matrix.prototype.beSymmetricPartOf = function(mat) {
	if (!mat.isSquare()) { /*jsmLogger.warn('Matrix.beSymmetricPartOf: argument matrix is not square');*/ return; }
	var r, c, temp, n = mat.nRows;
	for (r=0; r<n; r++) {
		if (r>=this.nRows) { this[r] = {}; }
		this[r][r] = mat[r][r];
	}
	for (r=0; r<n-1; r++) {
		for (c=r+1; c<n; c++) {
			temp = .5*(mat[r][c]+mat[c][r]);
			this[r][c] = temp;
			this[c][r] = temp;
		}
		for (c=n; c<this.nCols; c++) { delete this[r][c]; }
	}
	for (r=n; r<this.nRows; r++) { delete this[r]; }
	this.nRows = n;
	this.nCols = n;
};

/** Anti-symmetrize receiver ( this = this - .5*(this+this^T) ). Receiver must be square
* @example
* m = JSMatrix.mat([[1,2,3],[0,-1,5],[-1,1,7]])
* m.antiSymmetrize(); // m = Matrix([[0,1,2],[-1,0,2],[-2,-2,0]])
* b = m.isAntiSymmetric(); // b = true
*/
JSMatrix.Matrix.prototype.antiSymmetrize = function() {
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.antiSymmetrize: matrix is not square');*/ return; }
	var r, c, temp, n = this.nRows;
	for (r=0; r<n; r++) { this[r][r] = 0.; }
	for (r=0; r<n; r++) {
		for (c=r+1; c<n; c++) {
			temp = .5*(this[r][c]+this[c][r]);
			this[r][c] -= temp;
			this[c][r] -= temp;
		}
	}
};

/** Returns anti-symetric part of receiver ( ret = this - .5*(this+this^T) ). Receiver must be square
* @example
* m1 = JSMatrix.mat([[1,2,3],[0,-1,5],[-1,1,7]])
* m2 = m1.antiSymmetrized(); // m2 = Matrix([[0,1,2],[-1,0,3],[-2,0,0]])
* b = m2.isAntiSymmetric(); // b = true
*/
JSMatrix.Matrix.prototype.antiSymmetrized = function() {
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.antiSymmetrized: matrix is not square');*/ return null; }
	var ret = new JSMatrix.Matrix();
	var r, c, temp, n = this.nRows;
	for (r=0; r<n; r++) { ret[r] = {}; }
	for (r=0; r<n; r++) { ret[r][r] = 0.; }
	for (r=0; r<n-1; r++) {
		for (c=r+1; c<n; c++) {
			temp = .5*(this[r][c]+this[c][r]);
			ret[r][c] = this[r][c] - temp;
			ret[c][r] = this[c][r] - temp;
		}
	}
	ret.nRows = n;
	ret.nCols = n;
	return ret;
};

/** Alias for symetrized(), see {@JSMatrix.Matrix#antiSymmetrized}
* @function
*/
JSMatrix.Matrix.prototype.giveAntiSymmetricPart = JSMatrix.Matrix.prototype.antiSymmetrized;

/** Modifies receiver to become anti-symmetric part of given matrix (this = mat - .5*(mat+mat^T) ). mat has to be square
* @param {JSMatrix.Matrix} mat given matrix
* @example
* m1 = JSMatrix.mat([[1],[3]]);
* m2 = JSMatrix.mat([[1,2,3],[0,-1,5],[-1,1,7]])
* m1.beAntiSymmetricPartOf(m2); // m1 = Matrix([[0,1,2],[-1,0,3],[-2,0,0]])
* b = m1.isAntiSymmetric(); // b = true
*/
JSMatrix.Matrix.prototype.beAntiSymmetricPartOf = function(mat) {
	if (!mat.isSquare()) { /*jsmLogger.warn('Matrix.beAntiSymmetricPartOf: argument matrix is not square');*/ return; }
	var r, c, temp, n = mat.nRows;
	for (r=0; r<n; r++) {
		if (r>=this.nRows) { this[r] = {}; }
		this[r][r] = 0.;
	}
	for (r=0; r<n; r++) {
		for (c=r+1; c<n; c++) {
			temp = .5*(mat[r][c]+mat[c][r]);
			this[r][c] = mat[r][c] - temp;
			this[c][r] = mat[c][r] - temp;
		}
		for (c=n; c<this.nCols; c++) { delete this[r][c]; }
	}
	for (r=n; r<this.nRows; r++) { delete this[r]; }
	this.nRows = n;
	this.nCols = n;
};

/** Returns submatrix of receiver ( ret = this[rows][cols] )
* @param {Array.<number>|JSMatrix.Vector} rows array containing rows coefficients of desired submatrix
* @param {Array.<number>|JSMatrix.Vector} cols array containing columns coefficients of desired submatrix
* @returns {JSMatrix.Matrix} desired submatrix
* @example
* m1 = JSMatrix.mat([[11,12,13,14],[21,22,23,24],[31,32,33,34],[41,42,43,44]]);
* m2 = m1.getSubMatrix([1,2],[0,2]); // m2 = Matrix([[21,23],[31,33]])
*/
JSMatrix.Matrix.prototype.getSubMatrix = function(rows,cols) {
	var nr = rows.nElements==undefined? rows.length : rows.nElements;
	var nc = cols.nElements==undefined? cols.length : cols.nElements;
	if (nr>this.nRows || nc>this.nCols) { /*jsmLogger.warn('Matrix.getSubMatrix: wrong argument dimensions');*/ return null; }
	var r,c,rr,cc;
	var ret = new JSMatrix.Matrix();
	for (r=0; r<nr; r++) {
		ret[r] = {};
		for (c=0; c<nc; c++) {
			rr = rows[r]; cc = cols[c];
			if (rr<0 || rr>=this.nRows) { /*jsmLogger.warn('Matrix.getSubMatrix: wrong argument values');*/ return null; }
			if (cc<0 || cc>=this.nCols) { /*jsmLogger.warn('Matrix.getSubMatrix: wrong argument values');*/ return null; }
			ret[r][c] = this[rr][cc];
		}
	}
	ret.nRows = nr;
	ret.nCols = nc;
	return ret;
};

/** Alias for getSubMatrix, see {@JSMatrix.Matrix#getSubMatrix}
* @function
*/
JSMatrix.Matrix.prototype.getsm = JSMatrix.Matrix.prototype.getSubMatrix;

/** Sets receiver to be submatrix of given matrix ( this = mat[rows][cols] )
* @param {JSMatrix.Matrix} mat matrix to get submatrix from
* @param {Array.<number>|JSMatrix.Vector} rows array containing rows coefficients of desired submatrix
* @param {Array.<number>|JSMatrix.Vector} cols array containing columns coefficients of desired submatrix
* @example
* m1 = JSMatrix.mat([[11,12,13,14],[21,22,23,24],[31,32,33,34],[41,42,43,44]]);
* m2 = JSMatrix.mat([[1],[2]]);
* m2.beSubMatrixOf(m1,[1,2],[0,2]); // m2 = Matrix([[21,23],[31,33]])
*/
JSMatrix.Matrix.prototype.beSubMatrixOf = function(mat,rows,cols) {
	var nr = rows.nElements==undefined? rows.length : rows.nElements;
	var nc = cols.nElements==undefined? cols.length : cols.nElements;
	if (nr>mat.nRows || nc>mat.nCols) { /*jsmLogger.warn('Matrix.beSubMatrixOf: wrong argument dimensions');*/ return; }
	var r,c;
	for (r=0; r<nr; r++) { if (rows[r]<0 || rows[r]>=mat.nRows) { /*jsmLogger.warn('Matrix.beSubMatrixOf: wrong argument values');*/ return; } }
	for (c=0; c<nc; c++) { if (cols[c]<0 || cols[c]>=mat.nCols) { /*jsmLogger.warn('Matrix.beSubMatrixOf: wrong argument values');*/ return; } }
	for (r=0; r<nr; r++) {
		this[r] = {};
		for (c=0; c<nc; c++) { this[r][c] = mat[rows[r]][cols[c]]; }
		for (c=nc; c<this.nCols; c++) { delete this[r][c]; }
	}
	for (r=nr; r<this.nRows; r++) { delete this[r]; }
	this.nRows = nr;
	this.nCols = nc;
};

/** Sets submatrix at defined positions of receiver ( this[rows][cols] = mat )
* @param {Array.<number>|JSMatrix.Vector} rows array containing rows coefficients to be set
* @param {Array.<number>|JSMatrix.Vector} cols array containing columns coefficients to be set
* @param {JSMatrix.Matrix} mat matrix to be set on desired positions
* @example
* m = JSMatrix.mat([[11,12,13,14],[21,22,23,24],[31,32,33,34],[41,42,43,44]]);
* m.setSubMatrix([1,2],[0,2],JSMatrix.mat([[66,77],[88,99]]));
* // m = Matrix([[11,12,13,14],[66,22,77,24],[88,32,99,34],[41,42,43,44]])
*/
JSMatrix.Matrix.prototype.setSubMatrix = function(rows,cols,mat) {
	var nr = rows.nElements==undefined? rows.length : rows.nElements;
	var nc = cols.nElements==undefined? cols.length : cols.nElements;
	var mnr = mat.nRows, mnc=mat.nCols;
	if (mnr != nr || mnc != nc) { /*jsmLogger.warn('Matrix.setSubMatrix: wrong argument dimensions');*/ return; }
	if (nr>this.nRows || nc>this.nCols) { /*jsmLogger.warn('Matrix.setSubMatrix: wrong argument dimensions');*/ return; }
	var r,c;
	for (r=0; r<nr; r++) { if (rows[r]<0 || rows[r]>=this.nRows) { /*jsmLogger.warn('Matrix.setSubMatrix: wrong argument values');*/ return; } }
	for (c=0; c<nc; c++) { if (cols[c]<0 || cols[c]>=this.nCols) { /*jsmLogger.warn('Matrix.setSubMatrix: wrong argument values');*/ return; } }
	for (r=0; r<nr; r++) {
		for (c=0; c<nc; c++) { this[rows[r]][cols[c]] = mat[r][c]; }
	}
};

/** Alias for setSubMatrix, see {@JSMatrix.Matrix#setSubMatrix}
* @function
*/
JSMatrix.Matrix.prototype.setsm = JSMatrix.Matrix.prototype.setSubMatrix;

/** Increment submatrix at defined positions of receiver ( this[rows][cols] += mat )
* @param {Array.<number>|JSMatrix.Vector} rows array containing rows coefficients to be incrementes
* @param {Array.<number>|JSMatrix.Vector} cols array containing columns coefficients to be incremented
* @param {JSMatrix.Matrix} mat matrix to be incremented on desired positions
* @example
* m = JSMatrix.mat([[11,12,13,14],[21,22,23,24],[31,32,33,34],[41,42,43,44]]);
* m.incrSubMatrix([1,2],[0,2],JSMatrix.mat([[66,77],[88,99]]));
* // m = Matrix([11,12,13,14],[87,22,100,24],[119,32,132,34],[41,42,43,44]])
*/
JSMatrix.Matrix.prototype.incrSubMatrix = function(rows,cols,mat) {
	var nr = rows.nElements==undefined? rows.length : rows.nElements;
	var nc = cols.nElements==undefined? cols.length : cols.nElements;
	var mnr = mat.nRows, mnc=mat.nCols;
	if (mnr != nr || mnc != nc) { /*jsmLogger.warn('Matrix.incrSubMatrix: wrong argument dimensions');*/ return; }
	if (nr>this.nRows || nc>this.nCols) { /*jsmLogger.warn('Matrix.incrSubMatrix: wrong argument dimensions');*/ return; }
	var r,c;
	for (r=0; r<nr; r++) { if (rows[r]<0 || rows[r]>=this.nRows) { /*jsmLogger.warn('Matrix.incrSubMatrix: wrong argument values');*/ return; } }
	for (c=0; c<nc; c++) { if (cols[c]<0 || cols[c]>=this.nCols) { /*jsmLogger.warn('Matrix.incrSubMatrix: wrong argument values');*/ return; } }
	for (r=0; r<nr; r++) {
		for (c=0; c<nc; c++) { this[rows[r]][cols[c]] += mat[r][c]; }
	}
};

/** Alias for incrSubMatrix, see {@JSMatrix.Matrix#incrSubMatrix}
* @function
*/
JSMatrix.Matrix.prototype.incrsm = JSMatrix.Matrix.prototype.incrSubMatrix;

/** Decrement submatrix at defined positions of receiver ( this[rows][cols] -= mat )
* @param {Array.<number>|JSMatrix.Vector} rows array containing rows coefficients to be decrementes
* @param {Array.<number>|JSMatrix.Vector} cols array containing columns coefficients to be decremented
* @param {JSMatrix.Matrix} mat matrix to be decremented on desired positions
* @example
* m = JSMatrix.mat([[11,12,13,14],[21,22,23,24],[31,32,33,34],[41,42,43,44]]);
* m.decrSubMatrix([1,2],[0,2],JSMatrix.mat([[66,77],[88,99]]));
* // m = Matrix([11,12,13,14],[-45,22,-54,24],[-57,32,-66,34],[41,42,43,44]])
*/
JSMatrix.Matrix.prototype.decrSubMatrix = function(rows,cols,mat) {
	var nr = rows.nElements==undefined? rows.length : rows.nElements;
	var nc = cols.nElements==undefined? cols.length : cols.nElements;
	var mnr = mat.nRows, mnc=mat.nCols;
	if (mnr != nr || mnc != nc) { /*jsmLogger.warn('Matrix.decrSubMatrix: wrong argument dimensions');*/ return; }
	if (nr>this.nRows || nc>this.nCols) { /*jsmLogger.warn('Matrix.decrSubMatrix: wrong argument dimensions');*/ return; }
	var r,c;
	for (r=0; r<nr; r++) { if (rows[r]<0 || rows[r]>=this.nRows) { /*jsmLogger.warn('Matrix.decrSubMatrix: wrong argument values');*/ return; } }
	for (c=0; c<nc; c++) { if (cols[c]<0 || cols[c]>=this.nCols) { /*jsmLogger.warn('Matrix.decrSubMatrix: wrong argument values');*/ return; } }
	for (r=0; r<nr; r++) {
		for (c=0; c<nc; c++) { this[rows[r]][cols[c]] -= mat[r][c]; }
	}
};

/** Alias for decrSubMatrix, see {@JSMatrix.Matrix#decrSubMatrix}
* @function
*/
JSMatrix.Matrix.prototype.decrsm = JSMatrix.Matrix.prototype.decrSubMatrix;

/** assemble receiver to given matrix ( mat[rows][cols] += this )
* @param {JSMatrix.Matrix} mat matrix where receiver is assembled
* @param {Array.<number>|JSMatrix.Vector} rows array containing row coefficients of mat to be incremented
* @param {Array.<number>|JSMatrix.Vector} cols array containing column coefficients of mat to be incremented
* @example
* m1 = JSMatrix.mat([[11,12,13,14],[21,22,23,24],[31,32,33,34],[41,42,43,44]]);
* m2 = JSMatrix.mat([[66,77],[88,99]]);
* m2.assemble(m1,[1,2],[0,2]);
* // m1 = Matrix([11,12,13,14],[87,22,100,24],[119,32,132,34],[41,42,43,44]])
*/
JSMatrix.Matrix.prototype.assemble = function(mat,rows,cols) {
	var nr = rows.nElements==undefined? rows.length : rows.nElements;
	var nc = cols.nElements==undefined? cols.length : cols.nElements;
	var mnr = this.nRows, mnc = this.nCols;
	if (mnr != nr || mnc != nc) { /*jsmLogger.warn('Matrix.assemble: wrong argument dimensions');*/ return; }
	if (nr>mat.nRows || nc>mat.nCols) { /*jsmLogger.warn('Matrix.assemble: wrong argument dimensions');*/ return; }
	var r,c;
	for (r=0; r<nr; r++) { if (rows[r]<0 || rows[r]>=mat.nRows) { /*jsmLogger.warn('Matrix.assemble: wrong argument values');*/ return; } }
	for (c=0; c<nc; c++) { if (cols[c]<0 || cols[c]>=mat.nCols) { /*jsmLogger.warn('Matrix.assemble: wrong argument values');*/ return; } }
	for (r=0; r<nr; r++) {
		for (c=0; c<nc; c++) { mat[rows[r]][cols[c]] += this[r][c]; }
	}
};

/** Returns copy of receiver
* @returns {JSMatrix.Matrix} copy of receiver
* @example
* m1 = JSMatrix.mat([[11,12],[21,22]]);
* m2 = m1;
* m3 = m1.copied();
* m1.set(1,0,6);
* // m1 = Matrix([[11,12],[6,22]])
* // m2 = Matrix([[11,12],[6,22]])
* // m3 = Matrix([[11,12],[21,22]])
*/
JSMatrix.Matrix.prototype.copied = function() {
	var ret = new JSMatrix.Matrix();
	for (var r=0; r<this.nRows; r++) {
		ret[r] = {};
		for (var c=0; c<this.nCols; c++) {
			ret[r][c] = this[r][c];
		}
	}
	ret.nRows = this.nRows;
	ret.nCols = this.nCols;
	return ret;
};

/** Alias for copied(), see {@JSMatrix.Matrix#copied}
* @function
*/
JSMatrix.Matrix.prototype.copy = JSMatrix.Matrix.prototype.copied;

/** Modifies receiver to become copy of given matrix (this = mat). Receiver's size is adjusted
* @param {JSMatrix.Matrix} mat matrix to be copied to receiver
* @example
* m1 = JSMatrix.mat([[4],[5]]);
* m2 = JSMatrix.mat([[1,2,3],[4,5,6]]);
* m1.beCopyOf(m2); // m1 = Matrix([[1,2,3],[4,5,6]])
*/
JSMatrix.Matrix.prototype.beCopyOf = function(mat) {
	var r, c, nr = mat.nRows, nc = mat.nCols;
	for (r=0; r<nr; r++) {
		if (r >= this.nRows) { this[r] = {}; }
		for (c=0; c<nc; c++) { this[r][c] = mat[r][c]; }
		for (c=nc; c<this.nCols; c++) { delete this[r][c]; }
	}
	for (r=nr; r<this.nRows; r++) { delete this[r]; }
	this.nRows = nr;
	this.nCols = nc;
};

/** Returns vector containing receiver's diagonal elements
* @returns {JSMatrix.Vector} vector containing receiver's diagonal elements
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* v = m.diag(); // v = Vector([1,5,9])
*/
JSMatrix.Matrix.prototype.diagonalToVector = function() {
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.diagonalToVector: matrix is not square');*/ return null; }
	var ret = new JSMatrix.Vector()
	for (var r=0; r<this.nRows; r++) { ret[r] = this[r][r]; }
	ret.nElements = this.nRows;
	return ret;
};

/** Alias for diagonalToVector(), see {@JSMatrix.Matrix#diagonalToVector}
* @function
*/
JSMatrix.Matrix.prototype.diag = JSMatrix.Matrix.prototype.diagonalToVector;

/** Returns string representation of receiver
* @returns {string} string representation of receiver
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* s = m.toString(); // s = 'Matrix([[1,2,3],[4,5,6],[7,8,9]])'
*/
JSMatrix.Matrix.prototype.toString = function() {
	var nr = this.nRows;
	var nc = this.nCols;
	var r,c,val;
	var ret = "Matrix([[ ";
	if (nr == 0) { return ret+"]])"; }
	if (nc == 0) {
		for (r=0; r<nr-1; r++) {
			ret += " [], ";
		}
		return ret + "])";
	}
	for (r=0; r<nr-1; r++) {
		for (c=0; c<nc-1; c++) {
			val = this[r][c];
			ret += val.toExponential(3)+", ";
		}
		val = this[r][nc-1];
		ret += val.toExponential(3)+" ], [ ";
	}
	for (c=0; c<nc-1; c++) {
		val = this[nr-1][c];
		ret += val.toExponential(3)+", ";
	}
	val = this[nr-1][nc-1];
	ret += val.toExponential(3)+" ]])"
	return ret;
};

/** Resize receiver according to given size (delete extra elements or add zero elements)
* @param {number} nRows new number of rows
* @param {number} nCols new number of columns
* @example
* m1 = JSMatrix.mat([[11,12,13,14],[21,22,23,24],[31,32,33,34],[41,42,43,44]]);
* m2 = JSMatrix.mat([[11,12],[21,22]]);
* m1.resize(3,3); // m1 = Matrix([[11,12,13],[21,22,23],[31,32,33]])
* m2.resize(3,3); // m2 = Matrix([[11,12,0],[21,22,0],[0,0,0]])
*/
JSMatrix.Matrix.prototype.resize = function(nRows,nCols) {
	if (this.nRows < nRows) {
		for (var r=this.nRows; r<nRows; r++) {
			this[r] = {};
			for (var c=0; c<this.nCols; c++) {
				this[r][c] = 0.;
			}
		}
		this.nRows = nRows;
	}
	if (this.nRows > nRows) {
		for (var r=nRows; r<this.nRows; r++) { delete this[r]; }
		this.nRows = nRows;
	}
	if (this.nCols < nCols) {
		for (var r=0; r<this.nRows; r++) {
			for (var c=this.nCols; c<nCols; c++) {
				this[r][c] = 0.;
			}
		}
		this.nCols = nCols;
	}
	if (this.nCols > nCols) {
		for (var r=0; r<this.nRows; r++) {
			for (var c=nCols; c<this.nCols; c++) {
				delete this[r][c];
			}
		}
		this.nCols = nCols;
	}
};

/** Returns resized copy of receiver according to given size (delete extra elements or add zero elements)
* @param {number} nRows new number of rows
* @param {number} nCols new number of columns
* @returns {JSMatrix.Matrix} resized copy of receiver
* @example
* m1 = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* m2 = m1.resized(2,4); // m2 = Matrix([[11,12,13,0],[21,22,23,0]])
* m3 = m1.resized(4,2); // m3 = Matrix([[11,12],[21,22],[31,32],[0,0]])
*/
JSMatrix.Matrix.prototype.resized = function(nRows,nCols) {
	var r, c, ret = new JSMatrix.Matrix();
	for (r=0; r<nRows; r++) { ret[r] = {}; }
	if (nRows<=this.nRows) {
		if (nCols<=this.nCols) {
			for (r=0; r<nRows; r++) {
				for (c=0; c<nCols; c++) { ret[r][c] = this[r][c]; }
			}
		} else {
			for (r=0; r<nRows; r++) {
				for (c=0; c<this.nCols; c++) { ret[r][c] = this[r][c]; }
				for (c=this.nCols; c<nCols; c++) { ret[r][c] = 0.; }
			}
		}
	} else {
		if (nCols<=this.nCols) {
			for (r=0; r<this.nRows; r++) {
				for (c=0; c<nCols; c++) { ret[r][c] = this[r][c]; }
			}
			for (r=this.nRows; r<nRows; r++) {
				for (c=0; c<nCols; c++) { ret[r][c] = 0.; }
			}
		} else {
			for (r=0; r<this.nRows; r++) {
				for (c=0; c<this.nRows; c++) { ret[r][c] = this[r][c]; }
				for (c=this.nCols; c<nCols; c++) { ret[r][c] = 0.; }
			}
			for (r=this.nRows; r<nRows; r++) {
				for (c=0; c<nCols; c++) { ret[r][c] = 0.; }
			}
		}
	}
	ret.nRows = nRows;
	ret.nCols = nCols;
	return ret;
};

/** Appends vector(s)/matrix(s) as row(s) to receiver
* @param {JSMatrix.Vector|JSMatrix.Matrix|Array.<JSMatrix.Vector>|Array.<JSMatrix.Matrix>} rows new row(s) to be appended
* @example
* m = JSMatrix.mat([[1,2],[3,4]]);
* m.appendRows(JSMatrix.vec([7,6]));
* // m = Matrix([[1,2],[3,4],[7,6]])
* m.appendRows(JSMatrix.mat([[1,3],[2,4]]));
* // m = Matrix([[1,2],[3,4],[7,6],[1,3],[2,4]])
* m = JSMatrix.mat([[1,2],[3,4]]);
* m.appendRows([JSMatrix.vec([2,1]),JSMatrix.vec([33,44])]);
* // m = Matrix([[1,2],[3,4],[2,1],[33,44]])
* m = JSMatrix.mat([[1,2],[3,4]]);
* m.appendRows([JSMatrix.mat([[11,12],[21,22]]),JSMatrix.mat([[5,6],[7,8]])]);
* // m = Matrix([[1,2],[3,4],[11,12],[21,22],[5,6],[7,8]])
*/
JSMatrix.Matrix.prototype.appendRows = function(rows) {
	if (rows instanceof JSMatrix.Vector) {
		if (rows.nElements != this.nCols) { /*jsmLogger.warn('Matrix.appendRows: wrong argument dimensions');*/ return; }
		this[this.nRows] = {};
		for (var c=0; c<this.nCols; c++) { this[this.nRows][c] = rows[c]; }
		this.nRows++;
		return;
	}
	if (rows instanceof JSMatrix.Matrix) {
		if (rows.nCols != this.nCols) { /*jsmLogger.warn('Matrix.appendRows: wrong argument dimensions');*/ return; }
		for (var r=0; r<rows.nRows; r++) {
			this[this.nRows] = {};
			for (var c=0; c<this.nCols; c++) { this[this.nRows][c] = rows[r][c]; }
			this.nRows++;
		}
		return;
	}
	if (rows instanceof Array) {
		for (var i=0; i<rows.length; i++) {
			this.appendRows(rows[i]);
		}
	}
};

/** Appends vector(s)/matrix(s) as column(s) to receiver
* @param {JSMatrix.Vector|JSMatrix.Matrix|Array.<JSMatrix.Vector>|Array.<JSMatrix.Matrix>} cols new row(s) to be appended
* @example
* m = JSMatrix.mat([[1,2],[3,4]]);
* m.appendCols(JSMatrix.vec([7,6]));
* // m = Matrix([[1,2,7],[3,4,6]])
* m.appendCols(JSMatrix.mat([[1,3],[2,4]]));
* // m = Matrix([[1,2,7,1,3],[3,4,6,2,4]])
* m = JSMatrix.mat([[1,2],[3,4]]);
* m.appendCols([JSMatrix.vec([2,1]),JSMatrix.vec([33,44])]);
* // m = Matrix([[1,2,2,33],[3,4,1,44]])
* m = JSMatrix.mat([[1,2],[3,4]]);
* m.appendCols([JSMatrix.mat([[11,12],[21,22]]),JSMatrix.mat([[9,8],[7,6]])]);
* // m = Matrix([[1,2,11,12,9,8],[3,4,21,22,7,6]])
*/
JSMatrix.Matrix.prototype.appendCols = function(cols) {
	if (cols instanceof JSMatrix.Vector) {
		if (cols.nElements != this.nRows) { /*jsmLogger.warn('Matrix.appendCols: wrong argument dimensions');*/ return; }
		for (var r=0; r<this.nRows; r++) { this[r][this.nCols] = cols[r]; }
		this.nCols++;
		return;
	}
	if (cols instanceof JSMatrix.Matrix) {
		if (cols.nRows != this.nRows) { /*jsmLogger.warn('Matrix.appendCols: wrong argument dimensions');*/ return; }
		for (var r=0; r<this.nRows; r++) {
			for (var c=0; c<cols.nCols; c++) { this[r][this.nCols+c] = cols[r][c]; }
		}
		this.nCols += cols.nCols;
		return;
	}
	if (cols instanceof Array) {
		for (var i=0; i<cols.length; i++) {
			this.appendCols(cols[i]);
		}
	}
};

/** Returns receiver's maximal element
* @returns {number} receiver's maximal element
* @example
* a = JSMatrix.mat([[1,-2,3],[-4,5,-6],[-7,8,-9]]);
* b = a.max(); // b = 8
*/
JSMatrix.Matrix.prototype.max = function() {
	var ret = -1e20;
	for (var r=0; r<this.nRows; r++) {
		for (var c=0; c<this.nCols; c++) {
			ret = Math.max(ret,this[r][c]);
		}
	}
	return ret;
};

/** Returns receiver's minimal element
* @returns {number} receiver's minimal element
* @example
* a = JSMatrix.mat([[1,-2,3],[-4,5,-6],[-7,8,-9]]);
* b = a.min(); // b = -9
*/
JSMatrix.Matrix.prototype.min = function() {
	var ret = +1e20;
	for (var r=0; r<this.nRows; r++) {
		for (var c=0; c<this.nCols; c++) {
			ret = Math.min(ret,this[r][c]);
		}
	}
	return ret;
};

/** Returns receiver's elementwse absolute value
* @returns {JSMatrix.Matrix} receiver's elementwise absolute value
* @example
* a = JSMatrix.mat([[1,-2,3],[-4,5,-6],[-7,8,-9]]);
* b = a.abs(); // b = Matrix([[1,2,3],[4,5,6],[7,8,9]])
*/
JSMatrix.Matrix.prototype.abs = function() {
	var ret = new JSMatrix.Matrix(this.nRows,this.nCols);
	for (var r=0; r<this.nRows; r++) {
		for (var c=0; c<this.nCols; c++) {
			ret[r][c] = Math.abs(this[r][c]);
		}
	}
	return ret;
};

/** Testing matrix equality
* @param {JSMatrix.Matrix} mat matrix to be compared with receiver
* @returns {boolean} true if this==mat, false otherwise
* @example
* a = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* b = JSMatrix.mat([[1,2,2.999999999],[4.0000000000002,5,6],[7,8,8.9999999999]]);
* c = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9],[10,11,12]]);
* t1 = a.isEqualTo(b); // t1 = true
* t2 = a.isEqualTo(c); // t2 = false
*/
JSMatrix.Matrix.prototype.isEqualTo = function(mat) {
	if (!this.isSameSizeAs(mat)) { return false; }
	for (var r=0; r<this.nRows; r++) {
		for (var c=0; c<this.nCols; c++) {
			if (Math.abs(this[r][c] - mat[r][c]) > JSMatrix.TOL) { return false; }
		}
	}
	return true;
};

/** Testing matrix size equality
* @param {JSMatrix.Matrix} mat matrix to be tested
* @returns {boolean} true if mat and receiver has same size, false otherwise
* @example
* a = JSMatrix.mat([[1,2,3],[4,5,6]]);
* b = JSMatrix.mat([[5,6,7],[8,9,10],[11,12,13]]);
* c = JSMatrix.mat([[14,15,16],[17,18,19]]);
* t1 = a.isSameSizeAs(b); // t1 = false
* t2 = a.isSameSizeAs(c); // t2 = true
*/
JSMatrix.Matrix.prototype.isSameSizeAs = function(mat) {
	return (this.nRows==mat.nRows && this.nCols==mat.nCols);
};

/** Testing if receiver can multiply given matrix ( this * mat )
* @param {JSMatrix.Matrix} mat matrix to be tested
* @returns {boolean} true if the multiplication is possible, false otherwise
* @example
* m1 = JSMatrix.mat([[1,2,3],[4,5,6]]);
* m2 = JSMatrix.mat([[7,8,9],[10,11,12]]);
* m3 = JSMatrix.mat([[11,12],[21,22],[31,32]]);
* t1 = m1.canMultiplyMat(m2); // t1 = false
* t2 = m1.canMultiplyMat(m3); // t2 = true
*/
JSMatrix.Matrix.prototype.canMultiplyMat = function(mat) {
	return (this.nCols==mat.nRows);
};

/** Testing if receiver can multiply given vector ( this * vec )
* @param {JSMatrix.Vector} vec vector to be tested
* @returns {boolean} true if the multiplication is possible, false otherwise
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6]]);
* v1 = JSMatrix.vec([1,2,3]);
* v2 = JSMatrix.vec([4,5]);
* t1 = m.canMultiplyVec(v1); // t1 = true
* t2 = m.canMultiplyVec(v2); // t2 = false
*/
JSMatrix.Matrix.prototype.canMultiplyVec = function(vec) {
	return (this.nCols==vec.nElements);
};

/** Testing if receiver is square
* @returns {boolean} true if receiver is square matrix, false otherwise
* @example
* m1 = JSMatrix.mat([[1,2,3],[4,5,6]]);
* m2 = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* t1 = m1.isSquare(); // t1 = false
* t2 = m2.isSquare(); // t2 = true
*/
JSMatrix.Matrix.prototype.isSquare = function() {
	return (this.nRows==this.nCols);
};

/** Testing if receiver is symmetric
* @returns {boolean} true if receiver is symmetric, false otherwise
* @example
* m1 = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* m2 = JSMatrix.mat([[11,12,13],[12,22,23],[13,23,33]]);
* t1 = m1.isSymmetric(); // t1 = false
* t2 = m2.isSymmetric(); // t2 = true
*/
JSMatrix.Matrix.prototype.isSymmetric = function() {
	if (!this.isSquare()) { return false; }
	for (var r=0; r<this.nRows; r++) {
		for (var c=r+1; c<this.nCols; c++) {
			if (Math.abs(this[r][c] - this[c][r]) > JSMatrix.TOL) { return false; }
		}
	}
	return true;
};

/** Testing if receiver is anti-symmetric
* @returns {boolean} true if receiver is anti-symmetric, false otherwise
* @example
* m1 = JSMatrix.mat([[11,12,13],[21,22,23],[31,32,33]]);
* m2 = JSMatrix.mat([[0,12,13],[-12,0,23],[-13,-23,0]]);
* t1 = m1.isAntiSymmetric(); // t1 = false
* t2 = m2.isAntiSymmetric(); // t2 = true
*/
JSMatrix.Matrix.prototype.isAntiSymmetric = function() {
	if (!this.isSquare()) { return false; }
	for (var r=0; r<this.nRows; r++) {
		for (var c=r; c<this.nCols; c++) {
			if (Math.abs(this[r][c] + this[c][r]) > JSMatrix.TOL) { return false; }
		}
	}
	return true;
};

/** Testing if receiver is transposition of given matrix
* @param {JSMatrix.Matrix} mat given matrix
* @returns {boolean} true if this==mat^T, false otherwise
* @example
* m1 = JSMatrix.mat([[1,2,3],[4,5,6]]);
* m2 = JSMatrix.mat([[1,4],[2,5],[3,6]]);
* m3 = JSMatrix.mat([[1,2],[4,5]]);
* b1 = m1.isTranspositionOf(m2); // b1 = true
* b2 = m1.isTranspositionOf(m3); // b2 = false
*/
JSMatrix.Matrix.prototype.isTranspositionOf = function(mat) {
	if (this.nRows!=mat.nCols || this.nCols!=mat.nRows) { return false; }
	var r,c;
	for (r=0; r<this.nRows; r++) {
		for (c=0; c<this.nCols; c++) {
			if (Math.abs(this[r][c]-mat[c][r]) > JSMatrix.TOL) { return false; }
		}
	}
	return true;
};

/** Testing if receiver is lower triangular matrix
* @returns {boolean} true if receiver is lower triangular, false otherwise
* @example
* m1 = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* m2 = JSMatrix.mat([[1,0,0],[4,5,0],[7,8,9]]);
* t1 = m1.isLowerTriangular(); // t1 = false
* t2 = m2.isLowerTriangular(); // t2 = true
*/
JSMatrix.Matrix.prototype.isLowerTriangular = function() {
	for (var r=0; r<this.nRows; r++) {
		for (var c=r+1; c<this.nCols; c++) {
			if (Math.abs(this[r][c]) > JSMatrix.TOL) { return false; }
		}
	}
	return true;
};

/** Testing if receiver is upper triangular matrix
* @returns {boolean} true if receiver is upper triangular, false otherwise
* @example
* m1 = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* m2 = JSMatrix.mat([[1,2,3],[0,5,6],[0,0,9]]);
* t1 = m1.isUpperTriangular(); // t1 = false
* t2 = m2.isUpperTriangular(); // t2 = true
*/
JSMatrix.Matrix.prototype.isUpperTriangular = function() {
	for (var r=0; r<this.nRows; r++) {
		for (var c=0; c<r; c++) {
			if (Math.abs(this[r][c]) > JSMatrix.TOL) { return false; }
		}
	}
	return true;
};

/** Test receiver's singularity
* @returns {boolean} true if receiver is singular (if abs(this.det()) < JSMatrix.TOL), false otherwise
* @example
* a = JSMatrix.mat([[1,2,3],[2,4,6],[4,5,6]]);
* b = JSMatrix.mat([[4,2,1],[2,5,3],[1,3,6]]);
* s1 = a.isSingular(); // s1 = true
* s2 = b.isSingular(); // s2 = false
*/
JSMatrix.Matrix.prototype.isSingular = function() {
	var d = this.determinant();
	// matrix is singular <=> det(matrix)==0
	return isNaN(d)? true : (Math.abs(d) < JSMatrix.TOL);
};

/** Test if receiver is identity matrix
* @returns {boolean} true if receiver is identity, false otherwise
* @example
* m1 = JSMatrix.mat([[1,0,0],[0,1,0],[0,0,-1]]);
* m2 = JSMatrix.mat([[1,0,0],[0,1,0],[0,0,1]]);
* i1 = m1.isIdentity() // i1 = false
* i2 = m2.isIdentity() // i2 = true
*/
JSMatrix.Matrix.prototype.isIdentity = function() {
	if (!this.isSquare()) { return false; }
	var r, c;
	for (r=0; r<this.nRows; r++) {
		if (Math.abs(this[r][r]-1.) > JSMatrix.TOL) { return false; }
	}
	for (r=1; r<this.nRows; r++) {
		for (c=r+1; c<this.nCols; c++) {
			if (Math.abs(this[r][c]) > JSMatrix.TOL) { return false; }
			if (Math.abs(this[c][r]) > JSMatrix.TOL) { return false; }
		}
	}
	return true;
};

/** Alias for isIdentity, see {@JSMatrix.Matrix#isIdentity}
* @function
*/
JSMatrix.Matrix.prototype.isUnitMatrix = JSMatrix.Matrix.prototype.isIdentity;

/** Test if receiver is zero matrix (full of zeroes)
* @returns {boolean} true if receiver is zero, false otherwise
* @example
* m1 = JSMatrix.mat([[0,0,0],[0,1e-16,0],[0,0,0]]);
* m2 = JSMatrix.mat([[1,0,0],[0,1,0],[0,0,1]]);
* i1 = m1.isZero() // i1 = true
* i2 = m2.isZero() // i2 = false
*/
JSMatrix.Matrix.prototype.isZero = function() {
	for (var i=0; i<this.nRows; i++) {
		for (var j=0; j<this.nCols; j++) {
			if (Math.abs(this[i][j]) > JSMatrix.TOL) { return false; }
		}
	}
	return true;
};

/** Alias for isZero(), see {@JSMatrix.Matrix#isZero}
* @function
*/
JSMatrix.Matrix.prototype.containsOnlyZeros = JSMatrix.Matrix.prototype.isZero;

/** Test if receiver is matrix full of ones
* @returns {boolean} true if receiver is full of ones, false otherwise
* @example
* m1 = JSMatrix.mat([[1,1,1],[1,.999999999999,1],[1,1,1]]);
* m2 = JSMatrix.mat([[1,0,0],[0,1,0],[0,0,1]]);
* i1 = m1.isOne() // i1 = true
* i2 = m2.isOne() // i2 = false
*/
JSMatrix.Matrix.prototype.isOne = function() {
	for (var r=0; r<this.nRows; r++) {
		for (var c=0; c<this.nCols; c++) {
			if (Math.abs(this[r][c]-1.) > 2*JSMatrix.TOL) { return false; }
		}
	}
	return true;
};

/** Alias for isOne(), see {@JSMatrix.Matrix#isOne}
* @function
*/
JSMatrix.Matrix.prototype.containsOnlyOnes = JSMatrix.Matrix.prototype.isOne;

/** Tests receiver's orthogonality (i.e. this*this^T = I => this^-1 = this^T)
* @returns {boolean} true if receiver is orthogonal, false otherwise
* @example
* m1 = JSMatrix.Matrix.Householder(JSMatrix.vec([0,1,0]))
* // m1 = Matrix([[1,0,0],[0,-1,0],[0,0,1]]);
* m2 = JSMatrix.mat([[0,0,1],[1,0,0],[0,1,0]])
* m3 = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]])
* m4 = JSMatrix.mat([[0.8,0,0.6],[0,1,0],[-0.8,0,0.6]]);
* b1 = m1.isOrthogonal() // b1 = true
* b2 = m2.isOrthogonal() // b2 = true
* b3 = m3.isOrthogonal() // b3 = false
* b4 = m4.isOrthogonal() // b4 = true
*/
JSMatrix.Matrix.prototype.isOrthogonal = function() {
	var temp = new JSMatrix.Matrix();
	temp.beProductTOf(this,this);
	return temp.isIdentity();
};

/** Tests receiver's positive definitness
* @returns {boolean} true if receiver is positive definite (if Cholesky decomposition can be performed), false otherwise
* @example
* m1 = JSMatrix.mat([[9,2,1],[2,8,3],[1,3,7]])
* m2 = JSMatrix.mat([[-9,2,1],[2,8,3],[1,3,7]])
* b1 = m1.isPositiveDefinite() // b1 = true
* b2 = m2.isPositiveDefinite() // b2 = false
*/
JSMatrix.Matrix.prototype.isPositiveDefinite = function() {
	// matrix is positive definite <=> cholesky decomposition can be performed
	return this.choleskyDecomposition() != null;
};

/** Tests receiver's involutariness (if this*this=I -> this^-1=this)
* @returns {boolean} true if receiver is involutary, false otherwise
* @example
* m = JSMatrix.mat([[1,0,0],[0,-1,0],[0,0,-1]])
* t = m.isInvolutary() // t = true
*/
JSMatrix.Matrix.prototype.isInvolutary = function() {
	var temp = new JSMatrix.Matrix();
	temp.beProductOf(this,this);
	return temp.isIdentity();
};

/** Returns row permutation matrix using implicit partial pivoting (for each column find row with relatively highest element at that column and place that line to the position such that the highest element is on diagonal. Lines that are placed in this way are not considered for further steps)
* @param {boolean=} saveOrig =true] if true, returns permutated copy of receiver. If false, permutate receiver and return this
* @returns {{mat:JSMatrix.Matrix,coeffs:Array.<number>}|null} {mat:JSMatrix.Matrix,coeffs:Array.<number>} permutated matrix and array containing coefficients of original matrix according to position in the permutated one ( such that this.rowPermutation(coeffs) = ret )
* @example
* m1 = JSMatrix.mat([[1,2,4,9],[1,8,1,1],[9,4,5,3],[1,2,6,2]]);
* m2 = m1.implicitPartialPivotPermutation();
* // m2.mat = Matrix([[9,4,5,3],[1,8,1,1],[1,2,6,2],[1,2,4,9]])
* // m2.coeffs = [2,1,3,0]
* // m1 = Matrix([[1,2,4,9],[1,8,1,1],[9,4,5,3],[1,2,6,2]])
* m3 = m1.rowPermutation(m2.coeffs); // m3 = Matrix([[9,4,5,3],[1,8,1,1],[1,2,6,2],[1,2,4,9]])
* b = m3.isEqualTo(m2.mat); // b = true
* m4 = m2.mat.rowPermutation(m2.coeffs,true); // m4 = Matrix([[1,2,4,9],[1,8,1,1],[9,4,5,3],[1,2,6,2]])
* b = m4.isEqualTo(m1); // b = true
* m1.implicitPartialPivotPermutation(false);
* // m1 = Matrix([[9,4,5,3],[1,8,1,1],[1,2,6,2],[1,2,4,9]])
*/
JSMatrix.Matrix.prototype.implicitPartialPivotPermutation = function(saveOrig) {
	saveOrig = saveOrig==undefined? true : saveOrig;
	var n = this.nRows;
	var ret = saveOrig? this.copy() : this;
	var vv = [];
	var coeffs = JSMatrix.range(n);
	var big, i, j, k, temp, imax;
	for (i=0; i<n; i++) {
		big = 0.;
		for (j=0; j<n; j++) {
			if ((temp=Math.abs(ret[i][j])) > big) { big=temp; }
		}
		if (big == 0.0) { /*jsmLogger.warn('Matrix.implicitPartialPivotPermutation: singular matrix');*/ return null; }
		vv[i]=1.0/big;
	}
	for (k=0; k<n; k++) {
		big=0.0;
		for (i=k; i<n; i++) {
			temp = vv[i]*Math.abs(ret[i][k]);
			if (temp > big) {
				big = temp;
				imax = i;
			}
		}
		if (k != imax) {
			ret.swapRows(k,imax);
			vv[imax]=vv[k];
			temp = coeffs[k]
			coeffs[k] = imax;
			coeffs[imax] = temp;
		}
	}
	return {mat:ret,coeffs:coeffs};
};

/** Returns vector x as a solution of this*ret = vec, receiver is assumed to be lower triangular matrix
* @param {JSMatrix.Vector} vec right hand side
* @param {boolean=} saveOrig =true] if true, returns vec will be unchanged and new vector will be returned. If false, solution will be saved to vec and vec will be returned
* @returns {JSMatrix.Vector} solution of this*ret = vec
* @example
* a = JSMatrix.mat([[1,0,0,0],[2,3,0,0],[4,5,6,0],[7,8,9,10]]);
* y = JSMatrix.vec([4,17,43,80]);
* x = a.forwardSubstitution(y);
* // x = Vector([4,3,2,1])
* // y = Vector([4,17,43,80])
* check = a.mul(x); // check = Vector([4,17,43,80])
* yy = y.copy();
* // yy = Vector([4,17,43,80])
* a.forwardSubstitution(yy,false);
* // yy = Vector([4,3,2,1])
*/
JSMatrix.Matrix.prototype.forwardSubstitution = function(vec,saveOrig) {
	if (!this.isSquare() || !this.canMultiplyVec(vec)) { /*jsmLogger.warn('Matrix.forwardSubstitution: dimensions mismatched');*/ return null; }
	saveOrig = saveOrig==undefined? true : saveOrig;
	var temp;
	var n = this.nRows;
	var ret = saveOrig? vec.copy() : vec;
	var i, j;
	/* this is lower triangular
	* this*ret = vec
	*
	* |t00  0   0  0 ... 0 |   |ret0|   |vec0|
	* |t10 t11  0  0 ... 0 | * |ret1| = |vec1|
	* |t20 t21 t22 0 ... 0 |   |ret2|   |vec2|
	*
	* ret0 = vec0/t00
	* ret1 = (vec1-t10*ret0)/t11
	* ret2 = (vec2-t20*ret0-t21*ret1)/t22
	* ...
	* retI = (vecI - sum_{j=0,I-1}(tIj*retj))/tII
	*/
	for (i=0; i<n; i++) {
		temp = ret[i];
		for (j=0; j<i; j++) {
			temp -= this[i][j]*ret[j];
		}
		ret[i] = temp/this[i][i];
	};
	return ret;
};

/** Returns vector x as a solution of this*ret = vec, receiver is assumed to be upper triangular matrix
* @param {JSMatrix.Vector} vec right hand side
* @param {boolean=} saveOrig =true] if true, returns vec will be unchanged and new vector will be returned. If false, solution will be saved to vec and vec will be returned
* @returns {JSMatrix.Vector} solution of this*ret = vec
* @example
* a = JSMatrix.mat([[1,2,3,4],[0,5,6,7],[0,0,8,9],[0,0,0,10]]);
* y = JSMatrix.vec([20,34,25,10]);
* x = a.backwardSubstitution(y);
* // x = Vector([4,3,2,1])
* // y = Vector([20,34,25,10])
* check = a.mul(x); // check = Vector([20,34,25,10])
* yy = y.copy();
* // yy = Vector([20,34,25,10])
* a.backwardSubstitution(yy,false);
* // yy = Vector([4,3,2,1])
*/
JSMatrix.Matrix.prototype.backwardSubstitution = function(vec,saveOrig) {
	if (!this.isSquare() || !this.canMultiplyVec(vec)) { /*jsmLogger.warn('Matrix.backwardSubstitution: dimensions mismatched');*/ return null; }
	saveOrig = saveOrig==undefined? true : saveOrig;
	var temp;
	var n = this.nRows;
	var ret = saveOrig? vec.copy() : vec;
	var i, j;
	/* this is upper triangular
	* this*ret = vec
	*
	* | ...   (N, M=N-1, L=N-2)
	* |0 ... tLL tLM tLN|   |retL|   |vecL|
	* |0 ...  0  tMM tMN| * |retM| = |vecM|
	* |0 ...  0   0  tNN|   |retN|   |vecN|
	*
	* retN = vecN/tNN
	* retM = (vecM-tMN*retN)/tMM
	* retL = (vecL-tLN*retN-tLM*retM)/tMM
	* ...
	* retI = (vecI-sum_{j=I+1,N}(tIj*retj))/tII
	*/
	for (i=n-1; i>=0; i--) {
		temp = ret[i];
		for (j=i+1; j<n; j++) {
			temp -= this[i][j]*ret[j];
		}
		ret[i] = temp/this[i][i];
	}
	return ret;
};

/** Returns vector as a solution of system of linear equations using gaussian elimination method (this * ret = rhs --> ret )
* @param {JSMatrix.Vector} rhs vector of right hand sides
* @param {string=} pivoting ="default"] what type of pivoting to use. "default" and "nopivot" stands for no pivoting, "partpivot" for implicit partial pivoting (only row interchanges)
* @param {boolean=} saveOrig =true] if true, returns vec will be unchanged and new vector will be returned. If false, solution will be saved to vec and vec will be returned
* @returns {JSMatrix.Vector} vector of solution
* @example
* a = JSMatrix.mat([[1,2,9],[8,3,2],[3,7,3]]);
* b = JSMatrix.vec([32,20,26]);
* x = a.gaussianElimination(b); // x = Vector([1,2,3])
* x = a.gaussianElimination(b,'partpivot'); // x = Vector([1,2,3])
* a.gaussianElimination(b,'nopivot',false);
* // a = Matrix([[1,2,9],[8,-13,-70],[3,1,-29.384615384615387]])
* // b = Vector([1,2,3])
*/
JSMatrix.Matrix.prototype.gaussianElimination = function(rhs,pivoting,saveOrig) {
	if (!this.canMultiplyVec(rhs) || !this.isSquare()) { /*jsmLogger.warn('Matrix.gaussianElimination: dimensions mismatched');*/ return null; }
	pivoting = pivoting || "default";
	saveOrig = saveOrig==undefined? true : saveOrig;
	var n = this.nRows;
	var i, j, k, akk;
	if (pivoting == "nopivot" || pivoting=="default") {
		var a = saveOrig? this.copy() : this;
		var b = saveOrig? rhs.copy() : rhs;
		for (k=0; k<n-1; k++) {
			akk = a[k][k]; // "pivot"
			for (i=k+1; i<n; i++) {
				for (j=k+1; j<n; j++) {
					a[i][j] -= a[k][j]*a[i][k]/akk;
				}
				b[i] -= b[k]*a[i][k]/akk;
			}
		}
		return a.backwardSubstitution(b,false);
	}
	if (pivoting == "partpivot") {
		var a = this.implicitPartialPivotPermutation(saveOrig);
		var b = saveOrig? rhs.copy() : rhs;
		b.permuteElements(a.coeffs);
		return a.mat.gaussianElimination(b,"nopivot",false);
	}
	return this.gaussianElimination(rhs,"default",saveOrig);
};

/** Returns matrix, whose columns are solutions of system of equations this*ret = rhs
* @param {JSMatrix.Matrix|Array.<JSMatrix.Vector>|JSMatrix.Vector} rhs matrix/vector/array of vectors representing right hand sides
* @param {string=} pivoting ="default"] what type of pivoting to use. "default" and "nopivot" stands for no pivoting, "partpivot" for implicit partial pivoting (only row interchanges)
* @param {boolean=} saveOrig =true] if true, receiver is not changed. If false, solution will be saved to rhs and rhs will be returned
* @returns {JSMatrix.Vector|Array.<JSMatrix.Vector>|JSMatrix.Matrix} matrix, whose columns are solution for particular right hand sides
* @example
* a = JSMatrix.mat([[1,2,9],[8,3,2],[3,7,3]]);
* b1 = JSMatrix.vec([32,20,26]);
* b2 = JSMatrix.vec([16,32,26]);
* b3 = [b1.copy(),b2.copy()];
* b4 = JSMatrix.mat([b1.toArray(),b2.toArray()]).T();
* x1 = a.gaussJordanElimination(b1);
* // x1 = Vector([ 1, 2, 3 ])
* x2 = a.gaussJordanElimination(b2);
* // x2 = Vector([ 3, 2, 1 ])
* x3 = a.gaussJordanElimination(b3);
* // x3 = [Vector([ 1, 2, 3 ]) ,Vector([ 3, 2, 1 ])]
* x4 = a.gaussJordanElimination(b4);
* // x4 = Matrix([[ 1, 3 ], [ 2, 2 ], [ 3, 1 ]])
* x5 = a.gaussJordanElimination(b1,'nopivot');;
* // x5 = Vector([ 1, 2, 3 ])
* x6 = a.gaussJordanElimination(b2,'nopivot');
* // x6 = Vector([ 3, 2, 1 ])
* x7 = a.gaussJordanElimination(b3,'nopivot');
* // x7 = Vector([ 1, 2, 3 ]) ,Vector([ 3, 2, 1 ])
* x8 = a.gaussJordanElimination(b4,'nopivot');
* // x8 = Matrix([[ 1, 3 ], [ 2, 2 ], [ 3, 1 ]])
* a.gaussJordanElimination(b4,'nopivot',false);
* // b4 = Matrix([[ 1, 3 ], [ 2, 2 ], [ 3, 1 ]])
*/
JSMatrix.Matrix.prototype.gaussJordanElimination = function(rhs,pivoting,saveOrig) {
	if (!this.isSquare) { /*jsmLogger.warn('Matrix.gaussJordanElimination: dimensions mismatched');*/ return null; }
	var rhsInstance = (rhs instanceof JSMatrix.Vector)? 1 : ((rhs instanceof Array)? 2 : ((rhs instanceof JSMatrix.Matrix)? 3 : 4));
	if (!(rhsInstance==1 || rhsInstance==2 || rhsInstance==3)) { /*jsmLogger.warn('Matrix.gaussJordanElimination: wrong argument');*/ return null; }
	pivoting = pivoting || "default";
	saveOrig = saveOrig==undefined? true : saveOrig;
	var i, j, k, akk;
	var n = this.nRows;
	var nrhs = rhsInstance==2? rhs.length : (rhsInstance==3? rhs.nCols : 1);
	if (pivoting=="default" || pivoting=="nopivot") {
		var a = saveOrig? this.copy() : this;
		var b;
		if (rhsInstance==1 || rhsInstance==3) { b = saveOrig? rhs.copy() : rhs; }
		else {
			b = [];
			for (i=0; i<nrhs; i++) { b[i] = saveOrig? rhs[i].copy() : rhs[i]; }
		}
		for (k=0; k<n; k++) {
			akk = a[k][k];
			for (i=k+1; i<n; i++) {
				for (j=k+1; j<n; j++) {
					a[i][j] -= a[k][j]*a[i][k]/akk;
				}
				switch (rhsInstance) {
					case 1: b[i] -= b[k]*a[i][k]/akk; break;
					case 2: for (j=0; j<nrhs; j++) { b[j][i] -= b[j][k]*a[i][k]/akk; } break;
					case 3: for (j=0; j<nrhs; j++) { b[i][j] -= b[k][j]*a[i][k]/akk; } break;
				}
			}
		}
		for (k=n-1; k>=0; k--) {
			akk = a[k][k];
			for (i=k-1; i>=0; i--) {
				for (j=n-1; j>k; j--) {
					a[i][j] -= a[k][j]*a[i][k]/akk;
				}
				switch (rhsInstance) {
					case 1: b[i] -= b[k]*a[i][k]/akk; break;
					case 2: for (j=0; j<nrhs; j++) { b[j][i] -= b[j][k]*a[i][k]/akk; } break;
					case 3: for (j=0; j<nrhs; j++) { b[i][j] -= b[k][j]*a[i][k]/akk; } break;
				}
			}
		}
		for (k=0; k<n; k++) {
			akk = a[k][k];
			switch (rhsInstance) {
				case 1: b[k] /= akk; break;
				case 2: for (j=0; j<nrhs; j++) { b[j][k] /= akk; } break;
				case 3: for (j=0; j<nrhs; j++) { b[k][j] /= akk; } break;
			}
		}
		return b;
	}
	if (pivoting == "partpivot") {
		var a = this.implicitPartialPivotPermutation(saveOrig);
		var temp;
		switch (rhsInstance) {
			case 1: temp = saveOrig? rhs.copy() : rhs; temp.permuteElements(a.coeffs); break;
			case 2:
				temp = [];
				for (i=0; i<nrhs; i++) {
					temp[i] = saveOrig? rhs[i].copy() : rhs[i];
					temp[i].permuteElements(a.coeffs);
				}
				break;
			case 3:
				temp = saveOrig? rhs.copy() : rhs;
				temp.permuteRows(a.coeffs);
				break;
		}
		return a.mat.gaussJordanElimination(temp,"nopivot",false);
	}
	return this.gaussJordanElimination(rhs,"default",saveOrig);
};

/** Returns lower triangular matrix of Cholesky decomposition of receiver. Receiver must be square, symmetric and positive definite. ( ret * ret^T = this )
* @param {boolean=} saveOrig =true] if true, receiver is not changed. If false, solution will be saved to receiver
* @param {boolean=} checkSymm =false] if true, test of receiver symmetricity is performed, otherwise not
* @returns {JSMatrix.Matrix} lower triangular Cholesky matrix
* @example
* a = JSMatrix.mat([[1,2,4],[2,13,23],[4,23,77]]);
* l = a.choleskyDecomposition(); // l = Matrix([[1,0,0],[2,3,0],[4,5,6]])
* check = l.x(l.T()); // check = Matrix([[1,2,4],[2,13,23],[4,23,77]]);
*/
JSMatrix.Matrix.prototype.choleskyDecomposition = function(saveOrig,checkSymm) {
	checkSymm = checkSymm==undefined? true : checkSymm;
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.choleskyDecomposition: matrix is not square');*/ return null; }
	if (checkSymm) {
		if (!this.isSymmetric()) { /*jsmLogger.warn('Matrix.choleskyDecomposition: matrix is not symmetric');*/ return null; }
	}
	saveOrig = saveOrig==undefined? true : saveOrig;
	var n = this.nRows;
	var a = saveOrig? this.copy() : this;
	var aii, temp, i, j, k;
	/*
	* |t00 t01 t02 ..| |l00  0   0  ..| |l00 l10 l20 ..|
	* |t01 t11 t12 ..|=|l10 l11  0  ..|*| 0  l11 l21 ..|
	* |t02 t12 t22 ..| |l20 l21 l22 ..| | 0   0  l22 ..|
	* |  ..........  | |  ..........  | |  ..........  |
	*
	* t00 = l00*l00  ->  l00 = sqrt(t00)
	* t01 = l00*l10  ->  l10 = t01/l00
	* t02 = l00*l20  ->  l20 = t02/l00
	* t11 = l10*l10+l11*l11 -> l11 = sqrt(t11-l00*l00)
	* t12 = l10*l20+l11*l21 -> l21 = (t12-l10*l20)/l11
	* t22 = l20*l20+l21*l21+l22*l22 -> l22=sqrt(t22-l20*l20-l21*l21)
	* t32 = l30*l20+l31*l21+l32*l22 -> l32 = (t32-l30*l20-l31*l21)/l22
	* ...
	* lii = sqrt(tii - sum_{k=0,i-1}(lik*lik))
	* lij = (tij-sum_{k=0,j-1}(lik*ljk))/lii
	*/
	for (i=0; i<n; i++) {
		aii = a[i][i];
		for (k=0; k<i; k++) { aii -= a[i][k]*a[i][k]; }
		if (aii<=0.) { /*jsmLogger.warn('Matrix.choleskyDecomposition: matrix is not positive definite');*/ return null; }
		aii = Math.sqrt(aii);
		a[i][i] = aii;
		for (j=i+1; j<n; j++) {
			temp = a[j][i];
			for (k=0; k<i; k++) { temp -= a[i][k]*a[j][k]; }
			a[j][i] = temp/aii;
		}
		for (k=i+1; k<n; k++) { a[i][k] = 0.; }
	}
	return a;
};

/** Alias for Cholesky decomposition, see {@JSMatrix.Matrix#choleskyDecomposition}
* @function
*/
JSMatrix.Matrix.prototype.chol = JSMatrix.Matrix.prototype.choleskyDecomposition;

/** Alias for Cholesky decomposition, see {@JSMatrix.Matrix#choleskyDecomposition}
* @function
*/
JSMatrix.Matrix.prototype.lltDecomposition = JSMatrix.Matrix.prototype.choleskyDecomposition;

/** Returns lower (L) and upper (U) triangular matrix of the receiver such that L*U=this. Diagonal elements of L are all 1. Receiver must be square.
* @param {string=} pivoting ="default"] what type of pivoting to use. "default" and "nopivot" stands for no pivoting, "partpivot" for implicit partial pivoting (only row interchanges)
* @returns {{l:JSMatrix.Matrix,u:JSMatrix.Matrix}|null} {l,u}, lower and upper triangular matrices
* @example
* a = JSMatrix.mat([[6,5,4],[12,13,10],[18,27,21]]);
* lu = a.luDecomposition();
* l = lu.l; // l = Matrix([[1,0,0],[2,1,0],[3,4,1]])
* u = lu.u; // u = Matrix([[6,5,4],[0,3,2],[0,0,1]])
* check = l.x(u); // check = Matrix([[6,5,4],[12,13,10],[18,27,21]])
*/
JSMatrix.Matrix.prototype.luDecomposition = function(pivoting) {
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.luDecomposition: matrix is not square');*/ return null; }
	pivoting = pivoting || "default";
	var n = this.nRows;
	var l = JSMatrix.Matrix.Zeros(n,n);
	var u = JSMatrix.Matrix.Zeros(n,n);
	var i, j, k, imax;
	var temp;
	/*
	* |t00 t01 t02 ..| |l00  0   0  ..| |u00 u01 u02 ..|
	* |t10 t11 t12 ..|=|l10 l11  0  ..|*| 0  u11 u12 ..|
	* |t20 t21 t22 ..| |l20 l21 l22 ..| | 0   0  u22 ..|
	*
	* assumtion: lkk = 1.
	* t00 = l00*u00 -> l00=1, u00=t00
	* t01 = l00*u01 -> u01=t01/l00
	* t02 = l00*u02 -> u02=t02/l00
	* ...
	* t10 = l10*u00 -> l10=t10/u00
	* t11 = l10*u01+l11*u11 -> l11=1, u11=(t11-l10*u01)
	* t12 = l10*u02+l11*u12 -> u12=t12-l10*u02
	* t13 = l10*u03+l11*u13 -> u13=t13-l10*u03
	* ...
	* t20 = l20*u00 -> l20=t10/u00
	* t21 = l20*u01+l21*u11 -> l21=(t21-l20*u01)/u11
	* t22 = l20*u02+l21*u12+l22*u22 -> l22=1, u22=(t22-l20*u02-l21*u12)
	* ...
	* lkk = 1.
	* ukj = tkj - sum_{i=0,k-1}(lki*uij)
	* ljk = tjk - sum_{i=0,k-1}(lji*uik)
	*/
	if (pivoting=="default" || pivoting=="nopivot") {
		for (k=0; k<n; k++) {
			l[k][k] = 1.;
			for (j=k; j<n; j++) {
				temp = this[k][j];
				for (i=0; i<k; i++) {
					temp -= l[k][i]*u[i][j];
				}
				u[k][j] = temp;
			}
			for (j=k+1; j<n; j++) {
				temp = this[j][k]
				for (i=0; i<k; i++) {
					temp -= l[j][i]*u[i][k]
				}
				l[j][k] = temp/u[k][k]
			}
		}
		return {l:l,u:u};
	}
	if (pivoting == "partpivot") {
		var a = this.implicitPartialPivotPermutation()
		var ret = a.mat.luDecomposition("nopivot");
		ret.coeffs = a.coeffs;
		return ret;
	}
	return this.luDecomposition("default");
};

/** Returns lower triangular (L), diagonal (D) and upper triangular (U) matrix of the receiver such that L*D*U=this. Diagonal elements of both L and U are all 1. Receiver must be square.
* @param {string=} pivoting ="default"] what type of pivoting to use. "default" and "nopivot" stands for no pivoting, "partpivot" for implicit partial pivoting (only row interchanges)
* @returns {{l:JSMatrix.Matrix,d:JSMatrix.Vector,u:JSMatrix.Matrix}|null} {l,d,u}, lower, diagonal and upper triangular matrices
* @example
* a = JSMatrix.mat([[2,10,8],[4,23,22],[6,42,52]]);
* ldu = a.lduDecomposition();
* l = ldu.l; // l = Matrix([[1,0,0],[2,1,0],[3,4,1]])
* d = ldu.d; // d = Vector([2,3,4])
* u = ldu.u; // u = Matrix([[1,5,4],[0,1,2],[0,0,1]])
* check = l.x(d.diag()).x(u) // check = Matrix([[2,10,8],[4,23,22],[6,42,52]]);
*/
JSMatrix.Matrix.prototype.lduDecomposition = function(pivoting) {
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.lduDecomposition: matrix is not square');*/ return null; }
	pivoting = pivoting || "default";
	if (pivoting=="default" || pivoting=="nopivot") {
		var n = this.nRows;
		var l = JSMatrix.Matrix.Zeros(n,n);
		var d = JSMatrix.Vector.Zeros(n);
		var u = JSMatrix.Matrix.Zeros(n,n);
		var temp;
		for (var i=0; i<n; i++) {
			l[i][i] = 1.;
			u[i][i] = 1.;
			for (var j=0; j<i; j++) {
				temp = this[i][j];
				for (var k=0; k<j; k++) {
					temp -= l[i][k]*d[k]*u[k][j];
				}
				l[i][j] = temp/d[j];
			}
			temp = this[i][i];
			for (var k=0; k<i; k++) {
				temp -= l[i][k]*d[k]*u[k][i];
			}
			d[i] = temp;
			for (var j=i+1; j<n; j++) {
				temp = this[i][j];
				for (var k=0; k<i; k++) {
					temp -= l[i][k]*d[k]*u[k][j];
				}
				u[i][j] = temp/d[i];
			}
		}
		return {l:l,d:d,u:u};
	}
	if (pivoting = "partpivot") {
		var a = this.implicitPartialPivotPermutation();
		var ret = a.mat.lduDecomposition("nopivot");
		ret.coeffs = a.coeffs;
		return ret;
	}
	return this.lduDecomposition("default");
};

/* Returns orthogonal matrix Q and upper triangular matrix R of QR decomposition of receiver such that Q*R=this. Implemented for square matrices only
* @param {string=} method ="default"] method for solution. "default" and "grammschmidt" stand for Gramm-Schmidt process, "householder" for using Householder transformations
* @param {boolean=} saveOrig =true] if true, receiver is not changed. If false, solution (R matrix) will be saved to this and this will be returned
* @returns {{q:JSMatrix.Matrix,r:JSMatrix.Matrix}|null} {q,r} of QR decomposition of receiver
* @example
* a = JSMatrix.mat([[0.8,1.6,6],[0,4,5],[-0.6,-1.2,3]])
* qr = a.qrDecomposition();
* q = qr.q; // q = Matrix([[0.8,0,0.6],[0,1,0],[-0.6,0,0.8]])
* r = qr.r; // r = Matrix([[1,2,3],[0,4,5],[0,0,6]])
* check = q.x(r); // check = Matrix([[0.8,1.6,6],[0,4,5],[-0.6,-1.2,3]])
* o1 = q.isOrthogonal(); // o1 = true
* o2 = r.isUpperTriangular(); // o2 = true
* ///
* qr = a.qrDecomposition('householder');
* q = qr.q; // q = Matrix([[ 0.11624763874381927, 0.24071784839360483, 0.9636085325230572 ], [ 0.9299811099505543, -0.3670351351744072, -0.02050230920261827 ], [ 0.34874291623145787, 0.8985210776672174, -0.26653001963403705 ]])
* r = qr.r; // r = Matrix([[ 8.602325267042627, 5.463639020959507, 5.1148961047280475 ], [ 2.5252117732828282e-17, 5.669977834934511, 1.5968411725120326 ], [ -2.206040279013976e-16, 4.440892098500626e-16, 8.323937536263006 ]])
* check = q.x(r); // check = Matrix([[1,2,9],[8,3,4],[3,7,1]])
* o1 = q.isOrthogonal(); // o1 = true
* o2 = r.isUpperTriangular(); // o2 = true
*/
JSMatrix.Matrix.prototype.qrDecomposition = function(method,saveOrig) {
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.qrDecomposition: matrix is not square');*/ return null; }
	saveOrig = saveOrig==undefined? true : saveOrig;
	method = method || "default";
	var n = this.nRows;
	if (method=="default" || method=="grammschmidt" || method=="gramm-schmidt" || method=="gs") {
		//TODO
		var proj = function(ee,aa) { return ee.mulf(ee.dot(aa)/(ee.dot(ee))); }
		var r = saveOrig? this.copy() : this;
		var q = new JSMatrix.Matrix(n,n);
		var a, e=[];
		var j,k;
		for (k=0; k<n; k++) {
			a = r.getCol(k);
			e[k] = a.copy();
			for (j=0; j<k; j++) { e[k].isub(proj(e[j],a)); }
			e[k].normalize();
			for (j=0; j<=k; j++) { r[j][k] = e[j].dot(a); }
			for (j=k+1; j<n; j++) { r[j][k] = 0.; }
			q.setCol(k,e[k]);
		}
		return {q:q,r:r};
	}
	if (method=="householder") {
		//TODO
		var r = saveOrig? this.copy() : this;
		var e = JSMatrix.Vector.Zeros(n);
		e[0] = 1.;
		var a, x, rkn, qtemp, q=JSMatrix.Matrix.Identity(n), piv;
		// TODO operation reduction
		for (var k=0; k<n-1; k++) {
			e.resize(n-k);
			rkn = JSMatrix.range(k,n);
			a = r.getSubMatrix(rkn,rkn);
			x = a.getCol(0);
			qtemp = JSMatrix.Matrix.Identity(n);
			piv = x[k]>0.? -1. : 1.;
			qtemp.setSubMatrix(rkn,rkn,JSMatrix.Matrix.Householder(x.add(e.mulf(piv*x.norm()))));
			r = qtemp.mulm(r);
			q = q.mulm(qtemp.transposed());
		}
		return {q:q,r:r};
	}
	return this.qrDecomposition("default",saveOrig);
};

/** Returns vector/array of vectors/matrix as a solution of system of linear equations (this*ret=rhs --> ret). Receiver must be square. TODO: solution with LU decomposition with permutation
* @param {JSMatrix.Vector|Array.<JSMatrix.Vector>|JSMatrix.Matrix} rhs vector of right hand sides. Array of vectors or Matrix as rhs is supported only with "gaussjordan" method
* @param {string=} method ="default"] method of the solution. Implemented methods are "default" and "gauss" for gaussian elimination, "gaussjordan" for Gauss-Jordan elimination (supporting more than one right hand sides in form of array of vectors or matrix with columns as rhs), "lu" for using LU decomposition, "ldu" for using LDU decomposition, "cholesky" for using Cholesky decomposition, "qr" for using QR decomposition
* @param {boolean=} saveOrig =true] if true, receiver and rhs is not changed. If false, solution will be saved to rhs and rhs will be returned and receiver will be changed
* @param {Array.<JSMatrix.Matrix>=} precompDecomps =undefined] array of Matrices objects, used for solution using matrix decomposition as precomputed values (the decomposition is performed within linSolve function if args parameter is not specified)
* @returns {JSMatrix.Vector|Array.<JSMatrix.Vector>|JSMatrix.Matrix} vector/array of vectors/matrix of solution
* @example
* a = JSMatrix.mat([[1,2,9],[8,3,2],[3,7,3]]);
* b = JSMatrix.vec([32,20,26]);
* x1 = a.linSolve(b); // x1 = Vector([1,2,3])
* x2 = a.linSolve(b,'gauss'); // x2 = Vector([1,2,3])
* x3 = a.linSolve(b,'gaussjordan'); // x3 = Vector([1,2,3])
* x4 = a.linSolve(b,'lu'); // x4 = Vector([1,2,3])
* x5 = a.linSolve(b,'ldu'); // x5 = Vector([1,2,3])
* x7 = a.linSolve(b,'qr'); // x7 = Vector([1,2,3])
* a2 = JSMatrix.mat([[1,2,4],[2,13,23],[4,23,77]]);
* b2 = JSMatrix.vec([17,97,281]);
* x6 = a2.linSolve(b2,'cholesky'); // x6 = Vector([1,2,3])
*/
JSMatrix.Matrix.prototype.linSolve = function(rhs,method,saveOrig,precompDecomps) {
	saveOrig = saveOrig==undefined? true : saveOrig;
	method = method || "default";
	method = method.toLowerCase();
	if (method=="default" || method=="gauss" || method=="gaussian") { // gaussian elimination
		return this.gaussianElimination(rhs,"partpivot",saveOrig);
	}
	if (method=="gaussjordan" || method=="gauss-jordan" || method=="gj") { // gauss-jordan elimination
		return this.gaussJordanElimination(rhs,"partpivot",saveOrig);
	}
	if (method=="lu") { // use LU factorization
		var lu = precompDecomps || this.luDecomposition();
		var l = lu.l;
		var u = lu.u;
		var y;
		if (rhs instanceof JSMatrix.Vector) {
			y = l.forwardSubstitution(rhs,saveOrig);
			return u.backwardSubstitution(y,false);
		}
		if (rhs instanceof JSMatrix.Matrix) {
			var ret = saveOrig? JSMatrix.Matrix.Zeros(this.nRows,rhs.nCols) : rhs;
			for (var i=0; i<ret.nCols; i++) {
				y = l.forwardSubstitution(rhs.getCol(i),false);
				ret.setCol(i,u.backwardSubstitution(y,false));
			}
			return ret;
		}
		if (rhs instanceof Array) {
			return null;
			// TODO
		}
	}
	if (method=="ldu") { // use LDU decomposition
		var ldu = precompDecomps || this.lduDecomposition("nopivot");
		var l = ldu.l;
		var d = (ldu.d instanceof JSMatrix.Matrix)? ldu.d.diag() : ldu.d;
		var u = ldu.u;
		var y = l.forwardSubstitution(rhs,saveOrig);
		for (var r=0; r<this.nRows; r++) { y[r] /= d[r]; }
		return u.backwardSubstitution(y,false);
	}
	if (method=="cholesky" || method=="chol") { // use Cholesky's decomposition
		var l = precompDecomps || this.choleskyDecomposition();
		if (l==null) { /*jsmLogger.warn('Matrix.linSolve: matrix is not positive definite');*/ return null; }
		var y = l.forwardSubstitution(rhs,saveOrig);
		return l.transposed().backwardSubstitution(y,false);
	}
	if (method=="qr") { // use QR decomposition
		var qr = precompDecomps || this.qrDecomposition("grammschmidt",saveOrig);
		return qr.r.backwardSubstitution(qr.q.transposed().mulv(rhs),false);
	}
	return this.linSolve(rhs,"default",saveOrig,precompDecomps);
};

/** returns Schur form of receiver (upper triangular matrix with eigenvalues of the receiver on the diagonal). QR algorithm is used for calculation
* @param {boolean=} saveOrig =true] if true, returns new matrix, if false the Schur form is formed on original matrix
* @param {number=} maxiter =1000] maximum number of iterations
* @returns {JSMatrix.Matrix} Schur form of receiver (upper triangular matrix with eigenvalues on diagonal)
* @example
* m1 = JSMatrix.mat([[4,2,3,1],[3,6,5,7],[5,3,8,1],[2,3,4,5]]);
* m2 = m1.schurForm(); // m2 = Matrix([[ 15.6072, -3.6605, -0.7006, 2.3149 ], [ 0, 4.1916, 3.3807, 2.0775 ], [ 0, 0, 2.2012, 0.0769 ], [ 0, 0, 0, 1 ]])
* b = m2.isUpperTriangular(); // b = true
*/
JSMatrix.Matrix.prototype.schurForm = function(saveOrig,maxiter) {
	//return this.schurDecomposition(saveOrig,maxiter)[1];
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.schurForm: matrix is not square');*/ return null; }
	var n = this.nRows;
	saveOrig = saveOrig==undefined? true : saveOrig;
	maxiter = maxiter || 1000;
	var qr, s = saveOrig? this.copy() : this;
	var rhos, rhosNew=JSMatrix.Vector.Zeros(n), rhosRatio=JSMatrix.Vector.Zeros(n), maxRhosRatio, i, q, r, iter=0, qq=JSMatrix.Matrix.Identity(n);
	for (i=0; i<n; i++) { rhosNew[i] = 0.; }
	do {
		qr = s.qrDecomposition()
		q = qr.q;
		r = qr.r;
		s = r.mulm(q);
		rhos = rhosNew.copy();
		rhosNew = s.diagonalToVector();
		for (i=0; i<n; i++) { rhosRatio[i] = Math.abs((rhosNew[i] - rhos[i])/rhosNew[i]); }
		maxRhosRatio = 0.;
		for (i=0; i<n; i++) {
			if (rhosRatio[i] > maxRhosRatio) { maxRhosRatio = rhosRatio[i]; }
		}
		iter++;
		if (iter > maxiter) {
			if (maxRhosRatio > JSMatrix.TOL) { /*jsmLogger.warn('Matrix.schurForm: maximum number of iterations reached');*/ return null; }
			break;
		}
	} while (maxRhosRatio > JSMatrix.TOL2)
	return s;
};

/** returns orthogonal matrix Q and upper triangular matrix S as Schur decomposition of receiver such that this=Q*S*Q^-1. Columns of Q are eigenvectors of receiver and diagonal elements of S are eigenvalues of receiver 
* @param {boolean=} saveOrig =true] if true, returns new matrices, if false, S is formed into original matrix
* @param {number=} maxiter =1000] maximum number of iterations
* @returns {{q:JSMatrix.Matrix,s:JSMatrix.Matrix}|null} {q,s} - Scur decomposition of receiver
* @example
* m1 = JSMatrix.mat([[4,2,3,1],[3,6,5,7],[5,3,8,1],[2,3,4,5]]);
* qs = m1.schurDecomposition();
* q = qs.q; s = qs.s;
* // s = Matrix([[ 15.6072, -3.6605, -0.7006, 2.3149 ], [ 0, 4.1916, 3.3807, 2.0775 ], [ 0, 0, 2.2012, 0.0769 ], [ 0, 0, 0, 1 ]])
* qr = s.qrDecomposition();
* q1 = qr.q; r = qr.r;
* // q1 = Matrix([[ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ]])
* // r = Matrix([[ 15.6072, -3.6605, -0.7006, 2.3149 ], [ 0, 4.1916, 3.3807, 2.0775 ], [ 0, 0, 2.2012, 0.0769 ], [ 0, 0, 0, 1 ]])
* b1 = q.isOrthogonal() // b1 = true
* b2 = s.isUpperTriangular() // b2 = true
* c = q.mulm(s).mulm(q.inv()); // c = Matrix([[4,2,3,1],[3,6,5,7],[5,3,8,1],[2,3,4,5]])
*/
JSMatrix.Matrix.prototype.schurDecomposition = function(saveOrig,maxiter) {
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.schurDecomposition: matrix is not square');*/ return null; }
	var n = this.nRows;
	saveOrig = saveOrig==undefined? true : saveOrig;
	maxiter = maxiter || 1000;
	var qr, s = saveOrig? this.copy() : this;
	var rhos, rhosNew=JSMatrix.Vector.Zeros(n), rhosRatio=JSMatrix.Vector.Zeros(n), maxRhosRatio, i, q, r, iter=0, qq=JSMatrix.Matrix.Identity(n);
	for (i=0; i<n; i++) { rhosNew[i] = 0.; }
	do {
		// A(k)=Q(k)*R(k), A(k+1)=R(k)*Q(k)=Q(k)^T*A(k)*Q(k)
		// A(k+1) converges to Schur form
		// Furthermore Q=sum_k(Q(k)), S=A(N)=Q^T*this*Q -> this=Q*S*Q^T
		// (Q is othogonal)
		qr = s.qrDecomposition()
		q = qr.q;
		qq = qq.mulm(q);
		r = qr.r;
		s = r.mulm(q);
		rhos = rhosNew.copy();
		rhosNew = s.diagonalToVector();
		for (i=0; i<n; i++) { rhosRatio[i] = Math.abs((rhosNew[i] - rhos[i])/rhosNew[i]); }
		maxRhosRatio = 0.;
		for (i=0; i<n; i++) {
			if (rhosRatio[i] > maxRhosRatio) { maxRhosRatio = rhosRatio[i]; }
		}
		iter++;
		if (iter > maxiter) {
			if (maxRhosRatio > JSMatrix.TOL) { /*jsmLogger.warn('Matrix.schurDecomposition: maximum number of iterations reached');*/ return null; }
			break;
		}
	} while (maxRhosRatio > JSMatrix.TOL2)
	return {q:qq,s:s};
};

/** Returns orthogonal matrix U, diagonal matrix S and orthogonal matrix V of singular value decomposition (SVD) of receiver such that this = U*S*V^T. U is computed as orthonormalized iegenvectors of this^T*this, V as orhonormalized eigenvectors of this*this^T and diagonal components of S as square roots of eigenvalues of this^T*this (TODO faster algorithm). this*V = U*S, A^T*U = V*S. Current implementation allows only square matrices (TODO gegeral-shaped matrices)
* @returns {{u:JSMatrix.Matrix,s:JSMatrix.Vector,v:JSMatrix.Matrix}|null} {u,s,v} orthogonal matrix U, diagonal of matrix S and orthogonal matrix V of svd decomposition of receiver
* @example
* m1 = JSMatrix.mat([[1.6,0.6,0.],[-1.2,0.8,0.],[0.,0.,3.]])
* svd = m1.singularValueDecomposition();
* u = svd.u; // u = Matrix([[0.6,0.8,0],[0.8,-0.6,0],[0,0,1]])
* s = svd.s; // s = Vector([1,2,3])
* v = svd.v; // v = Matrix([[0,1,0],[1,0,0],[0,0,1]])
* b1 = u.isOrthogonal() // b1 = true
* b2 = v.isOrthogonal() // b2 = true
* c = u.mulm(s.diag()).mulm(v.transposed());
* // c = JSMatrix.mat([[1.6,0.6,0.],[-1.2,0.8,0.],[0.,0.,3.]])
*/
JSMatrix.Matrix.prototype.singularValueDecomposition = function() {
	if (!this.isSquare()) { return null; }
	var aat = this.mulm(this.transposed());
	var naat = aat.nRows;
	var ata = this.transposed().mulm(this);
	var nata = ata.nRows;
	var aatEig = aat.eigSymmetric({nEigModes:'all'});
	var ataEig = ata.eigSymmetric({nEigModes:'all'});
	var U = JSMatrix.Matrix.Zeros(aat.nRows,aat.nCols);
	var V = JSMatrix.Matrix.Zeros(ata.nRows,ata.nCols);
	for (var n=0; n<naat; n++) { U.setCol(n,aatEig.eigVecs[n]); }
	for (var n=0; n<nata; n++) { V.setCol(n,ataEig.eigVecs[n]); }
	var s = [];
	for (var n=0; n<nata; n++) { s[n] = Math.sqrt(ataEig.eigVals[n]); }
	return {u:U,s:JSMatrix.Vector.create(s),v:V};
};

/** Alias for singular value decomposition, see {@link JSMatrix.Matrix.singularValueDecomposition}
* @function
*/
JSMatrix.Matrix.prototype.svd = JSMatrix.Matrix.prototype.singularValueDecomposition;

/* * Returns pseudoinverse of receiver 
* @param {string=} method ="defualt"] method of evaluation. Options are 'defualt' for ( ret = (this^T*this)^-1 * this^T ), 'svd' for SVD decomposition ( this=U*S*V^T -> this^+ = V*S^+*U^T, A^+ means pseudoinverse)
* @param {{u:JSMatrix.Matrix,s:JSMatrix.Matrix,v:JSMatrix.Matrix}=} precompSvd =undefined] precomputed SVD decomposition of receiver
* @returns {JSMatrix.Matrix} pseudoinverse of receiver
* @ example
* m1 = JSMatrix.mat([[2,-4,5],[6,0,3],[2,-4,5],[6,0,3]]);
* b = JSMatrix.vec([1,3,-1,3]);
* m2 = m1.pseudoinverse();
* // m2 =
* x2 = m2.mulv(b);
* // x2 = 
* c2 = m1.mulv(x2)
* // c2 =
* /
pseudoinverse: function(method,precompSvd) {
	var method = method || 'defualt'
	method = method.toLowerCase();
	//if (method=='default') {
		//var t = this.inversed();
		//return (t.mulm(this)).inversed().mulm(t);
	//}
	if (method=='svd' || method=='default') {
		var svd = precompSvd? precompSvd : this.singularValueDecomposition();
		var u = svd[0];
		var s = svd[1];
		var v = svd[2];
		var n = s.nElements;
		var ss = new JSMatrix.Matrix(n,n);
		for (var i=0; i<n; i++) {
			if (Math.abs(s[i]) > JSMatrix.TOL) { ss[i][i] = 1./s[i]; }
		}
		return v.mulm(ss).mulm(u.transposed());
	}
	return this.pseudoinverse('default',precompSvd);
};
*/

/** Returns orthogonal matrix R and symmetric matrix U of polar decomposition of the receiver such that this=R*U. The decomposition is computed from SVD decomposition (this = U*S*V^T -> R = U*V^T, U = V*S*V^T
* @param {{u:JSMatrix.Matrix,s:JSMatrix.Matrix,v:JSMatrix.Matrix}=} precompSvd =undefined] precomputed SVD decomposition
* @returns {{r:JSMatrix.Matrix,u:JSMatrix.Matrix}|null} {r,u} orthogonal matrix R and symmetric matrix U of polar decomposition of receiver
* @example
* m = JSMatrix.mat([[1.6,.6,0],[-1.2,0.8,0],[0,0,3]])
* ru = m.polarDecomposition();
* r = ru.r; // r = Matrix([[0.8,0.6,0],[-0.6,0.8,0],[0,0,1]])
* u = ru.u; // u = Matrix([[2,0,0],[0,1,0],[0,0,3]])
* b1 = r.isOrthogonal(); // b1 = true
* b2 = u.isSymmetric(); // b2 = true
* c = r.mulm(u); // c = Matrix([[1.6,.6,0],[-1.2,0.8,0],[0,0,3]])
*/
JSMatrix.Matrix.prototype.polarDecomposition = function(precompSvd) {
	if (!this.isSquare) { /*jsmLogger.warn('Matrix.polarDecomposition: matrix is not square');*/ return null; }
	var svd = precompSvd? precompSvd : this.singularValueDecomposition();
	var u = svd.u;
	var s = svd.s;
	var v = svd.v;
	return {r:u.mulm(v.transposed()), u:v.mulm(s.diag()).mulm(v.transposed())};
};

/** Returns eigenvalues (L) and eigenvectors (P) of symmetric receiver ( (this-lambda*mat)*phi = 0 ). Eigenvectors are normalized with respect to mat
* @param {{
				mat: JSMatrix.Matrix,
				nEigModes: (number|string),
				method: string,
				highest: boolean,
				maxiter: number}=} params 
* @param {JSMatrix.Matrix=} params.mat =JSMatrix.Matrix.Identity] mat in (this-lambda*mat)*phi = 0 ). If mat==Matrix.Identity (default), then the "classical" eigenvalue problem is solved ( (this-lambda*I)*phi = 0 -> this*phi=lambda*phi)
* @param {number|string=} params.nEigModes ='all'] number of first eigenvalues and eigenvectors to be returned. if "all" is, then all eigenvalues are returned
* @param {string=} params.method ="default"] method of the solution. "default" and "invit" stands for Stodola's inverse iterations methods with Gramm-Schmidt orthogonalization, "subspace" for subspace iterations
@param {boolean=} params.highest =true] find params.nEigModes highest eigenvalues if true, params.nEigModes lowest eigenvalues otherwise
* @param {number=} params.maxiter =1000] maximum number of iterations
* @returns {{eigVals:Array.<number>,eigVecs:Array.<JSMatrix.Vector>}|null} [lambda1 ,lambda2,...,lambdaN],[phi1,phi2,...,phiN]], lambdai is i-th eigenvalue, phii is i-th eigenvector. N = nEigModes
* //
* @example
* k = JSMatrix.mat([[9,3,0,0],[3,8,2,0],[0,2,7,2],[0,0,2,8]]);
* m = JSMatrix.Matrix.Diagonal([3,2,4,3]);
* e = k.eigSymmetric({mat:m,nEigModes:2,highest:false});
* ls = e.eigVals;
* ps = e.eigVecs;
* l1 = ls[0]; // l1 = 1.2504038497315775
* l2 = ls[1]; // l2 = 2.279120228657416
* p1 = ps[0]; // p1 = Vector([ 0.1288306767139194, -0.22540165581364102, 0.4265175098731745, -0.20077135620035116 ])
* p2 = ps[1]; // p2 = Vector([ 0.4514687335612619, -0.32545467450354726, -0.11713473474653795, 0.2014979513549984 ])
* ///
* e = k.eigSymmetric({mat:m,nEigModes:2,highest:true,method:'subspace'});
* ls = e.eigVals;
* ps = e.eigVecs;
* l3 = ls[1]; // l3 = 2.948867853467429
* l4 = ls[0]; // l4 = 4.938274734810245
* p3 = ps[1]; // p3 = Vector([ 0.14681878659487543, -0.007507144503266177, -0.21233717286242818, -0.5016212772448967 ])
* p4 = ps[0]; // p4 = Vector([ 0.3022519091588433, 0.585847242190032, 0.09630780190740544, 0.028264211960396433 ])
* ///
* c1 = k.x(p1).sub(m.x(p1).x(l1)); // c1 = Vector([0,0,0,0])
* c2 = k.x(p2).sub(m.x(p2).x(l2)); // c2 = Vector([0,0,0,0])
* c3 = k.x(p3).sub(m.x(p3) .x(l3)); // c3 = Vector([0,0,0,0])
* c4 = k.x(p4).sub(m.x(p4).x(l4)); // c4 = Vector([0,0,0,0])
*/
JSMatrix.Matrix.prototype.eigSymmetric = function(params) {
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.eigSymmetric: matrix is not square');*/ return null; }
	params = params || {};
	var mat = params.mat || JSMatrix.Matrix.Identity(this.nRows);
	if (!mat.isSquare() || !this.isSameSizeAs(mat)) { /*jsmLogger.warn('Matrix.eigSymmetric: dimensions mismatched');*/ return null; }
	var nEigModes = params.nEigModes || this.nRows; // nEigModes first eigenvalues and eigenvectors will be returned
	if (nEigModes=="all") { nEigModes = this.nRows; }
	else if (nEigModes > this.nRows) { nEigModes = this.nRows; }
	var method = params.method || "default";
	var highest = params.highest || false;
	var maxiter = params.maxiter || 1000;
	if (method=="default" || method=="invit") { // Stodola's inverse iteration with Gramm-Schmidt orthogonalization
		var DM = highest? mat.inv().mulm(this) : this.inv().mulm(mat)
		var eigVals = [];
		var eigVecs = [];
		var eigVecsMulmM = [];
		var rho, rhoRatio, rhoNew, iter;
		var xk = JSMatrix.Vector.Zeros(this.nRows);
		for (var k=0; k<nEigModes; k++) {
			xk.one();
			iter = 0;
			rhoNew = 0.;
			do {
				xk = DM.mulv(xk);
				for (var j=0; j<k; j++) { xk.isub(eigVecs[j].mulf(eigVecsMulmM[j].dot(xk))); } // Gramm-Schmidt orthogonalization
				xk.imulf(1./xk.energyNorm(mat));
				rho = rhoNew;
				rhoNew = xk.squaredEnergyNorm(this) / xk.squaredEnergyNorm(mat)
				rhoRatio = Math.abs((rhoNew-rho)/rhoNew);
				iter++;
				if (iter > maxiter) {
					if (rhoRatio > JSMatrix.TOL) {
						/*jsmLogger.warn('Matrix.eigSymmetric: maximum number of iterations reached');*/
						return null;
					}
					break;
				}
			} while (rhoRatio > JSMatrix.TOL2)
			eigVals[k] = rhoNew;
			eigVecs[k] = xk.copy();
			eigVecsMulmM[k] = xk.mulm(mat);
		}
		return {eigVals:eigVals,eigVecs:eigVecs};
	}
	if (method == "subspace") { // subspace iteration
		var eigVals = [];
		var eigVecs = [];
		var n = this.nRows;
		var kk, mm, r, q;
		var rhos, rhosRatio, rhoNews, maxRhosRatio, iter;
		var xk = JSMatrix.Vector.Zeros(this.nRows);
		var x = JSMatrix.Matrix.Ones(n,nEigModes);
		for (var i=0; i<n; i++) { x.set(i,i,0.); }
		var lu = highest? mat.luDecomposition() : this.luDecomposition()
		var rhosNew = [];
		for (n=0; n<nEigModes; n++) { rhosNew[n] = 0.; }
		rhosRatio = rhosNew.slice();
		iter = 0;
		do {
			x = highest? mat.linSolve(this.mulm(x),"lu",false,lu) : this.linSolve(mat.mulm(x),"lu",false,lu)
			kk = x.transposed().mulm(this).mulm(x);
			mm = x.transposed().mulm(mat).mulm(x);
			r = kk.eigSymmetric({mat:mm,nEigModes:nEigModes,highest:highest});
			q = JSMatrix.Matrix.Zeros(nEigModes,nEigModes);
			for (var j=0; j<nEigModes; j++) {
				q.setCol(j,r.eigVecs[j]);
			}
			x = x.mulm(q);
			rhos = rhosNew.slice();
			rhosNew = r.eigVals.slice();
			for (n=0; n<nEigModes; n++) { rhosRatio[n] = Math.abs((rhosNew[n] - rhos[n])/rhosNew[n]); }
			maxRhosRatio = 0.;
			for (n=0; n<nEigModes; n++) {
				if (rhosRatio[n] > maxRhosRatio) { maxRhosRatio = rhosRatio[n]; }
			}
			iter++;
			if (iter > maxiter) {
				if (maxRhosRatio > JSMatrix.TOL) {
					return null;
					/*jsmLogger.warn('Matrix.eigSymmetric: maximum number of iterations reached');*/
				}
				break;
			}
		} while (maxRhosRatio > JSMatrix.TOL2)
		for (var i=0; i<nEigModes; i++) {
			eigVals.push(r.eigVals[i]);
			eigVecs.push(x.getCol(i));
		}
		return {eigVals:eigVals,eigVecs:eigVecs};
	}
	params.method = "default";
	return this.eigSymmetric(params);
};

/** Returns eigenvalues of receiver (this*x = ret*x)
* @param {string=} method ="default"] method of solution. "default" and "schur" stand for similar transformation of receiver into schur form
* @returns {Array.<number>} eigenvalues of receiver
* @example
* m1 = JSMatrix.mat([[4,2,3,1],[3,6,5,7],[5,3,8,1],[2,3,4,5]]);
* v1 = m1.schurForm().diag().toArray();
* // v1 = [15.6072, 4.1916, 2.2012, 1 ]
* v2 = m1.eigVals();
* // v2 = [15.6072, 4.1916, 2.2012, 1 ]
*/
JSMatrix.Matrix.prototype.eigVals = function(method) {
	method = method || "default";
	method = method.toLowerCase();
	if (method=="default" || method=="schur") {
		return this.schurForm().diagonalToVector().toArray();
	}
	return this.eigVals("default");
};

/** Returns orthogonal matrix P and diagonal matrix L of eigendecomposition of the receiver such that this=P*L*P^-1 (obtained from eigenvalue equation this*P = P*L). If receiver is symmetric, P is orthogonal. Columns of P are eigenvectors of this and diagonal elements of L are corresponding eigenvalues
* @param {string=} method (default="default") method of eigenvalues solution ["defaul","symm"]
* @param {{q:JSMatrix.Matrix,s:JSMatrix.Matrix}=} precompSchurDecomp precomputed Schur decomposition
* @returns {{p:JSMatrix.Matrix,l:JSMatrix.Matrix}|null} {p,l} matrix P and diagonal matrix L of eigendecomposition of receiver
* @example
* m1 = JSMatrix.mat([[4,2,3,1],[3,6,5,7],[5,3,8,1],[2,3,4,5]]);
* pl = m1.eigenDecomposition()
* p = pl.p
* l = pl.l
* // p = Matrix([[0.287,-0.05676,-0.7041, 0.4897],[0.6779,0.7997,-1.367,-0.9794],[0.5142,-0.6248,1.294,0],[0.4398,0.2640,0.1196,0.4897]])
* // l = Matrix([[15.61,0,0,0],[0,4.192,0,0],[0,0,2.201,0],[0,0,0,1]])
* c = p.mulm(l).mulm(p.inv())
* // c = 
*/
JSMatrix.Matrix.prototype.eigenDecomposition = function(method,precompSchurDecomp) {
	if (!this.isSquare) { /*jsmLogger.warn('Matrix.eigenDecomposition: matrix is not square');*/ return null; }
	method = method || "default"
	if (method=="default") {
		var i, j, n=this.nRows;
		var sd = precompSchurDecomp? precompSchurDecomp : this.schurDecomposition()
		var q = sd.q;
		var s = sd.s;
		var l = new JSMatrix.Matrix(n,n);
		for (i=0; i<n; i++) { l[i][i] = s[i][i]; }
		var p = JSMatrix.Matrix.Zeros(n,n);
		var t11 = new JSMatrix.Matrix(), nn, lambda, v = new JSMatrix.Vector(), y = new JSMatrix.Vector(n);
		p.setCol(0,q.getCol(0));
		var ii
		for (i=1; i<n; i++) {
			ii = JSMatrix.range(i);
			t11.beSubMatrixOf(s,ii,ii);
			lambda = s[i][i];
			v.resize(i);
			for (j=0; j<i; j++) {
				t11[j][j] -= lambda;
				v[j] = s[j][i];
			}
			y.beSolutionOf(t11,v);
			y.negate();
			y.appendElements(1.);
			y.resize(n);
			p.setCol(i,q.mulv(y));				
		}
		return {p:p,l:l};
	}
	if (method=="symm" || method=="symmeig" || method=="eigsymm" || method=="eigsymmetric") {
		var n = this.nRows;
		var e = this.eigSymmetric();
		var p = JSMatrix.Matrix.Zeros(n,n);
		var l = JSMatrix.Matrix.Zeros(n,n);
		for (var i=0; i<n; i++) {
			p.setCol(i,e.eigVecs[i]);
			l[i][i] = e.eigVals[i];
		}
		return {p:p,l:l};
	}
	return this.eigenDecomposition("default")
};


/** Returns spectral decomposition of receiver, identical to {@JSMatrix.Matrix#eigenDecomposition}
* @function
*/
JSMatrix.Matrix.prototype.spectralDecomposition = JSMatrix.Matrix.prototype.eigenDecomposition;

/* * Returns generalized eigenvalues of receiver (this*x = ret*mat*x)
* @param {JSMatrix.Matrix} mat matrix in generalized eigenvalue equation this*x = ret*mat*x
* @param {string=} method ="default"] method of solution. "default and "inv" stand for solution mat-1*this*x = ret*x (standard eigenvalue problem), "qz" for qz decomposition
* @returns {Array.<number>} eigenvalues of receiver
* @ example
* /
generalizedEigVals: function(mat,method) {
	var method = method || "default";
	method = method.toLowerCase();
	if (method=="default" || method=="inv") {
		return mat.inv().mulm(this).eigVals();
	}
	if (method=="qz" || method=="qzdecomposition") {
		var qz = this.qzDecomposition(mat);
		//TODO
		return null;
	}
	return this.generalizedEigVals(mat);
};
*/

/** Returns inversion of receiver
* @param {string=} method ="default"] method used for solution. Options are "default" and "gaussjordan" for Gauss-Jordan method, "lu" for solution with LU decomposition, "eigen" with eigen decomposition, "symmeigen" with eigen decomposition when receiver is symmetric, "svd" for SVD decomposition
* @param {Array.<JSMatrix.Matrix>=} precompDecomps =undefined] precomputed matrix decompositions if we want to use some
* @returns {JSMatrix.Matrix} inversion of receiver
* @example
* a = JSMatrix.mat([[1,2,2],[1,0,1],[1,2,1]]);
* i1 = a.inversed(); // i1 =
* c = a.mulm(i1) // c = Matrix([[1,0,0],[0,1,0],[0,0,1]])
* i2 = a.inversed('lu') // i2 = Matrix([[-1,1,1],[0,-0.5,0.5],[1,0,-1]])
* i3 = a.inversed('eigen') // i3 = Matrix([[-1,1,1],[0,-0.5,0.5],[1,0,-1]])
* i4 = a.inversed('svd') // i4 = Matrix([[-1,1,1],[0,-0.5,0.5],[1,0,-1]])
* check = a.x(i1); // check = Matrix([[1,0,0],[0,1,0],[0,0,1]])
*/
JSMatrix.Matrix.prototype.inversed = function(method,precompDecomps) {
	if (!this.isSquare) { /*jsmLogger.warn('Matrix.inversed: matrix is not square');*/ return null; }
	method = method || "default";
	method = method.toLowerCase();
	if (method=="default" || method=="gaussjordan") {
		/* solves this*ret = I with Gauss-Jordan ellimination*/
		return this.gaussJordanElimination(JSMatrix.Matrix.Identity(this.nRows));
	}
	if (method=="lu") {
		/* solves this*ret[:,i] = I[:,i] for i=0,1,..,this.nRows with the help of LU decomposition*/
		var lu = precompDecomps? precompDecomps : this.luDecomposition();
		var n = this.nRows;
		var ret = new JSMatrix.Matrix(n,n);
		var vec = JSMatrix.Vector.Zeros(n);
		for (var i=0; i<n; i++) {
			vec.zero();
			vec[i] = 1.;
			this.linSolve(vec,"lu",false,lu);
			ret.setCol(i,vec)
		}
		return ret;
	}
	if (method=="chol" || method=='cholesky') {
		/* solves this*ret[:,i] = I[:,i] for i=0,1,..,this.nRows with the help of Cholesky decomposition*/
		var l = precompDecomps? precompDecomps : this.choleskyDecomposition();
		if (l==null) { return null; }
		var n = this.nRows;
		var ret = new JSMatrix.Matrix(n,n);
		var vec = JSMatrix.Vector.Zeros(n);
		for (var i=0; i<n; i++) {
			vec.zero();
			vec[i] = 1.;
			this.linSolve(vec,"cholesky",false,l);
			ret.setCol(i,vec)
		}
		return ret;
	}
	if (method=="eigen" || method=="eig") {
		/* this = V*L*V^-1 -> this^-1=V*L^-1*V^-1. As L is diagonal, the inversion is easy to compute*/
		var e = precompDecomps? precompDecomps : this.eigenDecomposition();
		var p = e.p;
		var l = e.l;
		for (var i=0; i<this.nRows; i++) { l[i][i] = 1./l[i][i]; }
		return p.mulm(l).mulm(p.inversed());
	}
	if (method=="symmeig" || method=="eigsymm" || method=="eigsymmetric") {
		var e = precompDecomps? precompDecomps : this.eigenDecomposition("symm");
		var p = e.p;
		var l = e.l;
		for (var i=0; i<this.nRows; i++) { l[i][i] = 1./l[i][i]; }
		return p.mulm(l).mulm(p.transposed());
	}
	if (method=="svd") {
		var svd = precompDecomps? precompDecomps : this.singularValueDecomposition();
		var u=svd.u, s=svd.s, v=svd.v;
		var n = s.nElements;
		var ss = new JSMatrix.Matrix(n,n);
		for (var i=0; i<n; i++) { if (s[i] > JSMatrix.TOL) { ss[i][i] = 1./s[i]; } }
		return v.mulm(ss).mulm(u.transposed());
	}
	return this.inversed("default");
};

/** Alias for inversion, see {@JSMatrix.Matrix#inversed}
* @function
*/
JSMatrix.Matrix.prototype.inv = JSMatrix.Matrix.prototype.inversed;

/** Sets receiver as an inversion of given matrix
* @param {JSMatrix.Matrix} mat given matrix
* @param {string=} method ="default"] see {@JSMatrix.Matrix#inversed} for details}
* @param {Array.<JSMatrix.Matrix>=} precompDecomps =undefined] precomputed matrix decompositions if we want to use some
* @example
* m1 = JSMatrix.mat([[1],[2]]);
* m2 = JSMatrix.mat([[1,2,2],[1,0,1],[1,2,1]]);
* m1.beInverseOf(m2); // m1 = Matrix([[2.3,-1.4,-0.8],[-1.4,1.2,0.4],[-0.8,0.4,0.8]])
*/
JSMatrix.Matrix.prototype.beInverseOf = function(mat,method,precompDecomps) {
	if (!mat.isSquare()) { /*jsmLogger.warn('Matrix.beInverseOf: matrix is not square');*/ return; }
	method = method || "default";
	method = method.toLowerCase();
	if (method=="default" || method=="gaussjordan") {
		var temp = mat.copy();
		this.beUnitMatrix(mat.nRows);
		temp.gaussJordanElimination(this,"nopivot",false)
		return
	} else { // TODO
		this.beCopyOf(mat.inversed(method,precompDecomps));
	}
	this.beInverseOf(mat,"default",precompDecomps);
};

/** Returns determinant of receiver.
* @param {string=} method ="default"] used method. Options are "default" and "lu" for using LU decomposition, "schurForm" or "schur" for using schur form, "schurDecomposition" for using schur decomposition, "eigen" or "eigenDecomposition" for using eigendecomposition
* @param {Object=} precompDecomps precomputed decomposition (consistent with method parameter)
* @returns {number|null} determinant of receiver
* @example
* a = JSMatrix.mat([[4,2,1],[2,5,3],[1,3,6]]);
* d = a.determinant(); // d = 67
* d = a.determinant('lu'); // d = 67
* d = a.determinant('schur'); // d = 67
* d = a.determinant('schurdecomposition'); // d = 67
* d = a.determinant('eigen'); // d = 67
*/
JSMatrix.Matrix.prototype.determinant = function(method,precompDecomps) {
	if (!this.isSquare()) { /*jsmLogger.warn('Matrix.determinant: matrix is not square');*/ return null; }
	method = method || "default"
	method = method.toLowerCase();
	var mats = precompDecomps || undefined
	if (method=="default" || method=="lu") {
		// transform matrix to lower and upper triangular. Determinant is then product of diagonal elements. As L has 1. on diagonal, det is product of diagnal elements of U only
		var u = precompDecomps? precompDecomps[1] : this.luDecomposition().u;
		return u.diagProduct();
	}
	if (method=="schurform" || method=="schur") {
		// create schur form (upper triangular matrix with eigenvalues on diagonal). It is triangular, so det is product of diagonal elements (as well as product of eigenvalues)
		var s = precompDecomps? precompDecomps : this.schurForm();
		return s.diagProduct();
	}
	if (method=="schurdecomposition") {
		// see case:"schurform"
		var s = precompDecomps? precompDecomps : this.schurDecomposition().s;
		return s.diagProduct();
	}
	if (method=="eigen" || method=="eigendecomposition") {
		var l = precompDecomps? precompDecomps.l : this.eigenDecomposition().l
		return l.diagProduct();
	}
	return this.determinant("default",precompDecomps);
};

/** Alias for determinant, see {@JSMatrix.Matrix#determinant}
* @function
*/
JSMatrix.Matrix.prototype.det = JSMatrix.Matrix.prototype.determinant;

/** Alias for determinant, see {@JSMatrix.Matrix#determinant}
* @function
*/
JSMatrix.Matrix.prototype.giveDeterminant = JSMatrix.Matrix.prototype.determinant;

/** Constructs new Matrix from given 2D array
* @param {Array.<Array.<number>>} arry (default=[[]]) array containing elements of new matrix
* @returns {JSMatrix.Matrix} new Matrix object
* @example
* m = JSMatrix.Matrix.create([[1,2,3],[4,5,6],[7,8,9]]) // m = Matrix([[1,2,3],[4,5,6],[7,8,9]])
*/
JSMatrix.Matrix.create = function(arry) {
	if (arry==undefined) { return new JSMatrix.Matrix(); }
	if (arry.length==0) { return new JSMatrix.Matrix(); }
	var a = arry;
	var nRows = a.length, nCols = a[0].length, r, c;
	var ret = new JSMatrix.Matrix();
	for (r=0; r<nRows; r++) {
		ret[r] = {};
		for (c=0; c<nCols; c++) {
			ret[r][c] = arry[r][c];
		}
	}
	ret.nRows = nRows;
	ret.nCols = nCols;
	return ret;
};

/** Constructs new diagonal matrix from given 1D array
* @param {Array.<number>|JSMatrix.Vector} arry array or vector containing elements of new matrix
* @returns {JSMatrix.Matrix} new Matrix object 
* @example
* m1 = JSMatrix.Matrix.Diagonal([1,2,3]); // m1 = Matrix([[1,0,0],[0,2,0],[0,0,3]])
* m2 = JSMatrix.Matrix.Diagonal(JSMatrix.vec([1,2,3])); // m2 = Matrix([[1,0,0],[0,2,0],[0,0,3]])
*/
JSMatrix.Matrix.Diagonal = function(arry) {
	var s = arry.nElements==undefined? arry.length : arry.nElements;
	var ret = new JSMatrix.Matrix(s,s);
	for (var i=0; i<s; i++) {
		ret[i][i] = arry[i];
	}
	return ret;
};

/** Constructs new identity matrix of given size
* @param {number} size size of returned matrix
* @returns {JSMatrix.Matrix} identity matri
* @example
* m = JSMatrix.Matrix.Identity(4) // m = Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
*/
JSMatrix.Matrix.Identity = function(size) {
	var ret = new JSMatrix.Matrix(size,size);
	for (var i=0; i<size; i++) {
		ret[i][i] = 1.;
	}
	return ret;
};

/** Alias for JSMatrix.Matrix.Identity(), see {@JSMatrix.Matrix#Identity}
* @function
*/
JSMatrix.Matrix.I =  JSMatrix.Matrix.Identity;

/** Alias for Matrix.Identity(), see {@JSMatrix.Matrix#Identity}
* @function
*/
JSMatrix.Matrix.UnitMatrix = JSMatrix.Matrix.Identity;

/** Creates a matrix of given size full of zeros
* @param {number=} nRows =0] number of rows
* @param {number=} nCols =0] number of columns
* @returns {JSMatrix.Matrix} new matrix full of zero
* @example
* m = JSMatrix.Matrix.Zeros(2,3); // m = Matrix([[0,0,0],[0,0,0]])
*/
JSMatrix.Matrix.Zeros = function(nRows,nCols) {
	return new JSMatrix.Matrix(nRows,nCols);
};

/** Creates a matrix of given size full of ones
* @param {number=} nRows =0] number of rows
* @param {number=} nCols =0] number of columns
* @returns {JSMatrix.Matrix} new vector full of zero
* @example
* m = JSMatrix.Matrix.Ones(2,3); // m = Matrix([[1,1,1],[1,1,1]])
*/
JSMatrix.Matrix.Ones = function(nRows,nCols) {
	nRows = nRows || 0;
	nCols = nCols || 0;
	var ret = new JSMatrix.Matrix();
	for (var r=0; r<nRows; r++) {
		ret[r] = {};
		for (var c=0; c<nCols; c++) {
			ret[r][c] = 1.;
		}
	}
	ret.nRows = nRows;
	ret.nCols = nCols;
	return ret;
};

/** Creates a Householder matrix from given vector (ret = I - 2*v*v^T/(v^T*v). Householder matrix is symmetric, orthogonal and involutary
* @param {JSMatrix.Vector} vec=0 Householder vector to form Householder matrix
* @param {JSMatrix.Matrix=} precompI =undefined] precomputed identity matrix to save some time of creation
* @returns {JSMatrix.Matrix} Householder matrix
* @example
* m1 = JSMatrix.Matrix.Householder(JSMatrix.vec([0,1,0]))
* // m1 = Matrix([[1,0,0],[0,-1,0],[0,0,1]])
* b1 = m1.isSymmetric() // b1 = true
* b2 = m1.isOrthogonal() // b2 = true
* b3 = m1.isInvolutary() // b3 = true
*/
JSMatrix.Matrix.Householder = function(vec,precompI) {
	var nElements = vec.nElements;
	return (precompI? precompI : JSMatrix.Matrix.Identity(nElements)).sub(vec.dyadic(vec).mulf(2./vec.dot(vec)));
};
/*
* m2 = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]);
* v = m2.getCol(0);
* q = JSMatrix.Matrix.Householder(v.sub(JSMatrix.vec([1,0,0]).mulf(v.norm())));
* m3 = q.mulm(m2);
* // m3 = 
*/








/* **************************************************
*
* Initializing functions (shotcuts constructors)
*
***************************************************/

/** Returns new Matrix object from given array
* @param {Array.<Array.<number>>} arry (default=[[]]) 2D array containing matrix elements
* @returns {JSMatrix.Matrix} new Matrix object
* @example
* m = JSMatrix.mat([[1,2,3],[4,5,6],[7,8,9]]) // m = Matrix([[1,2,3],[4,5,6],[7,8,9]])
*/
JSMatrix.mat = JSMatrix.Matrix.create;

/** Returns new Vector object from given array
* @param {Array.<number>} arry (default=[]) array containing vector elements
* @returns {JSMatrix.Vector} new Vector object
* @example
* v = JSMatrix.vec([1,2,3,4]) // v = Vector([1,2,3,4])
*/
JSMatrix.vec = JSMatrix.Vector.create;

/** Returns new Vector or Matrix full of given size full of zeros
* @param {number} nRows number of rows of returned object
* @param {number=} nCols ] if specified, new Matrix of size (nRows,nCols) is returned, new vector of size nRows otherwise
* @returns {JSMatrix.Matrix|JSMatrix.Vector} new Matrix or Vector object
* @example
* m = JSMatrix.zeros(2,3); // m = Matrix([[0,0,0],[0,0,0]])
* v = JSMatrix.zeros(4) // v = Vector([0,0,0,0])
*/
JSMatrix.zeros = function(nRows,nCols) {
	if (nCols) { return JSMatrix.Matrix.Zeros(nRows,nCols); }
	if (nRows) { return JSMatrix.Vector.Zeros(nRows); }
	return null;
};

/** Returns new Vector or Matrix full of given size full of ones
* @param {number} nRows number of rows of returned object
* @param {number=} nCols ] if specified, new Matrix of size (nRows,nCols) is returned, new vector of size nRows otherwise
* @returns {JSMatrix.Matrix|JSMatrix.Vector} new Matrix or Vector object
* @example
* m = JSMatrix.ones(2,3); // m = Matrix([[1,1,1],[1,1,1]])
* v = JSMatrix.ones(4) // v = Vector([1,1,1,1])
*/
JSMatrix.ones = function(nRows,nCols) {
	if (nCols) { return JSMatrix.Matrix.Ones(nRows,nCols); }
	if (nRows) { return JSMatrix.Vector.Ones(nRows); }
	return null;
};

/** Returned identity matrix of given size, see {@JSMatrix.Matrix#Identity}
* @param {number} size number of rows and columns of returned matrix
* @returns {JSMatrix.Matrix} nRows x nRows identity matrix
* @example
* m = JSMatrix.identity(3); // m = Matrix([[1,0,0],[0,1,0],[0,0,1]])
*/
JSMatrix.identity = JSMatrix.Matrix.Identity;

/** Alias for dentity(), see {@link identity}
*/
JSMatrix.eye = JSMatrix.Matrix.Identity;

/** Alias for dentity(), see {@link identity}
*/
JSMatrix.unitMatrix = JSMatrix.Matrix.Identity;


/* ***********************************************************
*
* Other functions
*
************************************************************/

/** Linear system equation solver. Returns vector x as a solution of a*ret=rhs, see {@JSMatrix.Matrix#linSolve} for input desription
* @param {JSMatrix.Matrix} a
* @param {JSMatrix.Vector|Array.<JSMatrix.Vector>|JSMatrix.Matrix} rhs
* @param {string=} method ="default"]
* @param {boolean=} saveOrig =true]
* @param {Array.<JSMatrix.Matrix>=} precompDecomps ]
* @returns {JSMatrix.Vector|Array.<JSMatrix.Vector>|JSMatrix.Matrix} solution x of a*x=rhs
*/
JSMatrix.linSolve = function(a,rhs,method,saveOrig,precompDecomps) {
	return a.linSolve(rhs,method,saveOrig,precompDecomps);
};


/** Eigenvalues and eigenvectors solver for symmetric matrices. Returns first nEigModes eigenvales (lambda) and eigenvectors (phi) of problem mat1*phi=lambda*mat2*phi. See {@link JSMatrix.Matrix#eigSymmetric} for input description (except mat1 and mat2). Eigenvectors are normalized with respect to mat2
* @param {{
				mat1: JSMatrix.Matrix,
				mat2: JSMatrix.Matrix,
				nEigModes: number,
				method: string,
				highest: boolean,
				maxiter: number}=} params
* @param {JSMatrix.Matrix} params.mat1 mat1 in above equation
* @param {JSMatrix.Matrix} params.mat2 mat2 in above equation, mat in {@link JSMatrix.Matrix#eigSymmetric}
* @returns {{eigVals:Array.<number>,eigVecs:Array.<JSMatrix.Vector>}} [[lambda1,lambda2,...,lambdaN],[phi1,phi2,...,phiN]], lambdai is i-th eigenvale, phii is i-th eigenvector. N = nEigModes
*/
JSMatrix.eigSymmetric = function(params) {
	return params.mat1.eigSymmetric({
		mat:params.mat2,
		nEigModes:params.nEigModes,
		method:params.method,
		highest:params.highest,
		maxiter:params.maxiter
	});
};









/* ***********************************************************
*
* Array extending functions functions
*
************************************************************/

/** Returns sequence (array) with start, stop, step. Inspired by <a href="www.python.org">Python</a> syntax
* @param {number} a start/stop (according to b, see example)
* @param {number=} b stop
* @param {number=} c (default=1) step
* @returns {Array.<number>} array containing given sequence
* @example
* r = JSMatrix.range(8);     // r = [0,1,2,3,4,5,6,7]
* r = JSMatrix.range(2,8);   // r = [2,3,4,5,7]
* r = JSMatrix.range(2,8,2); // r = [2,4,6]
*/
JSMatrix.range = function(a,b,c) {
	var start;
	var step;
	var stop;
	if (b===undefined && c===undefined) {
		start = 0;
		stop  = a;
		step  = 1;
	} else if (c===undefined) {
		start = a;
		stop  = b;
		step  = 1;
	} else {
		start = a;
		stop  = b;
		step  = c;
	}
	var ret = []
	var val = start;
	while (val < stop) {
		ret.push(val);
		val += step;
	}
	return ret;
};

/**@ignore*/
if (Array.prototype.indexOf === undefined) {
	Array.prototype.indexOf = function(val) {
		for (var i=0; i<this.length; i++) {
			if (this[i] === val) {
				return i;
			}
		}
		return -1;
	};
}

/** Returns elements of receiver that are not in B (Matlab-like function), ret = this - B
* @param {Array} B
* @returns {Array} array of elements of receiver that are not in B
* @example
* a = [1,2,3,4,5,6];
* b = [2,5,1,99,54];
* c = a.setdiff(b) // c = [3,4,6]
*/
Array.prototype.setdiff = function(B) {
	var ret = [];
	for (var a=0; a<this.length; a++) {
		var noMatch = 1;
		for (var b=0; b<B.length; b++) {
			if (this[a] == B[b]) { noMatch = 0; break; }
		}
		if (noMatch) { ret.push(this[a]); }
	}
	return ret;
}
