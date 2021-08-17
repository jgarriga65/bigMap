/*
my licensing information
*/

#ifndef QTREE_H
#define QTREE_H

using namespace std;

// ++++

// +++ new-node-center placement factor (dx:cos(pi/4), dy:sin(pi/4))
// so: deltaCFactor = cos(pi/4) /4
// where factor 1/4 comes from: dc = node.sz /4 *cos(pi/4)
static const double deltaCFactor = .176776695;

// +++ maximum tree depth levels
// in case of duplicate data-points or dist(p, q) extremly small
// the tree would expand indefinitely trying to fit each single point in a different leaf
// We avoid this by setting a maximum depth value.
static const int treeMaxDepth = 12;

// ++++

struct Point{
	double x;                   // x coordinate
	double y;                   // y coordinate
};

static void showP(Point p) {
	cout << p.x << "," << p.y << endl;
};

static double sqDist(Point p, Point q) {
	return std::pow(p.x -q.x, 2) + std::pow(p.y -q.y, 2);
};

// +++

class Node {

public:

	double sz;                    // cell-size
	Point c;                      // cell-center
	Point m;                      // cell-mass-center
	Point p;                      // data-point in cell

	int Np;                       // number of data-points in cell
	std::vector<Node*> child;     // child nodes

public:

	Node(){};
	Node(double x_center, double y_center, double cell_size) {
		sz = cell_size;
		c = {x_center, y_center};
		m = {.0, .0};
		Np = 0;
	};
	~Node(){
		for (size_t k = 0; k < child.size(); k++) delete child[k];
	};

	string childName(int k) const {
		if (k == 0) return "upperLeft";
		if (k == 1) return "upperRight";
		if (k == 2) return "lowerRight";
		if (k == 3) return "lowerLeft";
		return "";
	};

	void show() const {
		if (Np != 1){
			cout << "    size: " << sz << ", #points: " << Np << endl;
		} else {
			cout << "    size:" << sz << ", (" << p.x << ", " << p.y << ")" << endl;
		}
		cout << "  center: " << c.x << ", " << c.y << endl;
		if (Np > 0){
			cout << "  c.mass: " << m.x /Np << ", " << m.y /Np << endl;
		} else {
			cout << "  c.mass: " << m.x << ", " << m.y << endl;
		}
	};

	void updM(Point q) {
		// update cell-mass center
		m.x = m.x *Np +q.x;
		m.y = m.y *Np +q.y;
		Np ++;
		m.x /= Np;
		m.y /= Np;
	};

	void split() {
		double childNodeSz = sz /2.0;
		double deltaC = sz * deltaCFactor;
		// upper-left child
		child.push_back(new Node(c.x -deltaC, c.y +deltaC, childNodeSz));
		// upper-right child
		child.push_back(new Node(c.x +deltaC, c.y +deltaC, childNodeSz));
		// lower-right child
		child.push_back(new Node(c.x +deltaC, c.y -deltaC, childNodeSz));
		// lower-left child
		child.push_back(new Node(c.x -deltaC, c.y -deltaC, childNodeSz));
		// if splitting it must be a non-empty leaf, thus push p
		int k = whereIs(p);
		child[k]->p = p;
		child[k]->updM(p);
	};

	int whereIs(Point q) const {
		if (q.x <= c.x & q.y >  c.y) return 0;
		if (q.x >  c.x & q.y >  c.y) return 1;
		if (q.x >  c.x & q.y <= c.y) return 2;
		if (q.x <= c.x & q.y <= c.y) return 3;
		// in case q = {nan, nan} (should never happen)
		return -1;
	};

	bool prune(Point q, double delta) const {
		return sz /(2 *sqrt(sqDist(m, q))) < delta;
	};

};

// ++++
class Quadtree {

public:

	size_t n;                   // number of data-points
	size_t m;                   // dataset dimension
	Node* root;                 // root node

	double minNodeSize;         // node minimum size
								// (avoid indefinite expansion in case of duplicates)

	double delta; 				// pt-SNE approximation factor (theta)

public:

	Quadtree();
	Quadtree(double* M, size_t rows, size_t cols, double theta) {

		n = rows;
		m = cols;
		delta = theta;

		double minx = .0, miny = .0, maxx = .0, maxy = .0;
		for (size_t i = 0; i < n; i++) {
			minx = std::min(minx, M[i *m +0]);
			maxx = std::max(maxx, M[i *m +0]);
			miny = std::min(miny, M[i *m +1]);
			maxy = std::max(maxy, M[i *m +1]);
		}
		double xc = (minx + maxx) /2.0;
		double yc = (miny + maxy) /2.0;
		double sz = std::sqrt(std::pow((maxx -minx), 2) +std::pow((maxy -miny), 2));

		root = new Node(xc, yc, sz);
		// limit the tree size up to treeMaxDepth levels
		minNodeSize = sz /std::pow(2, treeMaxDepth);

		// build the tree
		for (size_t i = 0; i < n; i++){
			Point q = {M[i *m +0], M[i *m +1]};
			if (root->whereIs(q) != -1) {
				insertP(root, q);
			}
		}
	};

	~Quadtree(){
		delete root;
	};

	void show(int depth) const {
		show(root, depth, -1, "root");
	};

	void insertP(Node* node, Point q) {
		// check if node is empty
		if (node->Np == 0 or node->sz < minNodeSize) {
			// node is empty or at treeMaxDepth
			// (in the latter case node->p is irrelevant)
			node->p = {q.x, q.y};
		} else {
			// if leaf split (and push p)
			if (node->child.empty()) node->split();
			// insert q
			insertP(node->child[node->whereIs(q)], q);
		}
		// update mass-center
		node->updM(q);
	};

	void repForces(double* Y, double* qForce, double* Q) const {
		for (size_t i = 0; i < n; i++) {
			Point q = {Y[i *m +0], Y[i *m +1]};
			repForce(q, root, qForce +i *m, Q);
		}
	};

	void repF(double* Y, double* qForce, double* Q, size_t i) const {
		Point q = {Y[i *m +0], Y[i *m +1]};
		qForce[0] = .0;
		qForce[1] = .0;
		repForce(q, root, qForce, Q);
	};

private:

	void show(Node* node, int depth, int lvl, string nodeName) const {
		if (lvl == depth -1 || lvl == depth) {
			cout << "+++ lvl. " << lvl << "." << nodeName << endl;
			node->show();
		}
		if (lvl < depth) {
			lvl ++;
			if (!node->child.empty()) {
				for (size_t k = 0; k < node->child.size(); k++) {
					show(node->child[k], depth, lvl, node->childName(k));
				}
			}
		}
	};

	void repForce(Point q, Node* node, double qForce[], double* Q) const {
		for (size_t k = 0; k < 4; k++) {
			if (node->child[k]->Np > 0) {
				if (node->child[k]->child.empty()) {
					// child-node is leaf
					// (Att.!! can hold more than one single point if at treeMaxDepth)
					addForce(q, node->child[k]->m, node->child[k]->Np, qForce, Q);
				} else if (node->child[k]->prune(q, this ->delta)) {
					// branch can be prunned: approximate repForce by its mass-center
					addForce(q, node->child[k]->m, node->child[k]->Np, qForce, Q);
				} else {
					// branch can not be prunned: explore branch deeper
					repForce(q, node->child[k], qForce, Q);
				}
			}
		}
	};

	void addForce(Point q, Point p, int Np, double qForce[], double* Q) const {
		double Qqp = 1.0 /(1.0 + sqDist(q, p));
		qForce[0] += Qqp *Qqp *(q.x - p.x) *Np;
		qForce[1] += Qqp *Qqp *(q.y - p.y) *Np;
		// update Q-normalization factor
		*Q += Qqp *Np;
	};

};

#endif
