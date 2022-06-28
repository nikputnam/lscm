#include "LSCM.h"
#include <assert.h>
#include <math.h>
#include <float.h>
#include <Eigen/Sparse>
#include <iostream>
#include <sstream>
#include <Eigen/IterativeLinearSolvers>

using namespace MeshLib;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::VectorXd VectorXd;

LSCM::LSCM(Mesh * mesh) {
	m_mesh = mesh;
}

LSCM::~LSCM(){}

void LSCM::set_coefficients() {
	for (MeshEdgeIterator eiter(m_mesh); !eiter.end(); ++eiter) {
		Edge * e = *eiter;
		e_l(e) = m_mesh->edge_length(e);
	}

	for (MeshFaceIterator fiter(m_mesh); !fiter.end(); ++fiter) {
		Point  p[3];
		Face *f = *fiter;

		double l[3];
		HalfEdge * he = f->halfedge();
		for (int j = 0; j < 3; j++) {
			Edge * e = he->edge();
			l[j] = e_l(e);
			he = he->he_next();
		}

		double a = acos((l[0]*l[0] + l[2]*l[2] - l[1]*l[1]) / (2*l[0]*l[2]));

		p[0] = Point(0, 0, 0);
		p[1] = Point(l[0], 0, 0);
		p[2] = Point(l[2]*cos(a), l[2]*sin(a), 0);

		Point n = (p[1]-p[0]) ^ (p[2]-p[0]);
		double area = n.norm() / 2.0;
		n /= area;

		he = f->halfedge();
		for (int j = 0; j < 3; j++) {
			Point s = (n ^ (p[(j + 1) % 3] - p[j])) / sqrt(area);
			c_s(he) = s;
			he = he->he_next();
		}
	}
}

void LSCM::project() {
	set_coefficients();

	std::vector<Vertex*> vertices;
    std::vector<Face*> faces;

	int n_zip=0;
	int n_extra_constrants=0;

	std::map<int,int> zid2count;
	std::map<int,std::pair<Vertex*,Vertex*>> zid2vertexPair;
	std::map<int,std::pair<int,int>> zid2vertexIdPair;

	for (MeshVertexIterator viter(m_mesh); !viter.end(); ++viter){
		Vertex * v = *viter;
		if (v->string().substr(0,3) == "fix") {
			m_fix_vertices.push_back(v);
		
		} else {

			if (v->string().substr(0,3) == "zip") {
				n_zip++;

				std::string tmp;
				int zid;
				std::stringstream(v->string()) >> tmp >> zid;

				if ( zid2count.find(zid) == zid2count.end() ) { zid2count[zid]=0; } zid2count[zid]+=1;
				if ( zid2vertexPair.find(zid) == zid2vertexPair.end() ) { 
					zid2vertexPair[zid]=std::pair<Vertex*,Vertex*>(v,0); 
				} else {
					zid2vertexPair[zid]=std::pair<Vertex*,Vertex*>(zid2vertexPair[zid].first,v);
				}

				std::cout << "zip: " << zid << "\n";
			}
			vertices.push_back(v);
		}
	}
	assert(m_fix_vertices.size()>=2);
	std::cout << "n zip: " << n_zip << "\n";

	for (int k = 0; k < (int)vertices.size(); k++ ){
		v_idx(vertices[k]) = k;
	}

    // build zid2vertexIdPair
	for(auto it=zid2vertexPair.begin();it!=zid2vertexPair.end();it++) {
		int zid = it->first;
		Vertex* v1 = it->second.first;
		Vertex* v2 = it->second.second;
		int v1id = v_idx(v1);
		int v2id = v_idx(v2);
		zid2vertexIdPair[zid] = std::pair<int,int>(v1id,v2id);
		n_extra_constrants += 2;
		std::cout << "zip pair " << zid << "\t" << v1id << "\t" << v2id << "\t" << v1 << "\t" << v2 << "\n";
	}

	for (int k = 0; k < (int)m_fix_vertices.size(); k++) {
		v_idx(m_fix_vertices[k]) = k;
		Vertex *v = m_fix_vertices[k];
		std::string tmp;
		double uv0, uv1;
		std::stringstream(v->string()) >> tmp >> uv0 >> uv1;
		v_uv(v) = Point2(uv0, uv1);
	}

	int fn = m_mesh->numFaces();
	int	vfn = m_fix_vertices.size();
	int	vn = m_mesh->numVertices() - vfn;
	
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList1;
	std::vector<T> tripletList2;
	tripletList1.reserve(fn+n_extra_constrants);
	tripletList2.reserve(fn+n_extra_constrants);
	VectorXd b(vfn*2 + 1);
	b[ vfn * 2]=1.0;
	int fid = 0;
	for (MeshFaceIterator fiter(m_mesh); !fiter.end(); ++fiter, ++fid) {
		Face *f = *fiter;
		HalfEdge *he = f->halfedge();

		for (int j = 0; j < 3; j++) {
			Point s = c_s(he);

			Vertex * v = he->he_next()->target();
			int vid = v_idx(v);

			if (v->string().substr(0,3) != "fix") {
				tripletList1.push_back(T(fid,vid,s[0]));
				tripletList1.push_back(T(fn + fid, vn + vid, s[0]));
				tripletList1.push_back(T(fid, vn + vid, -s[1]));
				tripletList1.push_back(T(fn + fid, vid, s[1]));

				assert( fid < (2*fn + n_extra_constrants) );
				assert( (fn + fid) < (2*fn + n_extra_constrants) );
				assert( (vid) < (2*vn) );
				assert( (vn+vid) < (2*vn) );

			}
			else {
				Point2 uv = v_uv(v);

				tripletList2.push_back(T(fid, vid, s[0]));
				tripletList2.push_back(T(fn + fid, vfn + vid, s[0]));
				tripletList2.push_back(T(fid, vfn + vid, -s[1]));
				tripletList2.push_back(T(fn + fid, vid, s[1]));

				b[vid] = uv[0];
				b[vfn + vid] = uv[1];
			}
			he = he->he_next();
		}
	}


	for(auto it=zid2vertexIdPair.begin();it!=zid2vertexIdPair.end();it++) {
		int zid = it->first;
		std::cout << "zid: " << zid << "\n";
		int v1id = it->second.first;
		int v2id = it->second.second;

		// u coordinate constrained to have difference of 1:
		tripletList1.push_back(T(2*fn + zid,v1id,-1.0));
		tripletList1.push_back(T(2*fn + zid,v2id, 1.0));
		tripletList2.push_back(T(2*fn + zid, (vfn * 2), -1.0));

			if (!((2*fn + n_zip/2 +  zid) < (2*fn + n_extra_constrants+1))) {
		std::cout << "zid: " << zid << "\n";

				std::cout << "index out of range:" 
				<< "2*fn+n_zip/2+zid=" << 2 <<"*" << fn <<"+"<<n_zip/2 << "+"<<zid << " = "
				<< 2*fn + n_zip +  zid << ">=" << 2*fn + n_extra_constrants << "=" << 2 << "*"<<fn<<"+"<<  n_extra_constrants << "\n";
			}
				assert( 2*fn + n_zip/2 +  zid < (2*fn + n_extra_constrants+1) );
				assert( (v1id) < (2*vn) );
				assert( (v2id) < (2*vn) );


		// v coordinate constrained to have difference of 0:
		tripletList1.push_back(T(2*fn + n_zip/2 + zid, vn + v1id, 1.0));
		tripletList1.push_back(T(2*fn + n_zip/2 + zid, vn + v2id,-1.0));
		//tripletList2.push_back(T(2*fn + n_zip+ zid, (vfn * 2), 0.0));  <- sparse matrix;  it's implicity zero anyway

	}


	SpMat A(2*fn + n_extra_constrants+1, 2*vn);
	SpMat B(2*fn + n_extra_constrants+1, 2*vfn + 1);

	A.setFromTriplets(tripletList1.begin(), tripletList1.end());
	B.setFromTriplets(tripletList2.begin(), tripletList2.end());

	VectorXd r, x;
	r = B * b;
	r = r * -1;

    // Solve the linear system Ax = r
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> lscg;
	lscg.compute(A);
	x = lscg.solve(r);

	for (int i = 0; i < vn; i++) {
		Vertex * v = vertices[i];
		v_uv(v) = Point2(x[i], x[i + vn]);
	}
	for (MeshVertexIterator viter(m_mesh); !viter.end(); viter++) {
		Vertex * v = *viter;
		Point2 p = v_uv(v);
		Point p3(p[0], p[1], 0);
		v->point() = p3;
	}
}
