#include <iostream>
#include <math.h>
#include <utility>
#include <vector>

#include "Eigen/Eigen"
#include "Eigen/Geometry"

using namespace Eigen;
using namespace std;

MatrixXd Euler2A(double alfa, double beta, double gama){
	cout << "Uglovi: " << alfa*180/M_PI << " " << beta*180/M_PI << " " << gama*180/M_PI << endl << endl;

	double alfasin = sin(alfa);
	double alfacos = cos(alfa);

	double betasin = sin(beta);
	double betacos = cos(beta);

	double gamasin = sin(gama);
	double gamacos = cos(gama);

	Matrix3d RXalfa;
	Matrix3d RYbeta;
	Matrix3d RZgamma;

	RXalfa << 1,       0,        0, 
	          0, alfacos, -alfasin,
	          0, alfasin,  alfacos;

	RYbeta << betacos, 0, betasin,
	                0, 1,       0,
	         -betasin, 0, betacos;

	RZgamma << gammacos, -gammasin, 0,
	           gammasin,  gammacos, 0,
	                  0,         0, 1;

	MatrixXd E2A;
	E2A = RZgamma * RYbeta * RXalfa;

	cout << "Euler2A: " << endl;
	cout << E2A << endl << endl;

	return E2A;
}

pair<Vector3d, double> AxisAngle(Matrix3d &A){

	//funkcija 2. : AxisAngle[A] - vraca jedinicni vektor p = (px, py, pz) i ugao alfa(0, pi) tako da A = Rp(alfa)

    Vector3d p;
    double ugao;

	//vrsimo provere da li je A*A.transpose() = E i da li je A.det = 1
	Matrix3d E;

	E << 1, 0, 0,
	     0, 1, 0,
	     0, 0, 1;

	Matrix3d ATA;

	ATA = A.transpose() * A;

	Matrix3d RATA;
	RATA << abs(round(ATA(0, 0))), abs(round(ATA(0, 1))), abs(round(ATA(0, 2))),
	            abs(round(ATA(1, 0))), abs(round(ATA(1, 1))), abs(round(ATA(1, 2))),
	            abs(round(ATA(2, 0))), abs(round(ATA(2, 1))), abs(round(ATA(2, 2)));

	if(RATA != E){
		cout << "Neodgovarajuca matrica" << endl;
		exit(EXIT_FAILURE);
	}

	if(round(A.determinant()) != 1){
		cout << "Neodgovarajuca matrica" << endl;
		exit(EXIT_FAILURE);
	}

	//vektorski proizvod -> p 
	Vector3d v1;

	v1 << A(0, 0)-1, A(0, 1), A(0, 2);

	Vector3d v2;

	//zavisni vektori?
	if( (A(0, 0)-1)*1.0/A(1, 0) == A(0, 1)*1.0/(A(1, 1)-1) && A(0, 1)*1.0/(A(1, 1)-1) == A(0, 2)*1.0/A(1, 2)){
		v2 << A(2, 0), A(2, 1), A(2, 2-1);
    }
	else
		v2 << A(1, 0), A(1, 1)-1, A(1, 2);

	//normalizovanje
	v1 = v1.normalized();
	v2 = v2.normalized();

    //vektorski proizvod
	p = v1.cross(v2);
	p = p.normalized();

	Vector3d v1rotiran;

	v1rotiran = A * v1;

	ugao = std::acos(v1.dot(v1rotiran));

	//mesoviti proizvod
	Matrix3d mix;
	mix << p, v1, v1rotiran;
	if(mix.determinant() < 0){
		p = -1 * p;
	}

	cout << "AxisAngle[A]: " << endl << endl 
              << "Ugao: " << ugao << endl << endl 
			  << "Vektor: " << endl << p << endl << endl;

	return std::make_pair(p, ugao);
}

Matrix3d Rodrigez(Vector3d &p, double ugao){
	
	Matrix3d Px;

	Px <<  0, -p(2, 0),  p(1, 0),
	       p(2, 0), 0, -p(0, 0),
	      -p(1, 0),  p(0, 0), 0;

	Matrix3d E;

	E << 1, 0, 0,
	     0, 1, 0,
	     0, 0, 1;

	Matrix3d MatricaR;

	MatricaR = p*p.transpose() + cos(ugao)*(E - p*p.transpose()) + sin(ugao)*Px;

	cout << "Rodrigez[p, fi]: " << endl << MatricaR << endl << endl;

	return MatricaR;
}

void A2Euler(Matrix3d &MatricaR){
	cout << "A2Euler[A]:" << endl;

	double psi, teta, fi; 
		if (MatricaR(2, 0) < 1){
			if (MatricaR(2, 0) > -1){ 
				psi = atan2(MatricaR(1, 0), MatricaR(0, 0));
				teta = asin(-MatricaR(2, 0));
				fi = atan2(MatricaR(2, 1), MatricaR(2, 2));
			} else { 
				psi = atan2(-MatricaR(0, 1), MatricaR(1, 1));
				teta = M_PI/2;
				fi = 0;
			}
			} else { 
				psi = atan2(-MatricaR(0, 1), MatricaR(1, 1));
				teta = -M_PI/2;
			    fi = 0;
			}

	cout << "Ugao fi: " << fi*180/M_PI << endl
		 << "Ugao teta: " << teta*180/M_PI << endl 
	     << "Ugao psi: " << psi*180/M_PI << endl 
		 << endl << endl;

	return;
}

int main(){

    double alfa = 30*M_PI/180;
	double beta = 60*M_PI/180;
	double gama = 90*M_PI/180;

	Matrix3d E2A;

	E2A = Euler2A(alfa, beta, gama);

    Vector3d p;
    double ugao;

	auto axisAngle = AxisAngle(E2A);

	p = axisAngle.first;
	ugao = axisAngle.second;

	Matrix3d MatricaR = Rodrigez(p, ugao);


	A2Euler(MatricaR);

	cout << "AxisAngle2Q[p, fi]: " << endl;

	MatrixXd q;

	q.resize(1, 4);

	q << p(0, 0)*sin(ugao/2), p(1, 0)*sin(ugao/2), p(2, 0)*sin(ugao/2), cos(ugao/2);

	cout << "Kvaternion: " << endl << q << endl << endl;
	
	cout << "Q2AxisAngle[q]: " << endl;

	double ugaoFiPom = acos(q(0, 3));
	Vector3d vektorQ2;

	vektorQ2 << q(0, 0)/sin(ugaoFiPom/2), q(0, 1)/sin(ugaoFiPom/2), q(0, 2)/sin(ugaoFiPom/2);

	cout << "Vektor: " << endl << vektorQ2.normalized() << endl << endl << "Ugao: " << ugaoFiPom*2 << endl << endl;

	return 0;
}
