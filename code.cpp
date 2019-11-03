#include <iostream>
#include "Eigen/Dense"
#include "Eigen/SVD"
#include "CImg-2.7.5/CImg.h"
#include <vector>

using namespace std;
using namespace cimg_library;
using namespace Eigen;

Matrix3d poslednja_kolona;	//DLT
Matrix3d rezultujuca_matrica; //DLT norm

//DLT i DLT_norm implementirani kao fje, naivni direktno u main-u

Matrix3d &DLT(vector<vector<double>> &originalne_tacke, vector<vector<double>> &tacke_slike, int broj_tacaka)
{
	//za svaku tacku pravimo dve jednacine i cuvamo u zasebnim matricama, 
    //koje spajamo u globalnu matricu nepoznatih
	MatrixXd globalna_matrica;
	globalna_matrica.resize(4, 9);

	MatrixXd pom1;
	pom1.resize(2, 9);

	MatrixXd pom2;
	pom2.resize(2, 9);

	// postoji bar 1 nepoznata, 2 jednacine
	// 0, 0, 0, -x1*x3', -x2*x3', -x3*x3', x1*x2', x2*x2', x3*x3'
	// x1*x3', x2*x3', x3*x3', 0, 0, 0, -x1*x1', -x2*x1', -x3*x1'
	// x1 = originalne_tacke[0][0], x1 = originalne_tacke[0][1], x3 = originalne_tacke[0][2]
	// x1' = tacke_slike[0][0], x1 = tacke_slike[0][1], x3 = tacke_slike[0][2]
	pom1 << 0, 0, 0,
		-originalne_tacke[0][0] * tacke_slike[0][2], -originalne_tacke[0][1] * tacke_slike[0][2], -originalne_tacke[0][2] * tacke_slike[0][2],
		originalne_tacke[0][0] * tacke_slike[0][1], originalne_tacke[0][1] * tacke_slike[0][1], originalne_tacke[0][2] * tacke_slike[0][1],
		originalne_tacke[0][0] * tacke_slike[0][2], originalne_tacke[0][1] * tacke_slike[0][2], originalne_tacke[0][2] * tacke_slike[0][2],
		0, 0, 0,
		-originalne_tacke[0][0] * tacke_slike[0][0], -originalne_tacke[0][1] * tacke_slike[0][0], -originalne_tacke[0][2] * tacke_slike[0][0];

	int i = 1;
	//iteriramo dok ima nepoznatih
	while (i < broj_tacaka)
	{
		//matrica ista kao pom1, menja se samo red matrice tacaka
		pom2 << 0, 0, 0,
			-originalne_tacke[i][0] * tacke_slike[i][2], -originalne_tacke[i][1] * tacke_slike[i][2], -originalne_tacke[i][2] * tacke_slike[i][2],
			originalne_tacke[i][0] * tacke_slike[i][1], originalne_tacke[i][1] * tacke_slike[i][1], originalne_tacke[i][2] * tacke_slike[i][1],
			originalne_tacke[i][0] * tacke_slike[i][2], originalne_tacke[i][1] * tacke_slike[i][2], originalne_tacke[i][2] * tacke_slike[i][2],
			0, 0, 0,
			-originalne_tacke[i][0] * tacke_slike[i][0], -originalne_tacke[i][1] * tacke_slike[i][0], -originalne_tacke[i][2] * tacke_slike[i][0];

		globalna_matrica << pom1, pom2; /*unesemo pomocne u globalnu*/
		pom1.resize((i + 1) * 2, 9);
		pom1 = globalna_matrica; /*pom1 sada cuva sve do sada azurirano*/
		globalna_matrica.resize((i + 2) * 2, 9); /*prosirujemo globalnu za nove jednacine*/
		i++;
	}

	//Ispis globalne matrice
	cout << "Globalna matrica (2*broj_tacaka, 9):" << endl;
	cout << pom1 << endl<< endl;

	//primenjujemo SVD
	//SVD dekompozicija pom1 = UDV.transponovano()
	BDCSVD<MatrixXd> svd(pom1, ComputeFullU | ComputeFullV);

	cout << "Matrica V dobijena SVD dekompozicijom: " << endl;
	cout << svd.matrixV() << endl << endl;

	//Matrix3d poslednja_kolona;
	MatrixXd V;
	V.resize(9, 9);
	V << svd.matrixV();

	//matrica transformacije je poslednja kolona matrice V
	poslednja_kolona << V(0, 8), V(1, 8), V(2, 8),
		                V(3, 8), V(4, 8), V(5, 8),
		                V(6, 8), V(7, 8), V(8, 8);

	return poslednja_kolona;
}

Matrix3d &DLT_normalizovani(vector<vector<double>> &originalne_tacke, vector<vector<double>> &tacke_slike, int broj_tacaka)
{
	//teziste originalnih tacaka i tacaka slike
	int zbir_x, zbir_y, zbir_x_slika, zbir_y_slika;

	//original - zbir po koordinatama
	for (int i = 0; i < broj_tacaka; i++)
	{
		zbir_x += originalne_tacke[i][0];
		zbir_y += originalne_tacke[i][1];
	}

	//slika - zbir po koordinatama
	for (int j = 0; j < broj_tacaka; j++)
	{
		zbir_x_slika += tacke_slike[j][0];
		zbir_y_slika += tacke_slike[j][1];
	}

	//original - centar
	vector<double> centar;
	double centar_x = zbir_x * 1.0 / broj_tacaka;
	double centar_y = zbir_y * 1.0 / broj_tacaka;
	centar.push_back(centar_x);
	centar.push_back(centar_y);
	cout << "Koordinate centra originalne slike: (" << centar[0] << ", " << centar[1] << ")" << endl << endl;

	//slika - centar
	vector<double> centar_s;
	double centar_x_s = zbir_x_slika * 1.0 / broj_tacaka;
	double centar_y_s = zbir_y_slika * 1.0 / broj_tacaka;
	centar_s.push_back(centar_x);
	centar_s.push_back(centar_y);
	cout << "Koordinate centra slike: (" << centar_s[0] << ", " << centar_s[1] << ")" << endl << endl;

	//prosecno rastojanje tacaka od centra
	double zbir_rastojanja = 0, zbir_rastojanja_s = 0;

	//original - zbir rastojanja
	for (int i = 0; i < broj_tacaka; i++)
	{
		zbir_rastojanja += sqrt((centar[0] - originalne_tacke[i][0]) * (centar[0] - originalne_tacke[i][0]) +
								(centar[1] - originalne_tacke[i][1]) * (centar[1] - originalne_tacke[i][1]));
	}
	
    //slika - zbir rastojanja
	for (int i = 0; i < broj_tacaka; i++)
	{
		zbir_rastojanja_s += sqrt((centar_s[0] - tacke_slike[i][0]) * (centar_s[0] - tacke_slike[i][0]) +
								  (centar_s[1] - tacke_slike[i][1]) * (centar_s[1] - tacke_slike[i][1]));
	}
	
    //prosek za original i sliku
	double prosecna_udaljenost = zbir_rastojanja / broj_tacaka;
	double prosecna_udaljenost_s = zbir_rastojanja_s / broj_tacaka;

	cout << "Prosecna udaljenost tacaka od centra (originalne tacke): " << prosecna_udaljenost << endl << endl;
	cout << "Prosecna udaljenost tacaka od centra (tacke slike): " << prosecna_udaljenost_s << endl << endl;

	//Norm i Norm_s su matrice koje normalizuju tacke originala i slike
	Matrix3d Norm, Norm_s;

	Norm << (sqrt(2) / prosecna_udaljenost), 0, -centar[0],
		    0, (sqrt(2) / prosecna_udaljenost), -centar[1],
		    0, 0, 1;

	Norm_s << (sqrt(2) / prosecna_udaljenost_s), 0, -centar_s[0],
		    0, (sqrt(2) / prosecna_udaljenost_s), -centar_s[1],
		    0, 0, 1;

	cout << "Matrica za normalizaciju originalnih tacaka: " << endl;
	cout << Norm << endl << endl;

	cout << "Matrica za normalizaciju tacaka slike: " << endl;
	cout << Norm_s << endl << endl;

	//normiranje

	//original
	MatrixXd pom_tacke;
	pom_tacke.resize(3, 1);
	MatrixXd pom_tacke_Norm;
	
	vector<vector<double>> norm_originalne;
	vector<double> pom;

	for (int i = 0; i < broj_tacaka; i++)
	{
		pom_tacke << originalne_tacke[i][0], originalne_tacke[i][1], originalne_tacke[i][2]; //uzmi tacku
		pom_tacke_Norm = Norm * pom_tacke;	//primeni Norm za originalne
		//mnozimo matricu 3x3 sa vektorom 3x1 -> dobijamo vektor 3x1
		double x, y, z;
		x = pom_tacke_Norm(0, 0);
		y = pom_tacke_Norm(1, 0);
		z = pom_tacke_Norm(2, 0);
		pom.push_back(x);
		pom.push_back(y);
		pom.push_back(z);
		norm_originalne.push_back(pom); //popunjavamo matricu normiranih originalnih tacaka
		pom.clear();
	}

	cout << "Normalizovane tacke originala" << endl;
	for (vector<double> beka : norm_originalne)
	{
		for (double b : beka)
		{
			cout << b << " ";
		}
		cout << endl;
	}
	cout << endl;
	
	//slika
	MatrixXd pom_tacke_s;
	pom_tacke_s.resize(3, 1);
	MatrixXd pom_tacke_Norm_s;

	vector<vector<double> > norm_tacke_slike;
	vector<double> pom_s;

	for (int i = 0; i < broj_tacaka; i++)
	{
		pom_tacke_s << tacke_slike[i][0], tacke_slike[i][1], tacke_slike[i][2]; //uzmi tacku
		pom_tacke_Norm_s = Norm_s * pom_tacke_s; //primeni Norm za tacke slike
		//kao i za tacke slike: matrica 3x3 * vektor 3x1 -> vektor 3x1
		double xs, ys, zs;
		xs = pom_tacke_Norm_s(0, 0);
		ys = pom_tacke_Norm_s(1, 0);
		zs = pom_tacke_Norm_s(2, 0);
		pom_s.push_back(xs);
		pom_s.push_back(ys);
		pom_s.push_back(zs);
		norm_tacke_slike.push_back(pom_s);
		pom_s.clear();
	}

	cout << "Normalizovane tacke originala" << endl;
	for (vector<double> beki : norm_tacke_slike)
	{
		for (double b : beki)
		{
			cout << b << " ";
		}
		cout << endl;
	}
	cout << endl;

	//Algoritam DLT nad normiranim tackama
	Matrix3d DLT_norm_matrica = DLT(norm_originalne, norm_tacke_slike, broj_tacaka);
	cout << "Matrica transformacije normiranih tacaka (DLT): " << endl;
	cout << DLT_norm_matrica << endl << endl;

	//Algoritam DLT norm
	//Matrix3d rezultujuca_matrica; --globalno deklarisana
	rezultujuca_matrica = Norm_s.inverse() * DLT_norm_matrica * Norm;

	cout << "Matrica transormacije normalizovanim DLT algoritmom: " << endl;
	cout << rezultujuca_matrica << endl << endl;

	return rezultujuca_matrica;
}

int main()
{
	cout << "Da biste pokrenuli algoritam, potrebno je da odaberete nacin unosa tacaka na originalnoj slici." << endl;
	cout << "Opcija 1: samostalni unos koordinata tacaka " << endl;
	cout << "Opcija 2: unos tacaka klikom " << endl;
	int opcija;
	cin >> opcija;

	if (opcija != 1 && opcija != 2)
	{
		cout << "Uneli ste pogresnu opciju. Moguce: 1 ili 2." << endl;
		return -1;
	}

	vector<vector<double>> originalne_tacke;
	int broj_tacaka = 4; /*default*/

	CImg<unsigned char> originalnaSlika("trijumfalna_kapija.bmp"),
		originalnaSlikaa("trijumfalna_kapija.bmp"),
		slikaSLike(500, 400, 1, 3, 0);

	int width;
	int height;

	width = 500;
	height = 400;
	originalnaSlika.resize(width, height);
	originalnaSlikaa.resize(width, height);


	vector<double> pomocni;
	if (opcija == 2)
	{
		CImgDisplay org_prikaz(originalnaSlika, "Odaberi 4 tacke klikom na sliku, a zatim zatvori prozor");
		while (!org_prikaz.is_closed())
		{

			org_prikaz.wait();
			if (org_prikaz.button() && org_prikaz.mouse_y() >= 0 && org_prikaz.mouse_x() >= 0)
			{
				const int x = org_prikaz.mouse_x();
				const int y = org_prikaz.mouse_y();

				pomocni.push_back(x);
				pomocni.push_back(y);
				pomocni.push_back(1);

				originalne_tacke.push_back(pomocni);
				pomocni.clear();
			}
		}

        cout << "Odabrali ste sledece tacke na slici:" << endl;
		for (vector<double> beks : originalne_tacke)
		{
			for (double b : beks)
			{
				cout << b << " ";
			}
			cout << endl;
		}
		broj_tacaka = originalne_tacke.size(); //moze da bude >= 4
	}

	//u slucaju opcije 1:
	if (opcija == 1)
	{
        cout << "Koliko tacaka zelite da unesete?" << endl;
		cin >> broj_tacaka;
		cout << "Unesite koordinate originalnih tacaka:" << endl;
		vector<double> pom1;
		double x;
		for (int i = 0; i < broj_tacaka; i++)
		{
			pom1.clear();
			for (int j = 0; j < 3; j++)
			{
				cin >> x;
				pom1.push_back(x);
			}
			originalne_tacke.push_back(pom1);
		}
	}
    cout << endl;
	cout << "Unesite koordinate tacaka slike:" << endl;

	vector<vector<double>> tacke_slike;
	vector<double> pom2;
	double y;
	for (int i = 0; i < broj_tacaka; i++)
	{
		pom2.clear();
		for (int j = 0; j < 3; j++)
		{
			cin >> y;
			pom2.push_back(y);
		}
		tacke_slike.push_back(pom2);
	}
    cout << endl;
	//kraj ucitavanja tacaka

	/*NAIVNI ALGORITAM*/

	//Kramerovo pravilo za P
	Matrix3d delta;
	delta << originalne_tacke[0][0], originalne_tacke[1][0], originalne_tacke[2][0],
		originalne_tacke[0][1], originalne_tacke[1][1], originalne_tacke[2][1],
		originalne_tacke[0][2], originalne_tacke[1][2], originalne_tacke[2][2];

	Matrix3d delta1;
	delta1 << originalne_tacke[3][0], originalne_tacke[1][0], originalne_tacke[2][0],
		originalne_tacke[3][1], originalne_tacke[1][1], originalne_tacke[2][1],
		originalne_tacke[3][2], originalne_tacke[1][2], originalne_tacke[2][2];

	Matrix3d delta2;
	delta2 << originalne_tacke[0][0], originalne_tacke[3][0], originalne_tacke[2][0],
		originalne_tacke[0][1], originalne_tacke[3][1], originalne_tacke[2][1],
		originalne_tacke[0][2], originalne_tacke[3][2], originalne_tacke[2][2];

	Matrix3d delta3;
	delta3 << originalne_tacke[0][0], originalne_tacke[1][0], originalne_tacke[3][0],
		originalne_tacke[0][1], originalne_tacke[1][1], originalne_tacke[3][1],
		originalne_tacke[0][2], originalne_tacke[1][2], originalne_tacke[3][2];

	double delta_det = delta.determinant();
	double delta1_det = delta1.determinant();
	double delta2_det = delta2.determinant();
	double delta3_det = delta3.determinant();

	double lambda1 = delta1_det * 1.0 / delta_det;
	double lambda2 = delta2_det * 1.0 / delta_det;
	double lambda3 = delta3_det * 1.0 / delta_det;

	//Provera lambdi
	cout << "Lambde za matricu P: (" << lambda1 << ", " << lambda2 << ", " << lambda3 << ")" << endl << endl;

    //(preslikavanje A0,B0,C0,D0  ---> A,B,C,D)*/
	Matrix3d matrica_preslikavanja_P;
	
	matrica_preslikavanja_P << originalne_tacke[0][0] * lambda1, originalne_tacke[1][0] * lambda2, originalne_tacke[2][0] * lambda3,
		originalne_tacke[0][1] * lambda1, originalne_tacke[1][1] * lambda2, originalne_tacke[2][1] * lambda3,
		originalne_tacke[0][2] * lambda1, originalne_tacke[1][2] * lambda2, originalne_tacke[2][2] * lambda3;
	
    cout << "Matrica preslikavanja P: " << endl;
	cout << matrica_preslikavanja_P << endl << endl;

	//Provera matrice P
	MatrixXd D0;
	D0.resize(3, 1);
	D0 << 1, 1, 1;
	MatrixXd D;

	D = matrica_preslikavanja_P * D0; //primena na kanonsku tacku D0(1,1,1)

	cout << "Provera da li matrica P preslikava tacku D0 u D: " << endl;
	cout << D << endl << endl; //dobijamo bas D

	/*******************************************************************/

	//Kramerovo pravilo za G
	Matrix3d deltag;
	deltag << tacke_slike[0][0], tacke_slike[1][0], tacke_slike[2][0],
		tacke_slike[0][1], tacke_slike[1][1], tacke_slike[2][1],
		tacke_slike[0][2], tacke_slike[1][2], tacke_slike[2][2];

	Matrix3d delta1g;
	delta1g << tacke_slike[3][0], tacke_slike[1][0], tacke_slike[2][0],
		tacke_slike[3][1], tacke_slike[1][1], tacke_slike[2][1],
		tacke_slike[3][2], tacke_slike[1][2], tacke_slike[2][2];

	Matrix3d delta2g;
	delta2g << tacke_slike[0][0], tacke_slike[3][0], tacke_slike[2][0],
		tacke_slike[0][1], tacke_slike[3][1], tacke_slike[2][1],
		tacke_slike[0][2], tacke_slike[3][2], tacke_slike[2][2];

	Matrix3d delta3g;
	delta3g << tacke_slike[0][0], tacke_slike[1][0], tacke_slike[3][0],
		tacke_slike[0][1], tacke_slike[1][1], tacke_slike[3][1],
		tacke_slike[0][2], tacke_slike[1][2], tacke_slike[3][2];

	double delta_detg = deltag.determinant();
	double delta1_detg = delta1g.determinant();
	double delta2_detg = delta2g.determinant();
	double delta3_detg = delta3g.determinant();

	double lambda1g = delta1_detg * 1.0 / delta_detg;
	double lambda2g = delta2_detg * 1.0 / delta_detg;
	double lambda3g = delta3_detg * 1.0 / delta_detg;

	//Provera lambdi
	cout << "Lambde za matricu G: (" << lambda1g << ", " << lambda2g << ", " << lambda3g<< ")" << endl << endl;

    //G je matrica preslikavanja A0,B0,C0,D0  ---> A',B',C',D'
	Matrix3d matrica_preslikavanja_G;

	matrica_preslikavanja_G << tacke_slike[0][0] * lambda1g, tacke_slike[1][0] * lambda2g, tacke_slike[2][0] * lambda3g,
		tacke_slike[0][1] * lambda1g, tacke_slike[1][1] * lambda2g, tacke_slike[2][1] * lambda3g,
		tacke_slike[0][2] * lambda1g, tacke_slike[1][2] * lambda2g, tacke_slike[2][2] * lambda3g;

    cout << "Matrica preslikavanja G: " << endl;
	cout << matrica_preslikavanja_G << endl << endl;

	//Provera matrice G
	MatrixXd Dg;
	Dg = matrica_preslikavanja_G * D0;

	cout << "Provera da li matrica G preslikava tacku D0 u D': " << endl;
	cout << Dg << endl << endl;

	//Primena naivnog algoritma
	//Naivni algoritam : Kompozicija preslikavanja P.inverz i G

	Matrix3d rezultujuce_preslikavanje;
	rezultujuce_preslikavanje = matrica_preslikavanja_G * matrica_preslikavanja_P.inverse();

	cout << "Rezultujuca matrica preslikavanja naivnog algoritma: " << endl;

	cout << rezultujuce_preslikavanje << endl << endl;

	//iscrtavanje prozora sa izmenjenom slikom
	int visina = 400, sirina = 500; //za sliku

	if (opcija == 2)
	{
		CImgDisplay izmenjeni_prikaz(slikaSLike, "Izmenjena slika"),
			org_prikazz(originalnaSlikaa, "Originalna slika");

		MatrixXd pom_org;
		pom_org.resize(3, 1);

		MatrixXd pom_slika;
		pom_slika.resize(3, 1);

		int i_org, j_org;

		//Pinv - invertovana matrica dobijena naivnim algoritmom
		Matrix3d P_inv;
		P_inv = rezultujuce_preslikavanje.inverse();

        //moramo biti u okvirima prozora za izmenjenu sliku
		for (int i = 0; i < sirina; i++)
			for (int j = 0; j < visina; j++)
			{
				pom_slika << i, j, 1;
				pom_org = P_inv * pom_slika; //primenimo inverz na piksel slike

				//pretvaramo homogene koordinate piksela u afine
				if (pom_org(2, 0) != 0.0)
				{

					i_org = ceil(pom_org(0, 0) / pom_org(2, 0));
					j_org = ceil(pom_org(1, 0) / pom_org(2, 0));
				}
				else
				{
					continue;
				}
				if (j_org > height || j_org < 0 || i_org > width || i_org < 0) //ako smo iskocili nastavi
					continue;

				slikaSLike(i, j, 0, 0) = (int)originalnaSlika(i_org, j_org, 0, 0); // 0 0 - crvena
				slikaSLike(i, j, 0, 1) = (int)originalnaSlika(i_org, j_org, 0, 1); // 0 1 - zelena
				slikaSLike(i, j, 0, 2) = (int)originalnaSlika(i_org, j_org, 0, 2); // 0 2 - plava
			}
		
		slikaSLike.display(izmenjeni_prikaz);
		originalnaSlika.display(org_prikazz);
		//dok ne zatvorimo - prozori su aktivni
		while (!izmenjeni_prikaz.is_closed() && !org_prikazz.is_closed())
		{
			izmenjeni_prikaz.wait();
			org_prikazz.wait();
		}
	}

	
	/*DLT algoritam*/

	Matrix3d matrica_dlt = DLT(originalne_tacke, tacke_slike, broj_tacaka);
	cout << "Matrica transformacije dobijena DLT algoritmom: " << endl;
	cout << matrica_dlt << endl << endl;

	//provera ispravnosti DLT-a
	Matrix3d DLT_provera;
	Matrix3d P;
	P << rezultujuce_preslikavanje; //u P smestamo rezultat naivnog alg
	// matricu iz DLT algoritma mnozimo sa P(0, 0) i delimo sa matrica_dlt(0,0)
	// P i DLT matrice su proporcionalne, pa bi ovim skaliranjem trebalo da dobijemo istu matricu
	DLT_provera << matrica_dlt(0, 0) * P(0, 0) / matrica_dlt(0, 0),
		matrica_dlt(0, 1) * P(0, 0) / matrica_dlt(0, 0),
		matrica_dlt(0, 2) * P(0, 0) / matrica_dlt(0, 0),
		matrica_dlt(1, 0) * P(0, 0) / matrica_dlt(0, 0),
		matrica_dlt(1, 1) * P(0, 0) / matrica_dlt(0, 0),
		matrica_dlt(1, 2) * P(0, 0) / matrica_dlt(0, 0),
		matrica_dlt(2, 0) * P(0, 0) / matrica_dlt(0, 0),
		matrica_dlt(2, 1) * P(0, 0) / matrica_dlt(0, 0),
		matrica_dlt(2, 2) * P(0, 0) / matrica_dlt(0, 0);

	cout << "Provera proporcionalnosti DLT i naivnog:" << endl;
	cout << DLT_provera << endl << endl;

	/*Modifikovani DLT algoritam (sa normalizacijom)*/

	//moramo da prebacimo iz homogenih u homogene sa trecom koordinatom 1,
    // za slucaj da treca koordinata nije 1 nego neki drugi broj (kao vid afinih jer su nam bitni x i y)
	vector<vector<double>> originalne_tacke_afino, tacke_slike_afino;
	vector<double> pom_o, pom_s;

	for (int i = 0; i < broj_tacaka; i++)
	{
			double x_afino = originalne_tacke[i][0] / originalne_tacke[i][2];
			double y_afino = originalne_tacke[i][1] / originalne_tacke[i][2];
			pom_o.push_back(x_afino);
			pom_o.push_back(y_afino);
			pom_o.push_back(1);
			
			originalne_tacke_afino.push_back(pom_o);
			pom_o.clear();
		
			double x_afino_s = tacke_slike[i][0] / tacke_slike[i][2];
			double y_afino_s = tacke_slike[i][1] / tacke_slike[i][2];
			pom_s.push_back(x_afino_s);
			pom_s.push_back(y_afino_s);
			pom_s.push_back(1);

			tacke_slike_afino.push_back(pom_s);
			pom_s.clear();
	
	}
	cout << "Originalne tacke sa 1 kao trecom koordinatom:" << endl;
	for(vector<double> provera : originalne_tacke_afino){
		for(double pro : provera){
			cout << pro << " ";
		}
		cout << endl;
	}
	cout << endl;

	cout << "Tacke slike sa 1 kao trecom koordinatom:" << endl;
	for(vector<double> proveras : tacke_slike_afino){
		for(double pros : proveras){
			cout << pros << " ";
		}
		cout << endl;
	}
	cout << endl;

    //primena DLT_norm
	Matrix3d dlt_norm = DLT_normalizovani(originalne_tacke_afino, tacke_slike_afino, broj_tacaka);
	Matrix3d dlt_norm_skaliran;

	dlt_norm_skaliran << dlt_norm(0, 0) * P(0, 0) / dlt_norm(0, 0),
		dlt_norm(0, 1) * P(0, 0) / dlt_norm(0, 0),
		dlt_norm(0, 2) * P(0, 0) / dlt_norm(0, 0),
		dlt_norm(1, 0) * P(0, 0) / dlt_norm(0, 0),
		dlt_norm(1, 1) * P(0, 0) / dlt_norm(0, 0),
		dlt_norm(1, 2) * P(0, 0) / dlt_norm(0, 0),
		dlt_norm(2, 0) * P(0, 0) / dlt_norm(0, 0),
		dlt_norm(2, 1) * P(0, 0) / dlt_norm(0, 0),
		dlt_norm(2, 2) * P(0, 0) / dlt_norm(0, 0);

	cout << "Provera proporcionalnosti DLT_norm i naivnog algoritma: " << endl;
	cout << dlt_norm_skaliran << endl << endl;

	//Provera da li DLT_normalizovani slika D u Dp

	MatrixXd provera_dlt_norm;
	provera_dlt_norm = dlt_norm_skaliran * D;

	cout << "Provera da li dlt_norm_skaliran slika D u Dp" << endl;
	cout << provera_dlt_norm << endl << endl;
	return 0;
}