#include <iostream>
#include "math.h"
#include <fstream>
/*#include "OH_fct.h"
*/
using namespace std;

void choix_initial (double &V0, double &largeur)
{
  cout << "Choisissez la profondeur du potentiel (eV) : ";
  cin >> V0;
}









// Fonctions de la methode RK4
// Calcul des coefficients K et L et des coefficients permettant le calcul de psi et phi
double K (double phi)
{
  double k = phi;
  return k;
}

double L (double psi, double puls, double potentiel, long double energie)
{
  double l;
  if ((potentiel - energie) < 0)
    {
      l = -pow (puls, 2) * psi;
    }
  if ((potentiel - energie) > 0)
    {
      l = +pow (puls, 2) * psi;
    }
  return l;
}

void coeff (double psi, double phi, double &a, double &b, double &c, double &d,
       double &m, double &n, double &o, double &p, double pas, double puls,
       double potentiel, long double energie)
{
  m = K (phi);			// k1
  a = L (psi, puls, potentiel, energie);	// l1
  n = K (phi + (pas / 2) * a);	// k2
  b = L (psi + (pas / 2) * m, puls, potentiel, energie);	// l2
  o = K (phi + (pas / 2) * b);	// k3
  c = L (psi + (pas / 2) * n, puls, potentiel, energie);	// l3
  p = K (phi + pas * c);	// k4
  d = L (psi + pas * o, puls, potentiel, energie);	// l4
}

// Renvoie la fonction d'onde
double
PSI (double psi, double m, double n, double o, double p, double pas)
{
  psi = psi + (pas / 6) * (m + 2 * n + 2 * o + p);
  return psi;
}

// Renvoie la derivee
double
PHI (double phi, double a, double b, double c, double d, double pas)
{
  phi = phi + (pas / 6) * (a + 2 * b + 2 * c + d);
  return phi;
}










// Calcul le potentiel et le vecteur d'onde en fonction de choix_pot
void calcul_potentiel (double x, double &v, double V0, double &vo, double masse,
		  long double e, double h_barre, double eV, double largeur)
{
  v = x*x;
  vo = eV * sqrt (2 * masse * abs (v - e)) / h_barre;
}

void calcul (double valeurs[][4], double x, double psi, double phi, double v,
	double V0, double vo, double masse, long double e, double h_barre,
	double eV, double h, double largeur, int points)
{
  // Coefficients pour la methodes RK4 :
  double a, b, c, d, m, n, o, p;

  // Boucle d'ecriture
  for (int i = 0; i < points; i++)
    {
      calcul_potentiel (x, v, V0, vo, masse, e, h_barre, eV, largeur);

      // On ecrit position, fonction d'onde, derivee et potentiel dans un tableau
      valeurs[points-1 - i][0] = x;
      valeurs[points-1 - i][1] = psi;
      valeurs[points-1 - i][2] = phi;
      valeurs[points-1 - i][3] = v;

      // Calcul coefficients RK4
      coeff (psi, phi, a, b, c, d, m, n, o, p, h, vo, v, e);
      // Calcul fonction d'onde et derivee
      psi = PSI (psi, m, n, o, p, h);
      phi = PHI (phi, a, b, c, d, h);

      // Increment position
      x = x + h;
    }
}




























// Ecriture de l'explosion en fonction de l'energie
void ecriture_tableau (long double explosion_min[][2])
{
  fstream f;

  // ON ECRIT L'EXPLOSION EN FONCTION DE L'ENERGIE
  f.open ("OH_explosion.txt", ios::out);
  for (int i = 0; i < 1000; i++)
    {
      f << explosion_min[i][0] << " " << explosion_min[i][1] << endl;
    }
  f.close ();

  // ON ECRIT LA VALEUR ABSOLUE DE L'EXPLOSION POUR FAIRE DES GRAPHES LISIBLES
  f.open ("OH_explosion_ABS.txt", ios::out);
  for (int i = 0; i < 1000; i++)
    {
      f << explosion_min[i][0] << " " << abs (explosion_min[i][1]) << endl;
    }
  f.close ();
}











// On ramene le maximum de la fonction d'onde a 0.25
void normalisation (double valeurs[][4], int points)
{
  double maximum = valeurs[0][1];
  for (int i = 0; i < points; i++)
    {
      if (abs (valeurs[i][1]) > abs (maximum))
	{
	  maximum = valeurs[i][1];
	}
    }
  for (int i = 0; i < points; i++)
    {
      valeurs[i][1] = valeurs[i][1] / (4 * maximum);
    }
}

// ecriture de la fonction d'onde dans un fichier
void ecriture (double valeurs[][4], long double e, int points)
{
  // On normalise avant d'ecrire dans un fichier
  normalisation (valeurs, points);
  // On ecrit dans un fichier
  fstream f;
  f.open ("OH_fct.txt", ios::out);
  for (int i = 0; i < 100; i++)
    {
      f << valeurs[i][0] << " " << valeurs[i][1] +
	e << " " << e << " " << valeurs[i][3] << endl;
    }
  f.close ();
}








void
balayage (double V0, double largeur, int points)
{
  double h = 0.01;		// Pas utilise pour la methode RK4 (negatif car on part de x=6 et on remonte jusqu'a x=0)

  //CONSTANTES
  double masse = 9.109 * pow (10, -31);
  double h_barre = 6.62607015 * pow (10, -34) / (2 * M_PI);
  double eV = 1.60217657 * pow (10, -19);


  // Creaction du tableau pour garder les valeurs position, fonction d'onde, derivee et potentiel en memoire
  double valeurs[points][4];

  double Emin, Emax;
  cout << endl << "De combien a combien souhaitez-vous balayer ?" << endl;
  cin >> Emin;
  cin >> Emax;

  long double balayage = (Emax - Emin) / 1000;
  long double e = Emin;

  long double explosion_min[1000][2];

  double vo;
  double x;
  double psi;
  double phi;
  double v;


  // BALAYAGE EN ENERGIE
  for (int i = 0; i < 1000; i++)
    {
      // On recalcule les conditions initiales a chaque iteration de la propagation
      vo = eV * sqrt (2 * masse * abs (V0 - e)) / h_barre;
      x = -5;
      psi = 0;
      phi = vo * exp (vo * x);

      // On ecrit la valeur de l'energie et l'explosion en x=0 dans un tableau
      calcul (valeurs, x, psi, phi, v, V0, vo, masse, e, h_barre, eV, h, largeur, points);
      explosion_min[i][0] = e;
      explosion_min[i][1] = valeurs[0][1];

      e += balayage;
    }
  cout << endl;


  // ON AFFICHE LE BALAYAGE ET ON L'ECRIT DANS UN FICHIER
  for (int i = 0; i < 1000; i++)
    {
      cout << "Energie = " << explosion_min[i][0] << "  Explosion = " <<
	explosion_min[i][1] << endl;
    }
  ecriture_tableau (explosion_min);




  // ON COMPTE LES VALEURS PROPRES DES ETATS LIES
  int j = 0;			// compte les changements de signe de l'explosion en x = 0
  double signe, change;

  if (explosion_min[0][1] < 0)
    {
      change = -1;
    }
  else
    {
      change = +1;
    }
  for (int i = 1; i < 1000; i++)
    {
      if (explosion_min[i][1] < 0)
	{
	  signe = -1;
	}
      else
	{
	  signe = +1;
	}

      if (signe != change)
	{
	  j++;
	  change = signe;
	}
    }
  cout << endl << "Nombre de passage par 0 : " << j << endl << endl;

  // ON ECRIT LES INTERVALLES DES VALEURS PROPRES DANS UN TABLEAU
  double vap[2 * j][2];
  int n = 0;
  if (explosion_min[0][1] < 0)
    {
      change = -1;
    }
  else
    {
      change = +1;
    }
  for (int i = 1; i < 1000; i++)
    {
      if (explosion_min[i][1] < 0)
	{
	  signe = -1;
	}
      else
	{
	  signe = +1;
	}

      if (signe != change)
	{
	  vap[n][0] = explosion_min[i - 1][0];
	  vap[n][1] = explosion_min[i - 1][1];
	  n++;
	  vap[n][0] = explosion_min[i][0];
	  vap[n][1] = explosion_min[i][1];
	  n++;
	  change = signe;
	}
    }
  // ON AFFICHE LES VALEURS PROPRES DES ETATS LIES ET ON LES ECRIT DANS UN TABLEAU
  fstream f;
  f.open ("Valeurs_propres_OH.txt", ios::out);
  for (int i = 0; i < 2 * j; i++)
    {
      cout << "Energie = " << vap[i][0] << "  Explosion = " << vap[i][1] <<
	endl;
      f << vap[i][0] << " " << vap[i][1] << endl;
    }
  f.close ();


  // tri du tableau pour trouver les valeurs propres de l'energie
  double a, b, c, d;

  for (int i = 0; i < 1000; i++)
    {
      for (int j = i + 1; j < 1000; j++)
	{
	  if (abs (explosion_min[i][1]) > abs (explosion_min[j][1]))
	    {
	      a = explosion_min[i][0];
	      b = explosion_min[i][1];
	      c = explosion_min[j][0];
	      d = explosion_min[j][1];

	      explosion_min[i][0] = c;
	      explosion_min[i][1] = d;
	      explosion_min[j][0] = a;
	      explosion_min[j][1] = b;
	    }
	}
    }

  // On connait maintenant la valeur de l'energie pour laquelle l'explosion est minimale, on genere les donnees pour cette energie
  e = explosion_min[0][0];
  cout << endl << "explosion minimale : e = " << e << endl;
  double choix;
  cout << endl << "Voulez-vous ecrire une fonction d'onde ?" << endl <<
    "1. Oui" << endl << "2. Non" << endl << endl;
  cin >> choix;
  if (choix == 1)
    {
      cout << "Entrer une valeur pour l'energie :" << endl << endl;
      cin >> e;
      vo = eV * sqrt (2 * masse * abs (V0 - e)) / h_barre;
      x = 6;
      psi = 0;
      phi = vo * exp (vo * x);

      // On ecrit la valeur de l'energie et l'explosion en x=0 dans un tableau
      calcul (valeurs, x, psi, phi, v, V0, vo, masse, e, h_barre, eV, h, largeur, points);
      ecriture (valeurs, e, points);
    }
}























int main ()
{
  int points = 1000;      // Discrétisation de l'axe X
  double V0;			// Valeur maximale du potentiel (en eV)
  double largeur;		// Largeur du potentiel

  /*choix_initial (V0, largeur);*/
  balayage (V0, largeur, points);




  return 0;
}
