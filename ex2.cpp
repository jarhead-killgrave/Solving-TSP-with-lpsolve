#include <lpsolve/lp_lib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <cstdlib>

using namespace std;

/**
 * Nom et Prenom : KITSOUKOU Manne Emile
 * Numero etudiant : 22013393
 * Parcours: L2 informatique 2B
 * 
 */


/**
 * La fonction permet de vider un tableau 1d
 * 
 * @param row le tableau à vider
 * @param nbVariable le nombre de case à vider
 */
void vider(REAL row[], int nbVariable)
{
    for (int i = 0; i < nbVariable; i++)
    {
        row[i] = 0.0;
    }
}

/**
 * Calcule la distance eucledienne au carré entre 2 points
 * 
 * @param x1 l'abscisse en x du premier point
 * @param y1 l'abscisse en y du premier point
 * @param x2 l'abscisse en x du deuxieme point
 * @param y2 l'abscisse en y du deuxieme point
 * @return double qui represente la distance au carré entre les 2 points
 */
double distance_euclidean_au_carre(double x1, double y1, double x2, double y2)
{   // On utilise la formule de la distance euclidienne au carré : (x1 - x2)^2 + (y1 - y2)^2
    return pow(x1 - x2, 2) + pow(y1 - y2, 2);
}

/**
 * Charge un fichier qui respecte les reglementation tsp
 *
 * @param nom_fichier le nom du fichier à charger
 * @param n  un pointeur qui stockera l'addresse de la valeur de la taille des données récupérés;
 * @return double** un tableau 2D contenant les données récupérés
 */
double **chargement_donnees(char *nom_fichier, size_t *n)
{

    double **result;

    ifstream fichier(nom_fichier, ios::in);

    if (fichier)
    {
        int position;
        double x;
        double y;
        int debut_donnees = 0;

        string ligne;
        string chaine;

        while (getline(fichier, ligne) && ligne != "EOF")
        {
            stringstream s(ligne);
            s >> chaine;

            if (debut_donnees == 1)
            {
                stringstream s(ligne);
                s >> position >> x >> y;

                result[position - 1] = (double *)malloc(sizeof(double) * 2);
                result[position - 1][0] = x;
                result[position - 1][1] = y;
            }
            else if (chaine == "NODE_COORD_SECTION")
            {
                debut_donnees = 1;
            }
            else if (chaine == "DIMENSION:")
            {
                stringstream s(ligne);
                s >> chaine >> (*n);
                result = (double **)malloc(sizeof(double *) * (*n));
            }
            
        }

        fichier.close();
    }

    return result;
}

int parcours_minimal(double **coords, int n)
{

    ///////////////////////////////// Initialisation des variables utiliser dans le probleme //////////////////////////////////////////////////
    
    // Declaration du nombre de variables du probleme
    int nbVariable = n * n;

    // Tableau qui stockent les coefficient des contraintes
    REAL row[nbVariable + 1];

    // Initialisation de la structure lprec
    lprec *lp = make_lp(0, nbVariable);


    /////////////////////////////////// Debut de l'implementation du probleme ///////////////////////////////////////////////////////////////

    // Tous les Xij sont compris entre 0 et 1
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                set_binary(lp, i * n + j + 1, TRUE);
            }
            else {
                // Si i = j
                set_int(lp, i * n + i + 1, true);
            }
        }
    }

    // Ajout de toute les contraintes 1 au probleme
    for (int i = 0; i < n; i++)
    {
        vider(row, nbVariable + 1);
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                row[i * n + j + 1] = 1;
            }
        }
        add_constraint(lp, row, EQ, 1);
    }

    // Ajout de toutes les contraintes 2 au probleme
    for (int i = 0; i < n; i++)
    {
        vider(row, nbVariable + 1);
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                row[j * n + i + 1] = 1;
            }
        }
        add_constraint(lp, row, EQ, 1);
    }



    // Ajout de toute les contraintes 3 au probleme
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < n; j++)
        {
            if (i != j)
            {
                vider(row, nbVariable + 1);
                row[i * n + i + 1] = 1;
                row[j * n + j + 1] = -1;
                row[i * n + j + 1] = n;
                add_constraint(lp, row, LE, n - 1);
            }
        }
    }

    // Ajout de toute les contraintes 4 au probleme
    // La contrainte 4 revient juste à fixer les bornes de des valeurs prises par les variables Ui
    for (int i = 1; i < n; i++)
    {

        set_lowbo(lp, i * n + i + 1, 1);
        set_upbo(lp, i * n + i + 1, n - 1);
    }

    // Fonction objective
    vider(row, nbVariable + 1);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                row[i * n + j + 1] = distance_euclidean_au_carre(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
            }
        }
    }
    set_obj_fn(lp, row);

    int max = 0;
    for (int i = 1; i <= n; i++)
    {
        max += distance_euclidean_au_carre(coords[i % n][0], coords[i % n][1], coords[i-1][0], coords[i-1][1]);
    }
    cout << "max = " << max << endl;

    // Fixe le nom du probleme
    char nom_model[] = "Raccordement des villes avec la fibre optique";
    set_lp_name(lp, nom_model);

    // Spécifie la branche à prendre en premier dans l'algorithme branch-and-bound.
    // On evaluera d'abord la borne minimal pour  ́eviter d’explorer inutilement certaines branches
    set_bb_floorfirst(lp, BRANCH_FLOOR);

    // On utilisera la methode d'evaluation PRICER_STEEPESTEDGE qui parait efficace pour cet probleme
    set_pivoting(lp, PRICER_STEEPESTEDGE);

    // On suppossera que la valeur de la fonction objective ne peut pas avoir une valeur superieur à max
    // max est la valeur du fonction objective si l'on va successivement de la ville 1 à la ville n dans unc cycle ferme
    set_obj_bound(lp, max);
    set_basiscrash(lp, CRASH_MOSTFEASIBLE);
    int res = solve(lp);
    


    ////////////////////////////////////////////// Affichage et Ecriture des solutions //////////////////////////////////////////////////

    // Si on a une solution
    if (res == 0)
    {
        string nomfichier("solution.tmp");

        ofstream f("solution.dot", ios::out | ios::trunc);
        if (f)
        {
            get_variables(lp, row);
            f << "digraph {\n" << endl;
            
            for (int i = 0; i < n; i++)
            {
                f << "\t" << (i + 1) << ";\n" << endl;
            }
            f << "\n" << endl;

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i != j && row[i * n + j] == 1.0)
                    {
                        f << "\t" << (i + 1) << " -> " << (j + 1) << ";\n" << endl;

                        cout << "Ville_" << i+1 << " -- Ville_" << j+1 << endl;
                    }
                }
            }
            f << "}" << endl;
            f.close();
            return 1;
        }
    }
    return 0;
}



int main(int argc, char *argv[])
{
    size_t n;
    int res;
    int nb_villes_a_etudier;
    double **tab;
    char * nom_fichier_donnees;
    char *nom_fichier_resultat;
    // Si nous avons 3 arguments de la forme <nom_excecutable> __ <fichier_entree_format_tsp> __ <fichier_sortie_format.dot> __ <nb_villes_a_etudier>
    if (argc == 4)
    {   
        nom_fichier_donnees = argv[1];
        nom_fichier_resultat = argv[2];
        nb_villes_a_etudier = stoi(argv[3]);
    }
    
    tab = chargement_donnees(nom_fichier_donnees, &n);

    if(nb_villes_a_etudier > 12){
        cout << "Vous avez lancé la resolution sur un nombre de ville superieur à 12" << endl;
        cout << "Cela peut prendre un certain temps. Pour une resolution rapide veuillez relancer avec au maximun 12 ville" << endl;
    }

    if (nb_villes_a_etudier <= n)
    {
        res = parcours_minimal(tab, nb_villes_a_etudier);
    }
    else
    {
        cout << "Vous avez spécifié un nombre de ville trop élévé par rapport au données" << endl;
    }

    if (res == 1) {
        cout << "Pas d'erreur lors de la resolution vous trouverez la solution également dans le fichier solution.dot ou celui dans vous avez proposer en parametre" << endl;
    }
    return 0;
}