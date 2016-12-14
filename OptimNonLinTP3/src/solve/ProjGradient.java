package solve;
import line.Dichotomy;
import line.LineSearch;

import java.util.ArrayList;

import func.RealFunc;
import util.Matrix;
import util.Singularity;
import util.Vector;


/**
 * The projected gradient algorithm for non-linear optimization under 
 * linear inequalities.
 * 
 * @author G Chabert
 *
 */
public class ProjGradient extends Algorithm {
	
	/**
	 *  Minimal step for the line search
	 */

	final static double MIN_ALPHA=1e-20;
	/**
	 *  Minimal norm of the gradient under which Khun-Tucker 
	 *  conditions are applied.
	 */
	public final static double MIN_GRAD_NORM=1e-10;
	
	/**
	 * Threshold under which an equality is considered as satisfied. 
	 *
	 * When for a constraint g<=0, |g(x)|<ALMOST_ZERO, the
	 * constraint is considered to be active. 
	 */
	public final static double ALMOST_ZERO=1e-08;

	// Number of variables
	int n;
	// Number of constraints
	int m;
	// Matrix of constraints
	Matrix A;
	// Right-hand vectors of constraints
	Vector b;
	// Function to minimize
	RealFunc f;
	// Line search
	LineSearch line;
	// Dichotomie pour alphaMax
	Dichotomy dich;
	// Le point courant
	Vector iter_vec;
	// Les contraintes actives
	ArrayList<Boolean> activeConstr;
	// Direction de descente
	Vector d;
		

	/**
	 * Check that x is feasible (the method raises and error if x is infeasible).
	 * 
	 * The projected gradient algorithm only yield feasible vectors.
	 * Call this method at each iteration, to check that your code is correct.
	 * 
	 * @param x - the current vector
	 * @throws Exception 
	 */
	private void check_feasible(Vector x) throws Exception {
		int i;
		boolean feasible = true;
		Vector res = new Vector(A.mult(x));
		
		/*System.out.println("A : \n" + A.toString());
		System.out.println("x : \n" + x.toString());
		System.out.println("res : \n" + res.toString());
		
		System.out.println("Comparaisons");*/
		i = 0;
		while (feasible && i < m) {
			//System.out.println(res.get(i) + " <= " + b.get(i));
			feasible = (res.get(i) <= (b.get(i) + ALMOST_ZERO));
			++i;
		}
		
		if (!feasible) {
			//System.out.println("Unfeasible");
			throw (new Exception("Le point n'est pas faisable : " + x.toString()));
		} else {
			//System.out.println("Feasible");
		}
	}
	
	
	/**
	 * This constructor:
	 * - build the projected gradient algorithm for minimizing f
	 *   under constraints Ax<=b. 
	 * - check that the initial point x0 is feasible.
	 * - determine the initial set of active constraints.
	 * @throws Exception 
	 * 
	 */
	public ProjGradient(Matrix A, Vector b, RealFunc f, LineSearch l, Vector x0) throws Exception {
		this.A = A;
		this.b = b;
		this.m = A.nb_rows();
		this.n = A.nb_cols();
		this.f = f;
		this.line = l;
		this.dich = new Dichotomy(f);
		this.iter_vec = x0;
		this.check_feasible(x0);
		activeConstr = new ArrayList<Boolean>(m);
		for (int i = 0; i < m; ++i) {
			activeConstr.add(i, false);
		}
		setActiveConstraints();
		
		
		
		// tests
		setDirection();
		
		
		System.out.println("direction : " + d);
		System.out.println("alpha : " + rechercheEnLigne());
	}
	
	/**
	 * Check the active constraints
	 * @return the list of booleans, true if the constraint is active, false otherwise
	 */
	private void setActiveConstraints() {
		Vector res = new Vector(A.mult(iter_vec));
		
		for (int i = 0; i < m; ++i) {
			activeConstr.set(i, (res.get(i) >= (b.get(i) - ALMOST_ZERO)));
		}
	}
	
	/**
	 * Set the direction of the steepest descent of the current function at the current point
	 */
	private void setDirection() {
		d = f.grad(iter_vec).minus();
	}
	
	private double rechercheEnLigne() {
		double alpha, alphaMax, alphaCourant, resultat;
		
		alpha = line.search(iter_vec, d);
		alphaMax = 0;
		
		for (int i = 0; i < m; ++i) {
			if (!activeConstr.get(i)) {
				alphaCourant = A.get_row(i).scalar(d);

				if (alphaCourant != 0) {	// Ou juste différent de 0 et alpha négatif signifie que la direction n'est pas admissible ?
					//System.out.print("alphaCourant = " +b.get(i) + " / " + alphaCourant + " = ");
					alphaCourant = b.get(i) / alphaCourant;
					//System.out.println(alphaCourant);
					if ((alphaCourant >= 0) && (alphaCourant < alphaMax)) {
						alphaMax = alphaCourant;
					}
				}
			}
		}
		
		System.out.println("alphaMax = " + alphaMax);
			
		if ((alpha >=0) && (alpha < alphaMax)) {
			System.out.println("Alpha < alphaMax");
			resultat = alpha;
		} else if (f.eval(iter_vec.add(d.leftmul(alphaMax))) < f.eval(iter_vec)) {
			System.out.println("alphaMax diminue f");
			resultat = alphaMax;
		} else {
			System.out.println("Dichotomie nécessaire");
			resultat = dich.search(iter_vec, d);
		}
		
		return resultat;
	}
	
	/**
	 * Calculate the next iterate.
	 * 
	 * @return x (the same reference) if the solution is reached.
	 */
	public void compute_next() throws EndOfIteration {
		/*
La boucle principale de la m ́ethode du gradient projet ́e (next) calculant x k+1 `a partir de x k consiste en les
 ́etapes suivantes :
1. Calculez une direction de descente admissible d
2. Si d ∼ 0, calculez les multiplicateurs de Lagrange (basez-vous sur le cours).
3. Si au moins l’un d’entre eux est n ́egatif, d ́esactivez la contrainte ayant le plus grand multiplicateur en
valeur absolue et recommencez la proc ́edure (nouveau choix de direction d).
4. Sinon, un minimum a  ́et ́e trouv ́e.
5. Sinon (si d est bien diff ́erent de 0) effectuez une recherche en ligne.
		 */
		setDirection();
		
	}

}
