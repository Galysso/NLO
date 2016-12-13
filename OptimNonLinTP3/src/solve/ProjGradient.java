package solve;
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
		this.iter_vec = x0;
		this.check_feasible(x0);
		activeConstr = new ArrayList<Boolean>(m);
		for (int i = 0; i < m; ++i) {
			activeConstr.add(i, false);
		}
		setActiveConstraints();
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
	
	/**
	 * Calculate the next iterate.
	 * 
	 * @return x (the same reference) if the solution is reached.
	 */
	public void compute_next() throws EndOfIteration {
		// 1)
		double alpha = line.search(iter_vec, d);
		
		// 2)
		
	}

}
