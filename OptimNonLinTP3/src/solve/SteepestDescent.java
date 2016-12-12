package solve;
import line.Dichotomy;
import line.LineSearch;
import func.RealFunc;
import util.Vector;

/**
 * Basic steepest descent algorithm for unconstrained minimization problem.
 * 
 * @author Gilles Chabert
 */

public class SteepestDescent extends Algorithm {
	
	public static int count;
	
	private RealFunc f;
	private LineSearch s;
	private Dichotomy dicho;
	
	/**
	 * Build the algorithm
	 * 
	 * @param f function to minimize
	 * @param s underlying line search algorithm
	 */
	public SteepestDescent(RealFunc f, LineSearch s) {
		this.f = f;
		this.s = s;
		this.dicho = new Dichotomy(f);
	}

	
	/**
	 * Calculate the next iterate.
	 * 
	 */
	public void compute_next() throws EndOfIteration {
	
		Vector mg=f.grad(iter_vec).leftmul(-1);
		double alpha=s.search(iter_vec, mg);
		
		if (f.eval(iter_vec)<=f.eval(iter_vec.add(mg.leftmul(alpha)))) {
			if (log) System.out.println("[steepest] slope iter failed: try dichotomy...");
			alpha=dicho.search(iter_vec, mg);
		}
		iter_vec=iter_vec.add(mg.leftmul(alpha));
	}
	
}
