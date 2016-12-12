package func;
import util.Vector;

/**
 * An example of function extracted from "Handbook of Test Problems in Local and Global Optimization" (Floudas et al.)
 * 
 * @author Gilles Chabert
 */
public class FuncTp3 implements RealFunc {

	/**
	 * @return f(x).
	 */
	public double eval(Vector x) {
		return (-0.5*(10*x.get(0)*x.get(0) + 10*x.get(1)*x.get(1) + 10*x.get(2)*x.get(2) + 10*x.get(3)*x.get(3) + 10*x.get(4)*x.get(4) + 10*
				x.get(5)*x.get(5) + 10*x.get(6)*x.get(6))) - 20*x.get(0) - 80*x.get(1) - 20*x.get(2) - 50*x.get(3) - 60*x.get(4) - 90*
				x.get(5) + 10*x.get(7) + 10*x.get(8) + 10*x.get(9);
	}


	/**
	 * @return the gradient of f at x.
	 */
	public Vector grad(Vector x) {		
		double[] _g = 
			{   -10*x.get(0)-20,
				-10*x.get(1)-80,
				-10*x.get(2)-20,
				-10*x.get(3)-50,
				-10*x.get(4)-60,
				-10*x.get(5)-90,
				-10*x.get(6),
				10,
				10,
				10
			};
		return new Vector(_g);
	}


	/**
	 * @return 10.
	 */
	@Override
	public int dim() {
		return 10;
	}

}
