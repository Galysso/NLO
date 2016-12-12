package main;
import line.SlopeIter;
import func.FuncTp3;
import util.Matrix;
import util.Vector;
import solve.ProjGradient;

public class Tp3 {


	public static void main(String[] args) throws Exception {

		double[][] _A = 
				// =============== constraints ===================
				{{-2 , -6 , -1 , 0 , -3 , -3 , -2 , -6 , -2 , -2},
				{6 , -5 , 8 , -3 , 0 , 1 , 3 , 8 , 9 , -3},
				{-5 , 6 , 5 , 3 , 8 , -8 , 9 , 2 , 0 , -9},
				{9 , 5 , 0 , -9 , 1 , -8 , 3 , -9 , -9 , -3},
				{-8 , 7 , -4 , -5 , -9 , 1 , -7 , -1 , 3 , -2},
				{-7 , -5 , -2 , 0 , -6 , -6 , -7 , -6 , 7 , 7},
				{1 , -3 , -3 , -4 , -1 , 0 , -4 , 1 , 6 , 0},
				{1 , -2 , 6 , 9 , 0 , -7 , 9 , -9 , -6 , 4},
				{-4 , 6 , 7 , 2 , 2 , 0 , 6 , 6 , -7 , 4},
				{1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1},
				{-1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1},
				// =============== for xi<=1 ======================
				{1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
				{0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
				{0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
				{0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
				{0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
				{0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
				{0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
				// =============== for xi>=0 ======================				
				{-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
				{0, -1, 0, 0, 0, 0, 0, 0, 0, 0},
				{0, 0, -1, 0, 0, 0, 0, 0, 0, 0},
				{0, 0, 0, -1, 0, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, -1, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, -1, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 0, -1, 0, 0, 0},
				{0, 0, 0, 0, 0, 0, 0, -1, 0, 0},
				{0, 0, 0, 0, 0, 0, 0, 0, -1, 0},
				{0, 0, 0, 0, 0, 0, 0, 0, 0, -1},
				};
		
		double[] _b={
				// ======= real constraints =======
				-4,22,-6,-23,-12,-3,1,12,15,9,-1,
				// ========== xi<=1 ===============
				1,1,1,1,1,1,1,1,1,1,
				// ========== xi>=0 ===============
				0,0,0,0,0,0,0,0,0,0
		};
		
		
		// Initial point (the solution actually minimizing x1-x10)
		double[] _x0={0 , 0 , 0 , 0.532416502935 , 0.801571709255 , 0.876227897868 , 0 , 1 , 0 , 1};
			
		Matrix A = new Matrix(_A);
		Vector b = new Vector(_b);
		FuncTp3 f = new FuncTp3();
		ProjGradient pg = new ProjGradient(A, b, f, new SlopeIter(f), new Vector(_x0));
		
		/* TODO */
	}
}
