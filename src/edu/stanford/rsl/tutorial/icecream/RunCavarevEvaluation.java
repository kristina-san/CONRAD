package edu.stanford.rsl.tutorial.icecream;

import java.io.File;
import java.util.Scanner;

public class RunCavarevEvaluation {

	public static void main(String[] args) throws Exception {
	
		CavarevEvaluation eval;
		System.out.println("Cavarev Evaluation v1.0");
		
		Scanner sc = new Scanner(System.in);
	    System.out.print("Please enter the path to your file that should be evaluated: ");
	    String evaluationFile = sc.next();
	    
	    File file = new File(evaluationFile);
	    if (!file.exists()) {
	       System.out.print("Your entered filepath doesn't exist. Abort program.");
	    }
	    else {
	    	System.out.println("Please enter if your file is card (c) oder cardbreath (b)");
	    	String mode = sc.next();
		    
		    switch (mode) {
		    case "b":	eval = new CavarevEvaluation(evaluationFile, true);
		    			eval.run();
		    			break;
		    case "c":	eval = new CavarevEvaluation(evaluationFile, false);
		    			eval.run();
		    			break;
		    default:	System.out.print("False input. Abort program.");
		    }
	    }
	    sc.close();
	}
}
