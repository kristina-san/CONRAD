package edu.stanford.rsl.tutorial.icecream;

import java.awt.Button;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JFileChooser;

import ij.gui.GenericDialog;

public class RunCav extends GenericDialog implements ActionListener {
	
	private Button button;

	public RunCav(String title) {
		super(title);
	}

	public void run() throws Exception {
		String evaluationString = "";
		String evaluationFile = "";
		String selection = "";
		Button button = new Button("open");
	    button.addActionListener(new ActionListener() {
	         public void actionPerformed(ActionEvent e) {
	 			JFileChooser chooser = new JFileChooser();
				chooser.showDialog(null, "Open");

				int returnVal = 0;
				if(returnVal == JFileChooser.APPROVE_OPTION)
				{
			         //txt.setText();
					System.out.println(chooser.getSelectedFile().getPath());
				}
	         }
	      });
	    this.addStringField("Filepath: ", null, 50);	// text von oben noch setzen!
	    this.add(button, null, 2);
	    
	    String[] labels = {"Cardbreath", "Card"};    
	    this.addRadioButtonGroup("Please check", labels, 1, 1, labels[0]);
	    
	    this.showDialog();
	    if (this.wasCanceled()) return;
	    
	    evaluationFile = this.getNextString();
	    selection = this.getNextRadioButton();
	    
	    boolean isOkay = check(evaluationFile, selection);
	    if(isOkay)
	    {
	    	evaluationString = evaluate(evaluationFile, selection);
		    this.addMessage(evaluationString);
	    }

	}
	
	public boolean check(String evaluationFile, String selection) {
		if(evaluationFile.isEmpty()) {
			//error message that string is empty
			return false;
		} else {
			File file = new File(evaluationFile);
			if (!file.exists()) {
				//error message that file doesn't exist
				return false;
			}
		}
		return true;
	}
	
	public String evaluate(String evaluationFile, String selection) throws Exception {
		CavarevEvaluation eval = null;
		String evaluationString = "";
    	switch (selection) {
			case "Card":
				eval = new CavarevEvaluation(evaluationFile, false);
				break;
			case "Cardbreath":	
				eval = new CavarevEvaluation(evaluationFile, true);
				break;
    	}
    	evaluationString = eval.run();
    	
    	return evaluationString;
	}
	
	public static void main(String [] args) throws Exception {
		
		RunCav rv = new RunCav("Run cavarev");
		rv.run();
		
	}
	
}
