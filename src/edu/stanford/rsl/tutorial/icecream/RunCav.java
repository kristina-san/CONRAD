package edu.stanford.rsl.tutorial.icecream;

import java.awt.Button;
import java.awt.Label;
import java.awt.Panel;
import java.awt.TextArea;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JFileChooser;

import ij.gui.GenericDialog;

public class RunCav extends GenericDialog implements ActionListener {
	
	private TextField tf;
	private TextArea resultTextArea;
	
	private String selection;
	public RunCav(String title) {
		super(title);
	}

	public void run() throws Exception {

		Panel panel = new Panel();
	    tf = new TextField("", 50);
		Button button = new Button("open");
	    button.addActionListener(new ActionListener() {
	         public void actionPerformed(ActionEvent e) {
	 			JFileChooser chooser = new JFileChooser();
				chooser.showDialog(null, "Open");

				int returnVal = 0;
				if(returnVal == JFileChooser.APPROVE_OPTION)
				{
					tf.setText(chooser.getSelectedFile().getPath());
				}
	         }
	      });
	    
		Label label = new Label("Filepath: ");
	    panel.add(label);
	    panel.add(tf);
	    panel.add(button, null, 2);
	    this.addPanel(panel);
	    
	    String[] labels = {"Cardbreath", "Card"};    
	    this.addRadioButtonGroup("Please check", labels, 1, 1, labels[0]);
	   	    
	    resultTextArea = new TextArea();
	    resultTextArea.setEditable(false);
	    panel.add(resultTextArea);
	    
	    Button run = new Button("Run");
	    run.addActionListener(new ActionListener() {
	         public void actionPerformed(ActionEvent e) {
	     	    String evaluationFile = tf.getText();
	    	    selection = getTextRadioButton();
	    	    
	    	    boolean isOkay = check(evaluationFile, selection);
    	    	String evaluationString = "";
	    	    if(isOkay)
	    	    {
					try {
						evaluationString = evaluate(evaluationFile, selection);
					} catch (Exception e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}
					resultTextArea.setText("Finished evaluation.\n" + evaluationString);;
	    	    }	        	 
	         }
	      });
	    panel.add(run);
	    
	    this.showDialog();
	    if (this.wasCanceled()) return;
	}
	
	public String getTextRadioButton() {
		return this.getNextRadioButton();
	}
	
	public boolean check(String evaluationFile, String selection) {
		if(evaluationFile.isEmpty()) {
			resultTextArea.setText("You forgot the filepath. Please insert and try again.");
			return false;
		} else {
			File file = new File(evaluationFile);
			if (!file.exists()) {
				resultTextArea.setText("The insert file doesn't exist. Please try again.");
				return false;
			}
		}
		return true;
	}
	
	public String evaluate(String evaluationFile, String selection) throws Exception {
		CavarevEvaluation eval = null;
		String evaluationString = "";
		resultTextArea.setText("Started evaluation... \nPlease see console for more information.");
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
