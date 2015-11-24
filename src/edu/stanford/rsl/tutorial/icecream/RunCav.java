package edu.stanford.rsl.tutorial.icecream;

import java.awt.Button;
import java.awt.Label;
import java.awt.Panel;
import java.awt.TextArea;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import ij.gui.GenericDialog;

public class RunCav extends GenericDialog {
	
	private TextField tf;				// text field for entering the path to the evaluation file
	private TextArea resultTextArea;	// text area for error or output messages 
	private String selection;			// get selection from RadioButton
	
	public RunCav(String title) {
		super(title);
	}

	public void run() throws Exception {

		Panel topPanel = new Panel();
	    tf = new TextField("", 50);	// text field for entering the path to the evaluation file
		Button button_chooseFile = new Button("open");	// button that opens a JFileChooser so the evaluation file could be choosed
		button_chooseFile.addActionListener(new ActionListener() {
	         public void actionPerformed(ActionEvent e) {
	 			JFileChooser chooser = new JFileChooser();
	 		    FileNameExtensionFilter filter = new FileNameExtensionFilter(
	 		           "Binary files", "bin");	//show only bin files
	 		    chooser.setFileFilter(filter);
				int returnVal = chooser.showOpenDialog(getParent());
				if(returnVal == JFileChooser.APPROVE_OPTION)
				{
					tf.setText(chooser.getSelectedFile().getPath());	//show the path to the file in text field
				}
	         }
	      });

		Label label = new Label("Filepath: ");
		topPanel.add(label);
		topPanel.add(tf);
		topPanel.add(button_chooseFile, null, 2);
	    this.addPanel(topPanel);	
	    
	    // Add radio buttons for choosing between cardbreath and card
	    String[] labels = {"Cardbreath", "Card"};    
	    this.addRadioButtonGroup("Please check", labels, 1, 1, labels[0]);
	   	
	    Panel bottomPanel = new Panel();
	    bottomPanel.add(Box.createHorizontalGlue());
	    resultTextArea = new TextArea();	// text area for error or output messages 
	    resultTextArea.setEditable(false);
	    bottomPanel.add(resultTextArea);
	    
	    // after clicking on button "Run" the insert file path is checked 
	    // and then the cavarev evaluation is started 
	    Button run = new Button("Run");
	    run.addActionListener(new ActionListener() {
	         public void actionPerformed(ActionEvent e) {
	        	String evaluationFile = tf.getText();
	    	    boolean isOkay = check(evaluationFile);

    	    	String evaluationString = "";
	    	    if(isOkay)
	    	    {
		        	selection = getTextRadioButton();
					try {
						evaluationString = evaluate(evaluationFile, selection);
					} catch (Exception e1) {
						e1.printStackTrace();
					}
					resultTextArea.setText("Finished evaluation.\n" + evaluationString);;
	    	    }	        	 
	         }
	      });
	    bottomPanel.add(run);
	    this.addPanel(bottomPanel);
	    this.showDialog();
	    if (this.wasCanceled()) return;
	}
	
	//get the text from the radio button
	public String getTextRadioButton() {
		return this.getNextRadioButton();
	}
	
	// check the insert file path isn't empty, if the file exists and if it's a binary file
	public boolean check(String evaluationFile) {
		if(evaluationFile.isEmpty()) {
			resultTextArea.setText("You forgot the filepath. Please insert and try again.");
			return false;
		} else {
			File file = new File(evaluationFile);
			if (!file.exists()) {
				resultTextArea.setText("The insert file doesn't exist. Please try again.");
				return false;
			}
			else if (!evaluationFile.endsWith(".bin"))
			{
				resultTextArea.setText("The insert file isn't a binary file. Please try again.");
				return false;
			}
		} 
		return true;
	}
	
	// start cavarev evaluation
	public String evaluate(String evaluationFile, String selection) throws Exception {
		CavarevEvaluation eval = null;
		String evaluationString = "";
		resultTextArea.setText("Started evaluation... \nProgram is running...\nPlease see console output for more information.");
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
