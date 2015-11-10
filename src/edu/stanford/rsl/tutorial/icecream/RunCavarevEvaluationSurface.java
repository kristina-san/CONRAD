package edu.stanford.rsl.tutorial.icecream;

import ij.gui.GenericDialog;

import java.awt.Checkbox;
import java.awt.CheckboxGroup;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

public class RunCavarevEvaluationSurface extends JFrame implements ActionListener, ItemListener  {
	
	private static final long serialVersionUID = 1L;
	JLabel label0;
	JLabel label;

	JPanel panel;

	JTextField txt;
	
	JButton button;
	JButton button1;
	
	CheckboxGroup checkboxgroupRadio;
	Checkbox card;
	Checkbox cardBreath;
	
	public RunCavarevEvaluationSurface() {
		this.setTitle("Run Cavarev");
		this.setSize(500,200);
		panel = new JPanel();
		
		label0 = new JLabel("Please insert the path to your file");
		label = new JLabel();
		
		txt = new JTextField(10);
		
		button = new JButton("Open file");
		button.addActionListener(this);
		button1 = new JButton("Run");
		button1.addActionListener(this);
		
		checkboxgroupRadio = new CheckboxGroup();
	    card = new Checkbox("Card", checkboxgroupRadio, false);
	    card.addItemListener((ItemListener) this);
	    cardBreath = new Checkbox("Cardbreath", checkboxgroupRadio, false);
	    cardBreath.addItemListener((ItemListener) this);		
	
		panel.add(label0);
		panel.add(txt);
		panel.add(button);
		
		panel.add(card);
		panel.add(cardBreath);
		
		panel.add(label);
		panel.add(button1);
		
		this.add(panel);
	}
	
	public void run() throws Exception {
		CavarevEvaluation eval;
		String evaluationFile = "";
		label.setText("");
		
		if(txt.getText() != "") {
			evaluationFile = txt.getText();
			File file = new File(evaluationFile);
		    if (!file.exists()) {
		       label.setText("Your entered filepath doesn't exist. Please use another one.");
		    } else {
			    Checkbox selected = checkboxgroupRadio.getSelectedCheckbox();
			    if (selected != null) {
				    String mode = selected.getLabel();
				    
				    switch (mode) {
				    case "Card":
				    			eval = new CavarevEvaluation(evaluationFile, false);
				    			eval.run();
				    			break;
				    case "Cardbreath":	
				    			eval = new CavarevEvaluation(evaluationFile, true);
				    			eval.run();
				    			break;
				    default:	System.out.print("False input. Abort program.");
				    }		    	
			    } else {
			    	label.setText("Please choose if your file is card oder cardbreath.");
			    }		    	
		    }
		}
	}
	public static void main(String [] args){

		RunCavarevEvaluationSurface cEval = new RunCavarevEvaluationSurface();
		cEval.setVisible(true);

	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if(e.getSource() == this.button) {
			JFileChooser chooser = new JFileChooser();
			chooser.showDialog(null, "Open");

			int returnVal = 0;
			if(returnVal == JFileChooser.APPROVE_OPTION)
			{
		         txt.setText(chooser.getSelectedFile().getPath());
			}
		}
		if(e.getSource() == this.button1) {
			try {
				run();
			} catch (Exception e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}
	}

	@Override
	public void itemStateChanged(ItemEvent e) {
	}
	
}
