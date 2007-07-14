
//============================================================================
// 
//  file: AmapPanel.java
// 
//  Copyright (c) 2007, Michael E. Smoot 
// 
//  This program is free software; you can redistribute it and/or modify it 
//  under the terms of the GNU General Public License as published by the 
//  Free Software Foundation; either version 2 of the License, or (at your 
//  option) any later version.
//  
//  This program is distributed in the hope that it will be useful, but 
//  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
//  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along 
//  with this program; if not, write to the Free Software Foundation, Inc., 
//  59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// 
//============================================================================

package amap;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.text.NumberFormatter;
import java.beans.*;
import java.text.DecimalFormat;
import java.util.Map;
import java.util.HashMap;

public class AmapPanel extends JPanel 
	implements ChangeListener, PropertyChangeListener, ActionListener {

	AlignmentPanel alignPanel;
    JFormattedTextField indexField;
    JFormattedTextField weightField;
    JSlider alignSlider;
    java.util.List<Alignment> aligns;
	JButton nextButton;
	JButton prevButton;
	int maxSliderVal;
	Map<JFormattedTextField,Boolean> ok2update;

    public AmapPanel(java.util.List<Alignment> aligns) {
		this.aligns = aligns;
		maxSliderVal = aligns.size()-1;

		ok2update = new HashMap<JFormattedTextField,Boolean>();

        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        //Create the label.
        JLabel sliderLabel = new JLabel("Alignment: ", JLabel.CENTER);
        sliderLabel.setAlignmentX(Component.CENTER_ALIGNMENT);

        //Create the alignment number field 
        java.text.NumberFormat intFormat = java.text.NumberFormat.getIntegerInstance();
        NumberFormatter intFormatter = new NumberFormatter(intFormat);
        intFormatter.setMinimum(new Integer(0));
        intFormatter.setMaximum(new Integer(maxSliderVal));
        indexField = new JFormattedTextField(intFormatter);
        indexField.setColumns(5); //get some space
        indexField.addPropertyChangeListener( this ); 
		handleEnterKeyStroke( indexField );
		ok2update.put(indexField,false);

        JLabel weightLabel = new JLabel("Weight: ", JLabel.CENTER);

        //Create the alignment weight field 
        java.text.NumberFormat numFormat = java.text.NumberFormat.getNumberInstance();
        NumberFormatter numFormatter = new NumberFormatter(numFormat);
        numFormatter.setMinimum(new Integer(0));
        weightField = new JFormattedTextField(numFormatter);
        weightField.setColumns(10); //get some space
        weightField.addPropertyChangeListener( this ); 
		handleEnterKeyStroke( weightField );
		ok2update.put(weightField,false);

		// create next button
		nextButton = new JButton("Next");
		nextButton.addActionListener(this);

		// create prev button
		prevButton = new JButton("Prev");
		prevButton.addActionListener(this);

        //Create the slider.
        alignSlider = new JSlider(JSlider.HORIZONTAL, 0, maxSliderVal, 0);
        alignSlider.addChangeListener(this);

        //Turn on labels at major tick marks.
        alignSlider.setMajorTickSpacing(calcTickSpacing(maxSliderVal));
        alignSlider.setPaintTicks(true);
        alignSlider.setPaintLabels(true);

        //Create the label that displays the animation.
        alignPanel = new AlignmentPanel(aligns,450);
        alignPanel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createLoweredBevelBorder(),
                BorderFactory.createEmptyBorder(10,10,10,10)));
        updatePicture(0); //display first frame

		// set up scrolling
		JScrollPane alignScroll = new JScrollPane(alignPanel);
		alignScroll.setPreferredSize(new Dimension(500,700));
		alignScroll.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		alignScroll.addComponentListener(alignPanel);

        //Create a subpanel for the label and text field.
        JPanel labelAndTextField = new JPanel(); //use FlowLayout
		labelAndTextField.add(prevButton);
		labelAndTextField.add(nextButton);
        labelAndTextField.add(sliderLabel);
        labelAndTextField.add(indexField);
		labelAndTextField.add(weightLabel);
		labelAndTextField.add(weightField);

        //Put everything together.
        add(alignScroll);
        add(alignSlider);
        add(labelAndTextField);
        setBorder(BorderFactory.createEmptyBorder(10,10,10,10));

		update(maxSliderVal);
    }

    public void stateChanged(ChangeEvent e) {
        JSlider source = (JSlider)e.getSource();
		if ( source == alignSlider )
			update((int)alignSlider.getValue());
	}

	private void update(int index) {
		updatePicture(index);
		alignSlider.setValue(index);
		indexField.setText(Integer.toString(index));
		weightField.setText(Double.toString(aligns.get(index).getWeight()));
	}

    public void propertyChange(PropertyChangeEvent e) {
       	if (!"value".equals(e.getPropertyName())) 
			return;
		//new Exception().printStackTrace();
		//System.out.println("pc weightfield: " + ok2update.get(weightField));
		//System.out.println("pc indexfield:  " + ok2update.get(indexField));
		//System.out.print("prop change: "  + e.getNewValue() + "  " );
		if ( e.getSource() == indexField && ok2update.get(indexField) ) {
		//	System.out.println("index" );
			Number value = (Number)e.getNewValue();
			if (value != null) {
				update(value.intValue());
				ok2update.put(indexField,false);
		//		System.out.println("index alignSilder val: " + alignSlider.getValue());
			}
		} else if ( e.getSource() == weightField && ok2update.get(weightField) ) {
			//System.out.println("weight" );
			Number weight = (Number)e.getNewValue();
			int index = getIndexFromWeight( weight.floatValue() );
			if (index > 0 && index <= maxSliderVal) {
				update(index);
				ok2update.put(weightField,false);
				//System.out.println("weight alignSilder val: " + alignSlider.getValue());
			}
		}
    }

	public void actionPerformed(ActionEvent e) {
		int curr = 	(int)alignSlider.getValue(); 
		//System.out.println("curr val: " + curr);
		if ( e.getSource() == nextButton ) {
			if ( curr == maxSliderVal )
				update(0);
			else
				update(curr+1);

		} else if ( e.getSource() == prevButton ) {
			if ( curr == 0 )
				update(maxSliderVal);
			else
				update(curr-1);
		}
	}

    protected void updatePicture(int frameNum) {
		alignPanel.setIndex(frameNum);
		repaint();
    }

	protected int calcTickSpacing(int val) {
		return ((val/10) - (val%10));
	}

	private int getIndexFromWeight(float w) {
		return bSearch(0,aligns.size()-1,w); 
	}

	private int bSearch(int low, int high, float weight) {
		if ( high < low )
			return -1;
		if ( low == high - 1 )
			return high;
		int mid = (high + low) / 2;
		if ( weight < aligns.get(mid).getWeight() ) 
			return bSearch(mid,high,weight);
		else if ( weight > aligns.get(mid).getWeight() ) 
			return bSearch(low,mid,weight);
		else 
			return mid;
	}

	//React when the user presses Enter.
	private void handleEnterKeyStroke( final JFormattedTextField field ) {
        field.getInputMap().put(KeyStroke.getKeyStroke( KeyEvent.VK_ENTER, 0), "check");
        field.getActionMap().put("check", new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
				//The text is invalid.
                if (!field.isEditValid()) { 
                    Toolkit.getDefaultToolkit().beep();
                    field.selectAll();
				//The text is valid, so use it.
                } else try {
					//System.out.println("commitEdit");
					ok2update.put(field,true);
					//System.out.println("weightfield: " + ok2update.get(weightField));
					//System.out.println("indexfield:  " + ok2update.get(indexField));
                    field.commitEdit();     
                } catch (java.text.ParseException exc) { }
            }
        });
	}

}
