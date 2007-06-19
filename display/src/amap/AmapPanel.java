
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

public class AmapPanel extends JPanel 
	implements ChangeListener, PropertyChangeListener {

	AlignmentPanel alignPanel;
    JFormattedTextField textField;
    JTextField weightField;
    JSlider alignSlider;
    java.util.List<Alignment> aligns;

    public AmapPanel(java.util.List<Alignment> aligns) {
		this.aligns = aligns;
		int maxSliderVal = aligns.size()-1;

        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        //Create the label.
        JLabel sliderLabel = new JLabel("Alignment: ", JLabel.CENTER);
        sliderLabel.setAlignmentX(Component.CENTER_ALIGNMENT);


        //Create the formatted text field and its formatter.
        java.text.NumberFormat numberFormat = java.text.NumberFormat.getIntegerInstance();
        NumberFormatter formatter = new NumberFormatter(numberFormat);
        formatter.setMinimum(new Integer(0));
        formatter.setMaximum(new Integer(maxSliderVal));
        textField = new JFormattedTextField(formatter);
        textField.setValue(new Integer(0));
        textField.setColumns(5); //get some space
        textField.addPropertyChangeListener(this);

        JLabel weightLabel = new JLabel("Weight: ", JLabel.CENTER);

        weightField = new JTextField();
		weightField.setText(getWeight(0));
        weightField.setColumns(10); //get some space
        weightField.addPropertyChangeListener(this);


        //React when the user presses Enter.
        textField.getInputMap().put(KeyStroke.getKeyStroke(
                                        KeyEvent.VK_ENTER, 0),
                                        "check");
        textField.getActionMap().put("check", new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                if (!textField.isEditValid()) { //The text is invalid.
                    Toolkit.getDefaultToolkit().beep();
                    textField.selectAll();
                } else try {                    //The text is valid,
                    textField.commitEdit();     //so use it.
                } catch (java.text.ParseException exc) { }
            }
        });

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

		JScrollPane alignScroll = new JScrollPane(alignPanel);
		alignScroll.setPreferredSize(new Dimension(500,700));
		alignScroll.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		alignScroll.addComponentListener(alignPanel);

        //Create a subpanel for the label and text field.
        JPanel labelAndTextField = new JPanel(); //use FlowLayout
        labelAndTextField.add(sliderLabel);
        labelAndTextField.add(textField);
		labelAndTextField.add(weightLabel);
		labelAndTextField.add(weightField);

        //Put everything together.
        add(alignScroll);
        add(alignSlider);
        add(labelAndTextField);
        setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
    }

    /** Listen to the slider. */
    public void stateChanged(ChangeEvent e) {
        JSlider source = (JSlider)e.getSource();
        int ind = (int)source.getValue();
        updatePicture(ind); //display the next picture
		if (!source.getValueIsAdjusting()) {
			textField.setValue(new Integer(ind));
			weightField.setText(getWeight(ind));
		} else {
			textField.setText(String.valueOf(ind));
			weightField.setText(getWeight(ind));
		}
    }

    public void propertyChange(PropertyChangeEvent e) {
        if ("value".equals(e.getPropertyName())) {
            Number value = (Number)e.getNewValue();
            if (alignSlider != null && value != null) {
                alignSlider.setValue(value.intValue());
            }
        }
    }

    protected void updatePicture(int frameNum) {
		alignPanel.setIndex(frameNum);
		repaint();
    }

	protected String getWeight(int index) {
		DecimalFormat df = new DecimalFormat("0.00000");
		return df.format( aligns.get(index).getNewWeight() );
	}

	protected int calcTickSpacing(int val) {
		return ((val/10) - (val%10));
	}
}
