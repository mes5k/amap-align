
//============================================================================
// 
//  file: SaveAsFastaAction.java
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

import java.util.*;
import java.beans.*;
import java.io.*;
import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.Component;


public class SaveAsFastaAction extends AbstractAction implements PropertyChangeListener {

	Alignment a;
	Component parent;
	JFileChooser fc;
	
	SaveAsFastaAction(Component parent ) {
		super("Save as Mulitple FASTA");
		this.parent = parent;
		fc = new JFileChooser();
		PropertyChangeHandler.getPropertyChangeSupport()
		                     .addPropertyChangeListener(PropertyChangeIDs.CHANGE_ALIGNMENT.toString(), this);
		a = null;
	}

	public void actionPerformed(ActionEvent e) {
		System.out.println("saving as multiple fasta");	
		if ( a != null ) {
			try {
				int returnVal = fc.showSaveDialog(parent);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					File file = fc.getSelectedFile();
					FileWriter fw = new FileWriter(file);
					String mfa = a.toMultiFasta();
					fw.write(mfa,0,mfa.length());
					fw.close();
				}
			} catch (Exception ioe) {
				JOptionPane.showMessageDialog(parent,
				                              "Failed to Save Multiple FASTA file:\n\n" + ioe.getMessage(),
				                              "ERROR Saving File!",
				                              JOptionPane.ERROR_MESSAGE);
				ioe.printStackTrace();
			}
		} else
			System.err.println("no alignment found");
			
	}

	public void propertyChange(PropertyChangeEvent e) {
		a = (Alignment)e.getNewValue();	
	}
}
