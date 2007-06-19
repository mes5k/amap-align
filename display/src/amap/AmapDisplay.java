
//============================================================================
// 
//  file: AmapDisplay.java
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

import java.io.*;
import java.util.*;
import javax.swing.*;
import java.awt.Dimension;

public class AmapDisplay {
	
	public static void main(String[] args) {
		
		if ( args.length != 1 ) {
			System.out.println("USAGE: java -jar amapDisplay.jar <input file>");
			System.exit(1);
		}

		try {
			FileInputStream fis = new FileInputStream(args[0]);
			AmapReader ar = new AmapReader(fis);

			List<Alignment> aligns = ar.getAlignments();

        	final JFrame frame = new JFrame("AMAP Display");
        	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        	AmapPanel panel = new AmapPanel(aligns);
        	panel.setPreferredSize( new Dimension(500,700) );

        	//Display the window.
        	frame.setContentPane(panel);
        	frame.pack();
        	frame.setVisible(true);
		} catch (Exception ioe) {
			ioe.printStackTrace();
		}

	}
}
