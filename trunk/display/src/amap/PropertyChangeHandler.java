
//============================================================================
// 
//  file: PropertyChangeHandler.java
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


public class PropertyChangeHandler {

	private static PropertyChangeSupport pcs;
	private static Object o;

	static {
		o = new Object();
		pcs = new PropertyChangeSupport( o );
	}

	public static PropertyChangeSupport getPropertyChangeSupport() {
		return pcs;
	}

	public static void firePropertyChange(String id, Object oldValue, Object newValue) {
		try {
			PropertyChangeIDs.valueOf(id);
		} catch (IllegalArgumentException ex) {
			System.err.println("Illegal PropertyChangeEvent ID: " + id + "   ignoring!");
			return;
		}
		PropertyChangeEvent e = new PropertyChangeEvent(pcs, id, oldValue, newValue);
		pcs.firePropertyChange(e);
	}
}
