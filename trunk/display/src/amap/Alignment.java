
//============================================================================
// 
//  file: Alignment.java
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

public class Alignment {

	Map<String,String> seqs;
	Map<String,String> colors;
	List<String> keys;
	double nWeight;

	public Alignment(List<String> keys, Map<String,String> seqs, Map<String,String> colors, double nWeight) {
		this.keys = keys;
		this.seqs = seqs;
		this.colors = colors;
		this.nWeight = nWeight;
	}

	public List<String> getOrderedKeys() {
		return keys;
	}

	public Map<String,String> getSequences() {
		return seqs;
	}

	public Map<String,String> getColors() {
		return colors;
	}

	public double getWeight() {
		return nWeight;
	}

	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append("Weight      = ");
		sb.append((new Double(nWeight)).toString());
		sb.append("\n");
		for ( String key : seqs.keySet() ) {
			sb.append(key);
			sb.append("\t");
			sb.append(seqs.get(key));
			sb.append("\n");
		}
		sb.append("\n");
		sb.append("\n");

		return sb.toString();
	}

	public String toMultiFasta() {
		StringBuffer sb = new StringBuffer();
		for ( String key : keys ) {
			sb.append(">");
			sb.append(key);
			sb.append("\n");
			sb.append(seqs.get(key));
			sb.append("\n");
		}
		return sb.toString();
	}
}
