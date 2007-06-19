
//============================================================================
// 
//  file: AmapReader.java
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

import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.*;
import java.io.*;

public class AmapReader {

	private Map<String,String> nameSeqMap;
	List<Alignment> alignments;
	InputStream is;

	double pWeight;
	double nWeight;

	public AmapReader(InputStream is) {
		this.is = is;
		nameSeqMap = new HashMap<String,String>();
		alignments = new ArrayList<Alignment>();
		read();
	}

	public List<Alignment> getAlignments() {
		return alignments;
	}

	private void read() {
		Pattern pw = Pattern.compile("^previous weight= (\\d+\\.\\d+)\\s+.*");
		Matcher pwm = pw.matcher("");

		Pattern nw = Pattern.compile("^new weight= (\\d+\\.\\d+)$");
		Matcher nwm = nw.matcher("");

		Pattern s = Pattern.compile("^(\\S+)\\s+(\\S+)$");
		Matcher sm = s.matcher("");

		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
		
			String line = null;
			while ( (line = br.readLine()) != null ) {
				pwm.reset(line);
				if ( pwm.matches() ) {
					saveAlignment();
					pWeight = new Double( pwm.group(1) ).doubleValue();
					continue;
				}

				nwm.reset(line);
				if ( nwm.matches() ) {
					nWeight = new Double( nwm.group(1) ).doubleValue();
					continue;
				}

				sm.reset(line);
				if ( sm.matches() ) {
					String key = sm.group(1);
					String val = sm.group(2);
					if ( nameSeqMap.containsKey( key ) )
						val = nameSeqMap.get(key) + val;
					nameSeqMap.put(key,val);
					continue;
				}
			}
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

	private void saveAlignment() {
		if ( pWeight > 0 && nWeight > 0 ) {
			if ( !( nameSeqMap.containsKey("wrong") ||
			        nameSeqMap.containsKey("Visited:") ||
			        nameSeqMap.containsKey("weigth:") ) ) {
				Alignment a = new Alignment( nameSeqMap, pWeight, nWeight );
				alignments.add(a);
			}
			nameSeqMap = new HashMap<String,String>();
			pWeight = -1;
			nWeight = -1;
		}
	}
}


