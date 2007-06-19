
//============================================================================
// 
//  file: AlignmentPanel.java
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

import javax.swing.*;
import java.util.*;
import java.awt.event.*;
import java.awt.Dimension;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;

public class AlignmentPanel extends JPanel implements ComponentListener {

	List<Alignment> aligns;
	List<String> orderedKeys; 
	int index;
	int nameXOff; 
	int nameYOff;
	int len; 
	int groupOff;
	Font defaultFont;
	int height;
	int width;
	public AlignmentPanel(List<Alignment> aligns, int initWidth) {
		super();
		this.aligns = aligns;
		this.width = initWidth;
		Map<String,String> seqs = aligns.get(0).getSequences();
		orderedKeys = new ArrayList<String>( seqs.keySet() );
		Collections.sort(orderedKeys);
		defaultFont = new Font("Courier",Font.PLAIN,12);
		index = 0;
		sizeScreen( new Dimension(initWidth,0) );
	}

	private void sizeScreen(Dimension d) {
		width = (int)(d.getWidth()); 

		nameXOff = 30;
		nameYOff = 20;
		groupOff = nameYOff * (orderedKeys.size() + 1);

		int allowedWidth = width - (4*nameXOff) - nameYOff;
		Map<String,String> seqs = aligns.get(0).getSequences();

		String align = seqs.get(orderedKeys.get(0));
		FontMetrics fm = getFontMetrics(defaultFont);

		while ( fm.stringWidth(align) > allowedWidth ) 
			align = align.substring(0,align.length()-1);

		len = align.length();
		int seqLen = seqs.get(orderedKeys.get(0)).length();
		height = ((seqLen/len)+1) * groupOff;

		setPreferredSize( new Dimension(width,height) );
	}

	public void setIndex(int i) {
		if ( i >= 0 && i < aligns.size() )
			index = i;
		else
			System.out.println("invalid index: " + i);
	}	

	public void paint(Graphics g) {
		Graphics2D g2 = (Graphics2D)g;
		g2.setBackground(Color.white);
		g2.setPaint(Color.white);
		g2.fillRect(0,0,width,height);

		g2.setPaint(Color.black);
		g2.setFont( defaultFont ); 
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

		int i = 0;
		Map<String,String> seqs = aligns.get(index).getSequences();
		for ( String name : orderedKeys ) {
			i++;	
			String seq = seqs.get(name);
			if ( seq == null )
				continue;
			int begin = 0;
			int end = Math.min(len, seq.length()-1);
			String sub = seq.substring(begin,end);
			int y = 0;
			while ( sub.length() > 0 ) {
				g2.drawString(name,nameXOff,y*groupOff+nameYOff*i);
				g2.drawString(sub,nameXOff*3,y*groupOff+nameYOff*i);
				end = Math.min(end + len, seq.length()-1);
				begin = Math.min(begin + len, end ); 
				sub = seq.substring(begin,end);
				y++;
			}
		}
	}

    public void componentResized(ComponentEvent e) {
		sizeScreen( e.getComponent().getSize() );
    }

    public void componentShown(ComponentEvent e) {};
    public void componentHidden(ComponentEvent e) {};
    public void componentMoved(ComponentEvent e) {};
}
