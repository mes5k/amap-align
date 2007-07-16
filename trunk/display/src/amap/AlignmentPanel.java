
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
import java.awt.Shape;
import java.awt.geom.Point2D;
import java.awt.font.GlyphVector;
import java.awt.font.FontRenderContext;

public class AlignmentPanel extends JPanel implements ComponentListener {

	List<Alignment> aligns;
	List<String> orderedKeys; 
	int index;
	int xOff; 
	int yOff;
	int numCharsAcross;
	int nameWidth;
	int glyphWidth;
	int len; 
	int groupOff;
	Font defaultFont;
	int height;
	int width;
	int os;
	Map<String,Color> backColMap;
	Map<String,Color> textColMap;
	public AlignmentPanel(List<Alignment> aligns, int initWidth) {
		super();
		this.aligns = aligns;
		this.width = initWidth;
		Map<String,String> seqs = aligns.get(0).getSequences();
		orderedKeys = aligns.get(0).getOrderedKeys(); 
		defaultFont = new Font("Courier",Font.BOLD,12);
		index = 0;

		// Aaaaargh
		if ( System.getProperty("os.name").matches(".*[Mm][Aa][Cc].*") )
			os = -1;
		else
			os = 1;

		backColMap = new HashMap<String,Color>();
		backColMap.put("0",new Color(255,255,255)); 
		backColMap.put("1",new Color(0,128,255)); // CMYK 1, 0.5, 0, 0
		backColMap.put("2",new Color(64,25,255)); // CMYK 0.75, 0.9, 0, 0
		backColMap.put("3",new Color(255,0,255)); // CMYK 0, 1, 0, 0

		textColMap = new HashMap<String,Color>();
		textColMap.put("0",Color.black); 
		textColMap.put("1",Color.white); 
		textColMap.put("2",new Color(255,230,41)); // CMYK 0, 0.1, 0.84, 0
		textColMap.put("3",Color.black); 

		sizeScreen( new Dimension(initWidth,0) );
	}

	private void sizeScreen(Dimension d) {
		width = (int)(d.getWidth()); 

		glyphWidth = 10;

		// find the width of the sequence name
		FontMetrics fm = getFontMetrics(defaultFont);
		nameWidth = 0;	
		for ( String key : orderedKeys )
			nameWidth = Math.max(nameWidth,fm.stringWidth(key));

		// give the seq name some breathing room
		nameWidth += glyphWidth;

		// set margins
		xOff = 20;
		yOff = 20;

		// size of one multiple alignment row
		groupOff = yOff * (orderedKeys.size() + 1);

		// figure out how many multiple alignment rows there will be	
		int seqAllowedWidth = width - xOff - nameWidth - xOff - xOff;
		numCharsAcross = seqAllowedWidth/glyphWidth;
		len = aligns.get(0).getSequences().get(orderedKeys.get(0)).length();

		height = ((len/numCharsAcross)+1) * groupOff;
		height = Math.max(height, (int)d.getHeight());

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
		FontRenderContext frc = g2.getFontRenderContext();
		Map<String,String> seqs = aligns.get(index).getSequences();
		Map<String,String> cols = aligns.get(index).getColors();
		int i = 0;
		for ( String name : orderedKeys ) {
			i++;

			String seq = seqs.get(name);
			String colString = cols.get(name);
			if ( seq == null )
				continue;

			GlyphVector gv = defaultFont.createGlyphVector(frc,seq);
			
			for ( int j = 0; j < gv.getNumGlyphs(); j++ ) {

				int xLoc = xOff+nameWidth+((j%numCharsAcross)*glyphWidth);
				int yLoc = ((j/numCharsAcross)*groupOff)+(yOff*i);

				// set glyph position
				Point2D pos = gv.getGlyphPosition(j); 
				pos.setLocation(xLoc, yLoc);
				gv.setGlyphPosition(j,pos);

				// draw the background rect 
				String colCode = colString.substring(j,j+1);
				g2.setPaint(backColMap.get( colCode ));
				g2.fillRect(xLoc-1,yLoc-yOff+6,glyphWidth,yOff);

				// draw the sequence name
				if ( j%numCharsAcross == 0 ) {
					g2.setPaint(Color.black);
					g2.drawString(name,xOff,yLoc );
				}

				// draw the glyph
				Shape glyph = gv.getGlyphOutline(j);
				g2.setPaint(textColMap.get( colCode ));
				g2.fill(glyph);
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
