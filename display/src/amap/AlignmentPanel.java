
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
import java.awt.geom.Point2D;
import java.awt.font.GlyphVector;
import java.awt.font.FontRenderContext;

public class AlignmentPanel extends JPanel implements ComponentListener {

	List<Alignment> aligns;
	List<String> orderedKeys; 
	int index;
	int nameXOff; 
	int nameYOff;
	int lineWidth;
	int glyphWidth;
	int len; 
	int groupOff;
	Font defaultFont;
	int height;
	int width;
	int os;
	Map<String,Color> colMap;
	public AlignmentPanel(List<Alignment> aligns, int initWidth) {
		super();
		this.aligns = aligns;
		this.width = initWidth;
		Map<String,String> seqs = aligns.get(0).getSequences();
		orderedKeys = new ArrayList<String>( seqs.keySet() );
		Collections.sort(orderedKeys);
		defaultFont = new Font("Courier",Font.PLAIN,12);
		index = 0;
		// Aaaaargh
		if ( System.getProperty("os.name").matches(".*[Mm][Aa][Cc].*") )
			os = -1;
		else
			os = 1;
		colMap = new HashMap<String,Color>();
		colMap.put("0",new Color(255,255,255)); 
		colMap.put("1",new Color(0,128,255)); // CMYK 1, 0.5, 0, 0
		colMap.put("2",new Color(64,25,255)); // CMYK 0.75, 0.9, 0, 0
		colMap.put("3",new Color(255,0,255)); // CMYK 0, 1, 0, 0
		/*
		colMap.put("0",new Color(0,128,255)); // CMYK 1, 0.5, 0, 0
		colMap.put("1",new Color(255,230,41)); // CMYK 0, 0.1, 0.84, 0
		colMap.put("2",new Color(64,25,255)); // CMYK 0.75, 0.9, 0, 0
		colMap.put("3",new Color(255,0,255)); // CMYK 0, 1, 0, 0
		*/
		sizeScreen( new Dimension(initWidth,0) );
	}

	private void sizeScreen(Dimension d) {
		width = (int)(d.getWidth()); 

		nameXOff = 30;
		nameYOff = 20;
		glyphWidth = 12;
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

		lineWidth = (width - nameXOff*4)/glyphWidth;

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

				int xLoc = nameXOff*3+(j%lineWidth)*glyphWidth;
				int yLoc = (j/lineWidth)*groupOff+nameYOff*i;
				int glyphYLoc = os*((j/lineWidth)*groupOff-nameYOff*i);

				// set glyph position
				Point2D pos = gv.getGlyphPosition(j); 
				pos.setLocation(xLoc, glyphYLoc);
				gv.setGlyphPosition(j,pos);

				// draw the background rect 
				g2.setPaint(colMap.get( colString.substring(j,j+1) ));
				g2.fillRect(xLoc,yLoc-nameYOff+5,glyphWidth,nameYOff);

				// draw the sequence name
				if ( j%lineWidth == 0 ) {
					g2.setPaint(Color.black);
					g2.drawString(name,nameXOff,yLoc );
				}
			}

			g2.setPaint(Color.black);
			g2.drawGlyphVector(gv,0,0);

		}
	}

    public void componentResized(ComponentEvent e) {
		sizeScreen( e.getComponent().getSize() );
    }

    public void componentShown(ComponentEvent e) {};
    public void componentHidden(ComponentEvent e) {};
    public void componentMoved(ComponentEvent e) {};
}
