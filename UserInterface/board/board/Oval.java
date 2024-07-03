package board;

import java.awt.*;

public class Oval {
	Point center;
	int w;
	int h;
	double angle;
	
	public Oval(Point p) {
		center = p;
	}
	public void extend(Point p) {
		w = Math.abs(p.x-center.x)*2;
		h = Math.abs(p.y-center.y)*2;
	}
	
	public void paint(Graphics g) {
		g.setColor(Color.black);
		g.fillOval(center.x-w/2, center.y-h/2, w, h);
	}
	
}
